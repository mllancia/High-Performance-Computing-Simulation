/* MPI + OMP backend
 *
 * Design
 * ------
 * Each rank owns an Lx*Ly block of the global lattice (x-minor, y-major).
 * The 5-point stencil needs one ring of neighbour values ("halo"), kept in
 * recv_buf[] per field:
 *     recv_buf[buf_yl]  row y=-1   (top row of the y-minus neighbour)
 *     recv_buf[buf_yr]  row y=Ly   (bottom row of the y-plus neighbour)
 *     recv_buf[buf_xl]  column x=-1 (right column of the x-minus neighbour)
 *     recv_buf[buf_xr]  column x=Lx (left column of the x-plus neighbour)
 *
 * flag_recv says whether a field's halo still matches its data. Anything that
 * writes the field (update, copy, upload) clears it; the halo is exchanged
 * only when it is needed and stale. Each kernel_clupdate therefore does
 * exactly one exchange (for the field it reads), instead of two.
 *
 * kernel_clupdate overlaps communication and computation:
 *     1. post non-blocking receives/sends for phi1's halo
 *     2. update the interior (needs no halo)       -- OpenMP
 *     3. wait for the halo
 *     4. update the outer ring of the block         -- OpenMP
 *
 * Rows in y are contiguous in memory, so they are sent straight from the
 * field (no copy). Columns in x are strided and are packed into send_buf.
 *
 * Every site is computed with the same floating-point operations in the
 * same order as kernel_cpu.c, so the field evolution is bit-for-bit
 * identical to the serial code for any process layout.
 *
 * MPI is only called by the master thread (MPI_THREAD_FUNNELED).
 */
#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <err.h>

#include "kernel.h"
#include "vanilla.h"
#include "macros.h"
#include <omp.h>

#include "mpi.h"
#include "comm.h"
#include <errno.h>
#include <math.h>

extern MPI_Comm Comm_world;
extern MPI_Comm Comm_world_x;
extern MPI_Comm Comm_world_y;
extern MPI_Status Comm_status;

enum {
    /* 2d: four boundaries l:left, r:right, the buf id is alway
     * from the point of view of the recipient */
    buf_yl=0,buf_yr=1,buf_xl=2,buf_xr=3,buf_n=4
};

/* message tags: which of my edges is being sent */
enum { tag_row_bottom=10, tag_row_top=11, tag_col_left=12, tag_col_right=13 };

typedef double real;

struct field {
    real *data; /* x-minor, y major */
    /* recv_buf is filled in from remote content */
    /* send_buf is where I write my boundaries (only the x columns are used,
     * y rows are sent directly from data) */
    real *recv_buf[buf_n];
    real *send_buf[buf_n];
    /* flag to avoid unneccessary communication */
    int   flag_recv; /* 1 if recv_buf currently matches data */
};

/* outstanding non-blocking halo exchange */
struct halo_req {
    MPI_Request r[8];
};

static double Kernel_memory_counter;

int kernel_master(){return comm_master();}


void kernel_init(int nx,int ny , int argc, char**argv)
{
    comm_init(nx, ny, argc, argv);
    if(comm_master())
        printf("# MPI layout: %d ranks as %d x %d, global %d x %d, local %d x %d, %d OpenMP threads per rank\n",
               Ps,Px,Py,Gx,Gy,Lx,Ly,omp_get_max_threads());
}

void kernel_finish() {comm_finish();}

void kernel_broadcast_s(char *buf,size_t bytes) { comm_broadcast_s(buf,bytes); }
void kernel_broadcast_i(int x) { comm_broadcast_i(x); }
void kernel_broadcast_d(double *buf,size_t bytes) { comm_broadcast_d(buf,bytes); }
double kernel_globsum(double x) { return comm_globsum(x); }
int kernel_global_or(int x) { return comm_global_or(x); }
void kernel_barrier() { comm_barrier(); }
void kernel_abort(int err){comm_abort(err);}

void kernel_latfield_cpy(latfield dest,latfield src)/*{{{*/
{
    real * restrict data_dest=dest->data;
    const real * restrict data_src=src->data;
    #pragma omp parallel for
    for(int y=0;y<Ly;y++) {
        int o=y*Lx;
        for(int x=0;x<Lx;x++) data_dest[o+x]=data_src[o+x];
    }
    /* copy the halo too if it is valid, saves one exchange */
    if(src->flag_recv) {
        for(int i=0;i<buf_n;i++)
            memcpy(dest->recv_buf[i],src->recv_buf[i],sizeof(real)*((i<2)?Lx:Ly));
    }
    dest->flag_recv=src->flag_recv;
}/*}}}*/

void *kernel_alloc(size_t a)/*{{{*/
{
    void *ptr;
#ifdef _ISOC11_SOURCE
    if(!(ptr=aligned_alloc(128,a)))
#else
    if(posix_memalign(&ptr,128,a))
#endif
        errx(1,"Allocation failure. This allocation %g MB, previous allocations %g MB\n",1e-6*a,Kernel_memory_counter);
    Kernel_memory_counter+=1e-6*a;
    return ptr;
}/*}}}*/

void kernel_free(void *ptr,size_t bytes)/*{{{*/
{
    Kernel_memory_counter-=1e-6*bytes;
    free(ptr);
}/*}}}*/

void kernel_memcpy(void *dest,const void *src,size_t bytes)/*{{{*/
{
    memcpy(dest,src,bytes);
}/*}}}*/

latfield kernel_latfield_new(void)/*{{{*/
{
    latfield F;
    NEW(F);
    F->data=kernel_alloc(sizeof(F->data[0])*Ls);
    for(int i = 0; i < buf_n; i++) {
        int size = (i < 2) ? Lx : Ly;
        F->send_buf[i] = kernel_alloc(sizeof(real) * size);
        F->recv_buf[i] = kernel_alloc(sizeof(real) * size);
        memset(F->send_buf[i],0,sizeof(real)*size);
        memset(F->recv_buf[i],0,sizeof(real)*size);
    }
    F->flag_recv = 0;
    /* first touch: the threads that will update a row also zero it */
    #pragma omp parallel for
    for(int y = 0; y < Ly; y++) {
        int o = y * Lx;
        for(int x = 0; x < Lx; x++) F->data[o + x] = 0;
    }
    return F;
}/*}}}*/

void kernel_latfield_delete(latfield F)/*{{{*/
{
    for(int i = 0; i < buf_n; i++) {
        int size = (i < 2) ? Lx : Ly;
        kernel_free(F->send_buf[i], sizeof(real) * size);
        kernel_free(F->recv_buf[i], sizeof(real) * size);
    }
    kernel_free(F->data, sizeof(F->data[0]) * Ls);
    DELETE(F);
}/*}}}*/

void * kernel_host_alloc(size_t size)/*{{{*/
{
    return kernel_alloc(size);
}/*}}}*/

void kernel_host_free(void *x,size_t size)/*{{{*/
{
    kernel_free(x,size);
}/*}}}*/

/* ---------------------------------------------------------------------- */
/* halo exchange                                                          */
/* ---------------------------------------------------------------------- */

/* post the non-blocking exchange of phi's halo */
static void halo_start(latfield phi,struct halo_req *H)/*{{{*/
{
    const real *row_bottom=&phi->data[0];
    const real *row_top=&phi->data[(Ly-1)*Lx];
    /* pack the strided x-columns */
    for(int y=0;y<Ly;y++) {
        phi->send_buf[buf_xl][y]=phi->data[y*Lx];
        phi->send_buf[buf_xr][y]=phi->data[y*Lx+Lx-1];
    }
    /* receives first, so that messages can land directly */
    MPI_Irecv(phi->recv_buf[buf_yr],Lx,MPI_DOUBLE,Nyp,tag_row_bottom,Comm_world,&H->r[0]);
    MPI_Irecv(phi->recv_buf[buf_yl],Lx,MPI_DOUBLE,Nym,tag_row_top,   Comm_world,&H->r[1]);
    MPI_Irecv(phi->recv_buf[buf_xr],Ly,MPI_DOUBLE,Nxp,tag_col_left,  Comm_world,&H->r[2]);
    MPI_Irecv(phi->recv_buf[buf_xl],Ly,MPI_DOUBLE,Nxm,tag_col_right, Comm_world,&H->r[3]);
    /* my bottom row is the top halo of the y-minus neighbour, etc. */
    MPI_Isend(row_bottom,            Lx,MPI_DOUBLE,Nym,tag_row_bottom,Comm_world,&H->r[4]);
    MPI_Isend(row_top,               Lx,MPI_DOUBLE,Nyp,tag_row_top,   Comm_world,&H->r[5]);
    MPI_Isend(phi->send_buf[buf_xl], Ly,MPI_DOUBLE,Nxm,tag_col_left,  Comm_world,&H->r[6]);
    MPI_Isend(phi->send_buf[buf_xr], Ly,MPI_DOUBLE,Nxp,tag_col_right, Comm_world,&H->r[7]);
}/*}}}*/

static void halo_finish(latfield phi,struct halo_req *H)/*{{{*/
{
    MPI_Waitall(8,H->r,MPI_STATUSES_IGNORE);
    phi->flag_recv=1;
}/*}}}*/

/* blocking version, used by the measurements */
static void halo_update(latfield phi)/*{{{*/
{
    if(phi->flag_recv) return;
    struct halo_req H;
    halo_start(phi,&H);
    halo_finish(phi,&H);
}/*}}}*/

/* kept for compatibility with older code that called it directly */
void kernel_sendrecv_boundaries(latfield phi)
{
    phi->flag_recv=0;
    halo_update(phi);
}

/* upload (host->device) and download (device->host) */
void kernel_upload(latfield field,const double *host_mem)/*{{{*/
{
    memcpy(field->data,host_mem,Ls*sizeof(double));
    field->flag_recv=0;
}/*}}}*/

void kernel_download(double *host_mem,latfield field)/*{{{*/
{
    kernel_memcpy(host_mem,field->data,Ls*sizeof(double));
}/*}}}*/

/* ---------------------------------------------------------------------- */
/* I/O: gather the sublattices on the master and write the same format    */
/* as kernel_cpu.c, so field2image works unchanged                        */
/* ---------------------------------------------------------------------- */

int kernel_latfield_write(latfield phi,const char *file_name)/*{{{*/
{
    enum { header_size=4 };
    int error=0;
    float *local=kernel_host_alloc(sizeof(float)*Ls);
    for(int i=0;i<Ls;i++) local[i]=phi->data[i];
    float *all=NULL;
    if(comm_master()) all=kernel_host_alloc(sizeof(float)*Gs);
    MPI_Gather(local,Ls,MPI_FLOAT,all,Ls,MPI_FLOAT,0,Comm_world);
    if(comm_master()) {
        float *disk=kernel_host_alloc(sizeof(float)*(Gs+header_size));
        ((int32_t*)disk)[0]=826886224;
        ((int32_t*)disk)[1]=1; /* could be Lz */
        ((int32_t*)disk)[2]=Gy;
        ((int32_t*)disk)[3]=Gx;
        /* block r sits at offset (Ox,Oy)=((r%Px)*Lx,(r/Px)*Ly) */
        for(int r=0;r<Ps;r++) {
            int ox=(r%Px)*Lx, oy=(r/Px)*Ly;
            for(int y=0;y<Ly;y++)
                for(int x=0;x<Lx;x++)
                    disk[header_size+(oy+y)*Gx+ox+x]=all[r*Ls+y*Lx+x];
        }
        FILE *F=fopen(file_name,"w");
        if(F==NULL) { warn("Could not create %s\n",file_name); error=1; }
        else if(fwrite(disk,sizeof(disk[0]),Gs+header_size,F)!=(size_t)(Gs+header_size) || fclose(F)!=0)
            error=1;
        kernel_host_free(disk,sizeof(float)*(Gs+header_size));
        kernel_host_free(all,sizeof(float)*Gs);
    }
    kernel_host_free(local,sizeof(float)*Ls);
    MPI_Bcast(&error,1,MPI_INT,0,Comm_world);
    return error;
}/*}}}*/

int kernel_latfield_read(latfield phi,const char *file_name)
{
    /* optional, not implemented */
    return 1;
}

/* ---------------------------------------------------------------------- */
/* classical update:                                                      */
/*   phi0 := A*phi0 + B*phi1 + C*phi1^3 + D*(sum of phi1 neighbours) + E  */
/* ---------------------------------------------------------------------- */

/* one site; the neighbour sum is (down + up + n1 + n2) where n1,n2 follow
 * the order kernel_cpu.c uses at the same GLOBAL site, so the result is
 * bit-identical */
#define SITE_UPDATE(r0,r1,down,up,n1,n2) \
    ((r0)=A*(r0)+(r1)*(B+C*(r1)*(r1))+D*((down)+(up)+(n1)+(n2))+E)

/* full row y: row1m/row1p point to the rows below/above (or to the halo) */
static inline void update_row(real * restrict row0,const real * restrict row1,
        const real * restrict row1m,const real * restrict row1p,
        real left,real right,
        const double A,const double B,const double C,const double D,const double E)
{
    const int Lxm1=Lx-1;
    /* the serial code adds (right, left) at global x=0 and (left, right)
     * everywhere else; floating point addition is not associative, so we
     * follow the same order to stay bit-identical */
    if(Ox==0) SITE_UPDATE(row0[0],row1[0],row1m[0],row1p[0],row1[1],left);
    else      SITE_UPDATE(row0[0],row1[0],row1m[0],row1p[0],left,row1[1]);
    for(int x=1;x<Lxm1;x++)
        SITE_UPDATE(row0[x],row1[x],row1m[x],row1p[x],row1[x-1],row1[x+1]);
    SITE_UPDATE(row0[Lxm1],row1[Lxm1],row1m[Lxm1],row1p[Lxm1],row1[Lxm1-1],right);
}

void kernel_clupdate(latfield phi0,latfield phi1,const double A,const double B,const double C,const double D,const double E)/*{{{*/
{
    const int Lxm1=Lx-1;
    const int Lym1=Ly-1;
    struct halo_req H;
    int need_halo=!phi1->flag_recv;

    /* 1. start the exchange of phi1's halo (only if it is stale) */
    if(need_halo) halo_start(phi1,&H);

    /* 2. interior: y=1..Ly-2, x=1..Lx-2, no halo needed */
    #pragma omp parallel for
    for(int y=1;y<Lym1;y++) {
        int o=y*Lx;
        real * restrict row0=&phi0->data[o];
        const real * restrict row1=&phi1->data[o];
        const real * restrict row1m=&phi1->data[o-Lx];
        const real * restrict row1p=&phi1->data[o+Lx];
        for(int x=1;x<Lxm1;x++)
            SITE_UPDATE(row0[x],row1[x],row1m[x],row1p[x],row1[x-1],row1[x+1]);
    }

    /* 3. wait for the halo */
    if(need_halo) halo_finish(phi1,&H);

    const real *hyl=phi1->recv_buf[buf_yl], *hyr=phi1->recv_buf[buf_yr];
    const real *hxl=phi1->recv_buf[buf_xl], *hxr=phi1->recv_buf[buf_xr];

    /* 4a. bottom and top rows (full rows, including corners) */
    {
        const real *r1=&phi1->data[0];
        const real *r1p=(Ly>1)?&phi1->data[Lx]:hyr;
        update_row(&phi0->data[0],r1,hyl,r1p,hxl[0],hxr[0],A,B,C,D,E);
    }
    if(Lym1>0) {
        int o=Lym1*Lx;
        const real *r1=&phi1->data[o];
        update_row(&phi0->data[o],r1,&phi1->data[o-Lx],hyr,hxl[Lym1],hxr[Lym1],A,B,C,D,E);
    }

    /* 4b. left and right columns, y=1..Ly-2 */
    #pragma omp parallel for
    for(int y=1;y<Lym1;y++) {
        int o=y*Lx;
        real * restrict row0=&phi0->data[o];
        const real * restrict row1=&phi1->data[o];
        const real * restrict row1m=&phi1->data[o-Lx];
        const real * restrict row1p=&phi1->data[o+Lx];
        if(Ox==0) SITE_UPDATE(row0[0],row1[0],row1m[0],row1p[0],row1[1],hxl[y]);
        else      SITE_UPDATE(row0[0],row1[0],row1m[0],row1p[0],hxl[y],row1[1]);
        SITE_UPDATE(row0[Lxm1],row1[Lxm1],row1m[Lxm1],row1p[Lxm1],row1[Lxm1-1],hxr[y]);
    }

    /* phi0 changed, its halo is stale; phi1's halo stays valid */
    phi0->flag_recv=0;
}/*}}}*/

/* ---------------------------------------------------------------------- */
/* Energy. Conserved in the absence of the first order time derivative.  */
/* Same formula as kernel_cpu.c, now summed over ALL local sites,         */
/* using the halo for the sites on the edge of the block.                 */
/* ---------------------------------------------------------------------- */

static inline double site_energy(double r0,double r1,double neigh,
        double a,double b,double c,double d,double e)
{
    double diff=r0-r1;
    return a*diff*diff+r1*(e+r1*(b+c*r1*r1)+d*neigh);
}

double kernel_energy(latfield phi0,latfield phi1,const double A,const double B,const double C,const double D,const double E)/*{{{*/
{
    halo_update(phi1); /* no communication if the halo is still valid */
    const int Lxm1=Lx-1;
    const double a=1./(2.*D), b=(1-A-B)*a, c=-0.5*C*a;
    const double d=-0.5, e=-2*E*a;
    const real *hyl=phi1->recv_buf[buf_yl], *hyr=phi1->recv_buf[buf_yr];
    const real *hxl=phi1->recv_buf[buf_xl], *hxr=phi1->recv_buf[buf_xr];
    double energy_local=0;

    #pragma omp parallel for reduction(+:energy_local)
    for(int y=0;y<Ly;y++) {
        int o=y*Lx;
        const real * restrict row0=&phi0->data[o];
        const real * restrict row1=&phi1->data[o];
        const real * restrict row1m=(y>0)?&phi1->data[o-Lx]:hyl;
        const real * restrict row1p=(y<Ly-1)?&phi1->data[o+Lx]:hyr;
        double sum=0;
        sum+=site_energy(row0[0],row1[0],row1m[0]+row1p[0]+row1[1]+hxl[y],a,b,c,d,e);
        for(int x=1;x<Lxm1;x++)
            sum+=site_energy(row0[x],row1[x],row1m[x]+row1p[x]+row1[x-1]+row1[x+1],a,b,c,d,e);
        sum+=site_energy(row0[Lxm1],row1[Lxm1],row1m[Lxm1]+row1p[Lxm1]+row1[Lxm1-1]+hxr[y],a,b,c,d,e);
        energy_local+=sum;
    }
    return kernel_globsum(energy_local)/Gs;
}/*}}}*/


void kernel_zeromode(double *avg,double *fluct, latfield phi)/*{{{*/
{
    double sum1l=0,sum2l=0;
    #pragma omp parallel for reduction(+:sum1l,sum2l)
    for(int y=0;y<Ly;y++) {
        int o=y*Lx;
        const real * restrict row=&phi->data[o];
        double sum1=0;
        double sum2=0;
        for(int x=0;x<Lx;x++) {
            sum1+=row[x];
            sum2+=row[x]*row[x];
        }
        sum1l+=sum1;
        sum2l+=sum2;
    }
    *avg=kernel_globsum(sum1l)/Gs;
    *fluct=kernel_globsum(sum2l)/Gs-(*avg)*(*avg);
}/*}}}*/

void kernel_synchronize(){comm_barrier();}
