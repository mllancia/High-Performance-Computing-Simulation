/* No communication whatsoever, CPU only,can add openMP */
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


typedef double real;

struct field {
    real *data; /* x-minor, y major */
};

static double Kernel_memory_counter;

int kernel_master(){return 1;}

void kernel_init(int nx,int ny, int argc, char**argv)
{
    vanilla_lattice_init(nx,ny,0,1);
}

void kernel_finish() {}

void kernel_broadcast_s(char *buf,size_t bytes) {}
void kernel_broadcast_i(int x) {}
void kernel_broadcast_d(double *buf,size_t bytes) {}
double kernel_globsum(double x) { return x; }
int kernel_global_or(int x) { return x; }
void kernel_barrier() {}

void kernel_latfield_cpy(latfield dest,latfield src)/*{{{*/
{ /////////////////////////// Copy data between fields
    /* we have to apply the first touch principle */
    real * restrict data_dest=dest->data;
    const real * restrict data_src=src->data;
    #pragma omp parallel for  ////////////// Added omp pragma, Divide the loop into threads, each thread copies part of the lines (y)
    for(int y=0;y<Ly;y++) {
        int o=y*Lx;
        for(int x=0;x<Lx;x++) data_dest[o+x]=data_src[o+x];
    }
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
{ /////////////////////// Create a new field and initialize it with zeros.
    latfield F;
    NEW(F);
    F->data=kernel_alloc(sizeof(F->data[0])*Ls);
    /* first touch */
    #pragma omp parallel for  ////////////// Added omp pragma, Threads initialize different lines simultaneously (first-touch for NUMA).
    for(int y=0;y<Ly;y++) {
        int o=y*Lx;
        for(int x=0;x<Lx;x++) F->data[o+x]=0;
    }
    return F;
}/*}}}*/

void kernel_latfield_delete(latfield F)/*{{{*/
{
    kernel_free(F->data,sizeof(F->data[0])*Ls);
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

/* upload (host->device) and download (device->host) */
void kernel_upload(latfield field,const double *host_mem)/*{{{*/
{
    memcpy(field->data,host_mem,Ls*sizeof(double));
}/*}}}*/

void kernel_download(double *host_mem,latfield field)/*{{{*/
{
    memcpy(host_mem,field->data,Ls*sizeof(double));
}/*}}}*/

int kernel_latfield_write(latfield phi,const char *file_name)/*{{{*/
{
    FILE *F=fopen(file_name,"w");
    int error=0;
    if(F==NULL) { warn("Could not create %s\n",file_name); return error=1; }
    enum { header_size=4 };
    int len=Ls+header_size;
    float *disk;
    disk=kernel_host_alloc(sizeof(disk[0])*len);
    ((int32_t*)disk)[0]=826886224;
    ((int32_t*)disk)[1]=1; /* could be Lz */
    ((int32_t*)disk)[2]=Ly;
    ((int32_t*)disk)[3]=Lx;
    for(int y=0;y<Ly;y++)
        for(int x=0;x<Lx;x++) disk[header_size+y*Lx+x]=phi->data[y*Lx+x];
    if(fwrite(disk,sizeof(disk[0]),len,F)!=len || fclose(F)!=0)
        error=1;
    return error;
}/*}}}*/



int kernel_latfield_read(latfield phi,const char *file_name);

/* classical update:
 * phi1 is communicated, then
 *
 * phi0 : = A*phi0 + B*phi1 + C*phi1*phi1*phi1 + D*phi1_neighbours + E
 */

void kernel_clupdate(latfield phi0,latfield phi1,const double A,const double B,const double C,const double D,const double E)/*{{{*/
{ //////////////////////////// Classic update of the equation \phi^4
    /* this version uses no buffers but accesses the boundaries directly */
    const int Lxm1=Lx-1;
    #pragma omp parallel for    /////////////// Added omp pragma, Each thread calculates the time evolution of different rows in the field.
    for(int y=0;y<Ly;y++) {
        int o=y*Lx;
        real * restrict row0=&phi0->data[o];
        const real * restrict row1=&phi1->data[o];
        const real * restrict row1m=&phi1->data[(y?o:Ls)-Lx];
        const real * restrict row1p=&phi1->data[(o+Lx<Ls)?o+Lx:0];
        row0[0]=A*row0[0]+row1[0]*(B+C*row1[0]*row1[0])
            +D*(row1m[0]+row1p[0]+row1[1]+row1[Lxm1])+E;
        for(int x=1;x<Lxm1;x++) 
            row0[x]=A*row0[x]+row1[x]*(B+C*row1[x]*row1[x])
                +D*(row1m[x]+row1p[x]+row1[x-1]+row1[x+1])+E;
        row0[Lxm1]=A*row0[Lxm1]+row1[Lxm1]*(B+C*row1[Lxm1]*row1[Lxm1])
            +D*(row1m[Lxm1]+row1p[Lxm1]+row1[Lxm1-1]+row1[0])+E;
    }
}/*}}}*/

/* Energy. This is conserved in the absence of a 1st order time derivative */

double kernel_energy(latfield phi0,latfield phi1,const double A,const double B,const double C,const double D,const double E)/*{{{*/
{ ///////////////////////////// Calculate the total energy of the system
    /* this version uses no buffers but accesses the boundaries directly */
    const int Lxm1=Lx-1;
    double a=1./(2.*D);
    double b=(1-A-B)*a;
    double c=-0.5*C*a;
    double d=-0.5;
    double e=-2*E*a;
    double energy=0;
    #pragma omp parallel for reduction(+:energy) //////////////// Added omp pragma, Threads calculate energy from different parts and sum the results at the end.
    for(int y=0;y<Ly;y++) {
        int o=y*Lx;
        real * restrict row0=&phi0->data[o];
        const real * restrict row1=&phi1->data[o];
        const real * restrict row1m=&phi1->data[(y?o:Ls)-Lx];
        const real * restrict row1p=&phi1->data[(o+Lx<Ls)?o+Lx:0];
        double sum=0.;
        { /* x=0 */
            double diff=row0[0]-row1[0];
            double neigh=row1m[0]+row1p[0]+row1[1]+row1[Lxm1];
            double phi=row1[0];
            sum+=a*diff*diff+phi*(e+phi*(b+c*phi*phi)+d*neigh);
        }
        for(int x=1;x<Lxm1;x++) { 
            double diff=row0[x]-row1[x];
            double neigh=row1m[x]+row1p[x]+row1[x-1]+row1[x+1];
            double phi=row1[x];
            sum+=a*diff*diff+phi*(e+phi*(b+c*phi*phi)+d*neigh);
        }
        { /* x=Lx-1 */
            double diff=row0[Lxm1]-row1[Lxm1];
            double neigh=row1m[Lxm1]+row1p[Lxm1]+row1[Lxm1-1]+row1[0];
            double phi=row1[Lxm1];
            sum+=a*diff*diff+phi*(e+phi*(b+c*phi*phi)+d*neigh);
        }
        energy+=sum;
    }
    return kernel_globsum(energy)/Gs;
}/*}}}*/


void kernel_zeromode(double *avg,double *fluct, latfield phi)/*{{{*/
{ ////////////////////////// Calculates the average and fluctuation of the field.
    double sum1l=0,sum2l=0;
    #pragma omp parallel for reduction(+:sum1l,sum2l) //////////////// Added omp pragma, Threads calculate partial sums and combine them at the end.
    for(int y=0;y<Ly;y++) {
        int o=y*Lx;
        real * restrict row=&phi->data[o];
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

void kernel_display(latfield phi, int flag){
	for (int y=0;y<Gy;y++){
		for (int x=0;x<Gx;x++){
			(flag==0)?	printf("%.6e ", phi->data[ (y%Ly)*Lx+(x%Lx)]):printf("phi(%d,%d)=%.6e\n", x,y,phi->data[ (y%Ly)*Lx+(x%Lx)]);
		}
		printf("\n");
	}
}
void kernel_synchronize(){}
