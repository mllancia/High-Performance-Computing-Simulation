#ifndef KERNEL_H
#define KERNEL_H
#include <stdlib.h>

/* the lattice field type */

/* it is the private matter of the backend how struct field is defined */
typedef struct field * latfield;

/* Global variables  */
/* these are defined in kernel_global.c */

/* kernel_init fills these in, but kernel_init may chose to call
 * vanilla_lattice_init, which fills these from the main program */
extern int Rx,Ry,Rs; /* rank within the torus partition, Rs = global rank */
extern int Px,Py,Ps; /* size of the torus partition, Ps = number of nodes */
extern int Gx,Gy,Gs; /* global lattice size, Gs=Gx*Gy */
extern int Lx,Ly,Ls; /* local lattice size, Ls=Lx*Ly */
extern int Ox,Oy;    /* offset of the local lattice within the global lattice */
extern int Nxp,Nxm,Nyp,Nym; /* Neighbour rank in the directions
                        x-plus, x-minus, y-plus and y-minus */

/* kernel_init fires up the device, the communication, the whole backend.
 * The lattice size is given in the arguments.
 * It fills the 'Group A' global variables directly or indirectly through
 * vanilla_lattice_init */
int kernel_master();
void kernel_init(int nx,int ny, int argc, char**argv);

/* kernel_finish() is called to give a chance for the backend to
 * chill down the device, call MPI_Finalize etc.  */
void kernel_finish();

/* generic communication tools: */
void kernel_broadcast_s(char *buf,size_t bytes); /* in-place, from rank=0 */
void kernel_broadcast_i(int buf); /* in-place, from rank=0 */
void kernel_broadcast_d(double *buf,size_t bytes); /* in-place, from rank=0 */
double kernel_globsum(double x); /* returns the sum */
void kernel_barrier();
void kernel_abort(int err);

void kernel_sendrecv(void *receive_buf,const void *send_buf,int rank_I_get_from,int rank_I_send_to,int bytes);

/* kerne_latfield_cpy copies dest <- src.
 * the dest object must be an already initialized object */
void kernel_latfield_cpy(latfield dest,latfield src);

/* I/O */

/* these are blocking collective calls for writing and reading
 * a configuration 
 * return value: 0 if everyhing went well */
int kernel_latfield_write(latfield phi,const char *file_name);
int kernel_latfield_read(latfield phi,const char *file_name);


/* constructor and destructor for the latfield object */
latfield kernel_latfield_new(void);
void kernel_latfield_delete(latfield F);

/* allocate and free in the host memory. The returned pointer
 * can be used to access memory from the host code.
 * An efficient upload and download is assumed */
void * kernel_host_alloc(size_t size);
void kernel_host_free(void *x,size_t size);

/* upload or download the local lattice to a host-accessible memory
 * that was allocated through kernel_host_alloc. It does not have
 * to work with other host allocations (like malloc or static arrays). 
 * The host local lattice format for a local field is:
 * double *host_format;
 *  host_format[y*Lx+x] = field(x,y);
 *
 * kernel_upload and kernel_download ar local operations.
 * */
void kernel_upload(latfield field,const double *host_mem);

void kernel_download(double *host_mem,latfield field);



/* classical update:
 * phi1 is communicated, then
 *
 * phi0 : = A*phi0 + B*phi1 + C*phi1*phi1*phi1 + D*phi1_neighbours + E
 */
void kernel_clupdate(latfield phi0,latfield phi1,const double A,const double B,const double C,const double D,const double E);

/* Energy. This is conserved in the absence of a 1st order time derivative */

double kernel_energy(latfield phi0,latfield phi1,const double A,const double B,const double C,const double D,const double E);

/* zeromode:
 * avg[0] = sum_sites phi(site) / volume
 * fluct[0] = sum_sites phi^2(site) / volume - avg[0]^2
 *
 * It returns the same value at all nodes, the global sums.
 * */

void kernel_zeromode(double *avg,double *fluct, latfield phi);
void kernel_test(latfield phi);
void kernel_display(latfield phi,int flag);
void kernel_synchronize();
#endif
