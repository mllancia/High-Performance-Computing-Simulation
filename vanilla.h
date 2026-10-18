#ifndef VANILLA_H
#define VANILLA_H

/* vanilla_init finds the optimal layout for a given lattice and partition.
 * The global variables Rx,Ry,Rs,Px,Py,Ps,Gx,Gy,Gs,Lx,Ly,Ls,Ox,Oy are set. */
/* The kernel_init function is encouraged to call this */
void vanilla_lattice_init(int nx,int ny,int comm_rank,int comm_size);


#endif
