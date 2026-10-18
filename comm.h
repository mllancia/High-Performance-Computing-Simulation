#ifndef COMM_H
#define COMM_H

#include "macros.h"

/*

     -  ^  ^    *----*----*----*...............*----*----*----*
     |  |  |    |    |    |    |               |    |    |    |
  Ry=Py    Ly   *----*----*----*...............*----*----*----*
     |     |    |    |    |    |               |    |    |    |
     -     v    *----*----*----*...............*----*----*----*
                ...............................................
     .     .    ...............................................
     .  Gy .    ...............................................
     .     .    ...............................................
                ...............................................
     -     ^    *----*----*----*...............*----*----*----*
     |     |    |    |    |    |               |    |    |    |
  Ry=0     Ly   *----*----*----*...............*----*----*----*
     |  |  |    |    |    |    |               |    |    |    |
     -  v  v    *----*----*----*...............*----*----*----*
                <--    Lx   -->     ...        <--    Lx    -->
                <--                 Gx                      -->

                |----  Rx=0  --|     ...       |----  Rx=Px --| 
topology:
There are Ps processes arranged in Px times Py as torus (i.e. boundaries are periodic)
My rank is Rs, located in position Rx, Ry and I have neighbors with ranks Nxp,Nxm,Nyp,Nym

lattice splitting:
The global lattice of size Gx times Gy is splitted into local lattices of sizes Lx times Ly with offset coordinates Ox,Oy

*/
extern int Rx,Ry,Rs; /* rank within the torus partition, Rs = global rank */
extern int Px,Py,Ps; /* size of the torus partition, Ps = number of nodes */
extern int Gx,Gy,Gs; /* global lattice size, Gs=Gx*Gy */
extern int Lx,Ly,Ls; /* local lattice size, Ls=Lx*Ly */
extern int Ox,Oy;    /* offset of the local lattice within the global lattice */
extern int Nxp,Nxm,Nyp,Nym; /* Neighbour rank in the directions
                        x-plus, x-minus, y-plus and y-minus */



/* initialize MPI, assign rank id, split the lattice accordingly, i.e. determine Lx,Ly,Gx,Gy etc. do not call comm routines before this */
void comm_init(int px,int py, int argc, char**argv); /* px,py are the partition size. 0 val: guess */
void comm_finish(); /* destructor, no comm call afterwards */

void comm_barrier(void);     /* all PEs take part in the barrier */
int comm_broadcast_i(int);     /* broadcast one int */
int comm_broadcast_s(char *buf, size_t bytes);     /* broadcast string */
int comm_broadcast_d(double *buf, size_t bytes);     /* broadcast doubles */
double comm_globsum(double); /* global sum one double */
int comm_global_or(int); /* global or one int */
int comm_sendrecv(void *send_buf, size_t send_bytes, int dest, void *recv_buf, size_t recv_bytes, int src, int tag);/* send out data and receive data simultaneously */
int comm_send(void *send_buf, size_t send_bytes, int dest, int tag);/* send out data */
int comm_recv(void *recv_buf, size_t recv_bytes, int src, int tag);/* receive data */
void comm_abort(int err);/* abort */
int comm_master();/* returns whether I am in the master rank */
int comm_master_rank();/* tells which rank is master */
int comm_mine(int coords[2]);/* tells whether the coordinate is in my local lattice */
#endif
