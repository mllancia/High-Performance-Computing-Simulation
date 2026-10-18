/* the definitions of the global variables that
 * are shared between the backend and the frontend */
#include "kernel.h"


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




int Rx,Ry,Rs; /* rank within the torus partition, Rs = global rank */
int Px,Py,Ps; /* size of the torus partition, Ps = number of nodes */
int Gx,Gy,Gs; /* global lattice size, Gs=Gx*Gy */
int Lx,Ly,Ls; /* local lattice size, Ls=Lx*Ly */
int Ox,Oy;    /* offset of the local lattice within the global lattice */
int Nxp,Nxm,Nyp,Nym; /* Neighbour rank in the directions
                        x-plus, x-minus, y-plus and y-minus */

