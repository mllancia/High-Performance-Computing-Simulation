#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <err.h>
#include "vanilla.h"
#include "kernel.h"
#include "macros.h"

int vanilla_lattice_validate(int nx,int ny)/*{{{*/
{
    if(nx<2 || ny<2 ) return -1;
    if(nx%2!=0 || ny%2!=0) return -1;
    return 0;
}/*}}}*/

void vanilla_lattice_init(int nx,int ny,int comm_rank,int comm_size)/*{{{*/
{
    double best_ratio=0;
    int best_px=0;
    if(vanilla_lattice_validate(nx,ny))
        errx(1,"Invalid lattice size: %dx%d\n",nx,ny);
    Gx=nx;
    Gy=ny;
    /* small initial: temporary variable, capital initial: global variable */
    /* guess optimal partition size */
    for(int px=1;px<=comm_size;px++) if(comm_size%px==0) {
        int py=comm_size/px;
        int lx=Gx/px;
        int ly=Gy/py;
        if(lx*px==Gx && ly*py==Gy && lx>=2 && ly>=2) {
            /* Ok, this px,py pair is not per se forbidden. Is it efficient? */
            double ratio=0.5*lx*ly/(lx+ly);
            if(best_ratio<ratio) { best_ratio=ratio; best_px=px; }
        }
    }
    if((Px=best_px)==0 || Px*(Py=comm_size/Px)!=comm_size)
        errx(1,"Could not lay out the lattice (%dx%d) on %d PEs\n",
                Gx,Gy,comm_size);
    Lx=Gx/Px;
    Ly=Gy/Py;
    Ls=Lx*Ly;
    Gs=Gx*Gy;
    /* virtual topology */
    Ps=comm_size;
    Rx=comm_rank%Px;
    Ry=comm_rank/Px;
    Rs=comm_rank;
    Ox=Rx*Lx;
    Oy=Ry*Ly;
}/*}}}*/

#if 0
int vanilla_latfield_write(latfield phi,const char *file_name)
{
    char file[256];
    strncpy(file,file_name,sizeof(file));
    if(file[sizeof(file)-1]) {
        file[sizeof(file)-1]=0;
        errx(1,"The file name was too long: %s...\n",file);
    }
}

int vanilla_latfield_read(latfield phi,const char *file_name);
#endif
