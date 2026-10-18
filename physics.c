#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <err.h>
#include <sys/stat.h>

#include "kernel.h"
#include "random.h"
#include "timer.h"

/* physics parameters, you may liberate them into runtime inputs */
double Cooling=0.5;
double Coolingoff=15;
double Dt=0.05;
double Noise=1.0;
double Start_zeromode=0.0;
double Msqr=-1;
double Lambda=6.0;
double Magnetic=0.0;
double Current_cooling;
double Tmax=0.;/* maximum t you want to observe, defaulted to be the same number as Gx */
char Seed[256]="test";

/* command line inputs */
int OptNx,OptNy; /* lattice size */
int OptPx,OptPy; /* partition size */
int Test;
int Verbosity;

/* output */
char Outname[256]="bulk.dat";
char Sequence_dir[64]; /* default: no sequence output */
FILE *Out;

/* write the field into a file */
void dump_field(latfield phi,const char *file,...)/*{{{*/
{
    char fn[256];
    FILE *F;
    va_list L;
    double error=0;
    va_start(L,file);
    vsnprintf(fn,sizeof(fn),file,L);
    va_end(L);
    kernel_latfield_write(phi,fn);
}/*}}}*/

/* classical_equation(phi0,phi1) updates the field twice {{{
 *
 * phi0 is the input, phi(t-1)
 * phi1 is the input, phi(t)
 * phi0 is the output, phi(t+1)
 * phi1 is the output, phi(t+2)
 */

void classical_equation(latfield phi0,latfield phi1)
{
    const double A=(2.*Current_cooling*Dt-1.);
    const double B=2-Dt*Dt*(4.+Msqr)-2.*Current_cooling*Dt;
    const double C=-Lambda*(1./6.)*Dt*Dt;
    const double D=+Dt*Dt;
    const double E=-Dt*Dt*Magnetic;
    kernel_clupdate(phi0,phi1,A,B,C,D,E);
    kernel_clupdate(phi1,phi0,A,B,C,D,E);
}/*}}}*/


void measurements(latfield phi0,latfield phi1,double t) /*{{{*/
{
    const double A=(2.*Current_cooling*Dt-1.);
    const double B=2-Dt*Dt*(4.+Msqr)-2.*Current_cooling*Dt;
    const double C=-Lambda*(1./6.)*Dt*Dt;
    const double D=+Dt*Dt;
    const double E=-Dt*Dt*Magnetic;
    double avg,fluct;
    double energy=kernel_energy(phi0,phi1,A,B,C,D,E);
    kernel_zeromode(&avg,&fluct,phi1);
	static int counter;

	if(Sequence_dir[0]) {
		static int created;
		if((created==0)&&kernel_master()) {
			if(mkdir(Sequence_dir,0775) && errno!=EEXIST)
				err(1,"Could not create the directory %s\n",Sequence_dir);
			created=1;
		}
		kernel_broadcast_i(counter);
		dump_field(phi1,"%s/seq.%04d",Sequence_dir,counter);
	}


	if (kernel_master()){
		/* write the output */
		static char mode[4]="w";
		if(Out==NULL && (Out=fopen(Outname,mode))==NULL)
			errx(1,"Could not create the output file %s",Outname);
		fprintf(Out,"%g %g %g %g\n",t,avg,fluct,energy);

		if(++counter%100==0) { fclose(Out); Out=NULL; sprintf(mode,"a"); }
	}


}/*}}}*/

void cmdline(int argc,char **argv)/*{{{*/
{
    int o;
    const char *Help="\n"
        "Opions: -N <N>        use NxN global lattice \n"
        "        -N <Nx>x<Ny>  global lattice size\n"
        "        -S <dir>      write sequence images from the simulation\n"
        "        -t <test>     0 (default) domain wall simulation \n"
        "                      1 benchmark\n"
//        "        -P <py>       or\n"
//        "        -P <px>x<py>  use px * py nodes in this layout\n"
        ;
    while((o=getopt(argc,argv,"N:P:S:t:hv"))!=-1) switch(o) {
        case 'N': if(sscanf(optarg,"%dx%d",&OptNx,&OptNy)!=2 &&
                          ((OptNx=-1)) && sscanf(optarg,"%d",&OptNy)!=1)
                      errx(1,"Invalid use of the -N option");
                  break;
        case 'P': if(sscanf(optarg,"%dx%d",&OptPx,&OptPy)!=2 &&
                          ((OptPx=1)) && sscanf(optarg,"%d",&OptPy)!=1)
                      errx(1,"Invalid use of the -P option");
                  break;
        case 'S': if(strncpy(Sequence_dir,optarg,sizeof(Sequence_dir))
                          [sizeof(Sequence_dir)-1]) 
                      errx(1,"Sequence directory was too long\n");
                      break;
        case 't': Test=atoi(optarg); break;
        case 'h': puts(Help); break;
        case 'v': Verbosity++; break;
        default: errx(1,"You can ask for help -h\n");
    }
    if(OptNx==-1) OptNx=OptNy;
    if(Seed[0]==0) random_bytes(Seed,sizeof(Seed));
    if(Verbosity>0) printf("Random Seed=%s\n",Seed);
}/*}}}*/

void domain_wall_simulation()/*{{{*/
{
    latfield phi0,phi1;
    phi0=kernel_latfield_new();
    phi1=kernel_latfield_new();
    random_field(phi0,Start_zeromode,Noise);
    kernel_latfield_cpy(phi1,phi0);
    if(Tmax==0) Tmax=Gx;
	int timax=round(Tmax/Dt);
    int meas_freq=round(1./Dt);
    for(int ti=0;ti<timax;ti+=2 /*double steps */) {
        double t=ti*Dt;
        Current_cooling=(t>0 && t<Coolingoff)?Cooling:0;
//        printf("step %d\n",ti);
        if(ti%meas_freq==0) {
            fprintf(stderr,"\r %d         \r",(int)(t+0.5));
            measurements(phi0,phi1,t);
        }
        classical_equation(phi0,phi1);
    }
    kernel_latfield_delete(phi1);
    kernel_latfield_delete(phi0);
}/*}}}*/

void benchmark()
{
    timer T=TIMER_INITIALIZER;
    latfield phi0,phi1;
    phi0=kernel_latfield_new();
    phi1=kernel_latfield_new();
    random_field(phi0,Start_zeromode,Noise);
    kernel_latfield_cpy(phi1,phi0);
    Current_cooling=0;
    int timax=100;
    const int flop=13; 
    const int mem1=sizeof(double);
    const int mem=3*mem1;
    int ti=0;

    classical_equation(phi0,phi1);
    classical_equation(phi0,phi1);

    kernel_synchronize();
    timer_on(&T);
    for(ti=0;ti<timax;ti+=10)  {
        classical_equation(phi0,phi1);
        classical_equation(phi0,phi1);
        classical_equation(phi0,phi1);
        classical_equation(phi0,phi1);
        classical_equation(phi0,phi1);
    }
    kernel_synchronize();
    timer_off(&T);
    if (kernel_master()){
	    //printf("Gs=%d Ls=%d\n",Gs,Ls);
	    double gflop_total=1e-9*Gs*flop*ti;
	    double gflop_local=1e-9*Ls*flop*ti;
	    double gbyte_local=1e-9*Ls*mem*ti;
	    double gbyte_comm=1e-9*((Px>1?2*Ly*mem1:0)+(Py>1?2*Lx*mem1:0))*ti;
	    double sec=timer_tot(&T);
	    printf("Runtime = %.3f iterations = %d\n",sec,ti);
	    printf("Global performance:\t%.1f Gflop/s\n",gflop_total/sec);
	    printf("Local performance:\t%.1f Gflop/s\n",gflop_local/sec);
	    printf("Local effective BW:\t%.1f Gbyte/s\n",gbyte_local/sec);
	    if(Ps>1)
		printf("Local outgoing BW:\t%.1f Gbyte/s\n",gbyte_comm/sec);
    }
    kernel_latfield_delete(phi1);
    kernel_latfield_delete(phi0);
}

int main(int argc,char **argv)/*{{{*/
{
    cmdline(argc,argv);
    kernel_init(OptNx,OptNy, argc, argv);
    random_init(Seed);
    switch(Test) {
        case 1: benchmark(); break;
        default: domain_wall_simulation(); break;
    }
    kernel_finish();
}/*}}}*/
