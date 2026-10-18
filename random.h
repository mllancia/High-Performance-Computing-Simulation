#ifndef RANDOM_H
#define RANDOM_H

/* Call random_init before you call any other random calls.
 * Do kernel_init before this. If 'seed' is not the same on all
 * nodes, all nodes use the seed parameter of the master node.
 */
void random_init(const char *seed); 

/* this function creates a random strings (different on each node),
 * but it may be used in random_init, since that uses the master seed only */
void random_bytes(unsigned char *bytes,int n);

double random_number(int y); /* uniform distrib on [0,1] , y=0..Ly-1*/
_Complex double random_complex_gaussian(int y); /* <imag^2>=<real^2>=0.5 */

/* creates gaussian white noise 
 * for all x: F(x) = shift + amplitude * normal();
 */
void random_field(latfield F,double shift,double amplitude);

#endif
