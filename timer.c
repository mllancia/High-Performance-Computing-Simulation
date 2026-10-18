#include "timer.h"
#include <stdlib.h>

void timer_on(timer *T)/*{{{*/
{
    T->n++;
    gettimeofday(&T->start,NULL);
}/*}}}*/

void timer_off(timer *T)/*{{{*/
{
    gettimeofday(&T->stop,NULL);
    T->last_usec=1000000ULL*(T->stop.tv_sec-T->start.tv_sec)
        +(T->stop.tv_usec-T->start.tv_usec);
    T->total_usec+=T->last_usec;
}/*}}}*/

void timer_nano_add(timer *T,uint64_t nsec) { T->total_nsec+=nsec; }

void timer_reset(timer *T) { T->total_usec=0; T->total_nsec=0; }

double timer_tot(timer *T) { 
    T->total_usec+=T->total_nsec/1000; 
    T->total_nsec=0; 
    return T->total_usec*1e-6;
}

double timer_avg(timer *T) 
{
    timer_tot(T);
    return T->total_usec*1e-6/T->n; 
}

int timer_n(timer *T) { return T->n; }

double timer_sec(timer *T) { return T->last_usec*1e-6; }
