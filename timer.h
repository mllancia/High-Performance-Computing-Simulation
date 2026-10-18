#ifndef TIMER_H
#define TIMER_H

#include <stdint.h>
#include <sys/time.h>

struct timer {
    struct timeval start;
    struct timeval stop;
    uint64_t total_usec;
    uint64_t total_nsec;
    uint64_t last_usec;
    int n;
};

#define TIMER_INITIALIZER { .total_usec=0,.last_usec=0,.n=0,.total_nsec=0}
typedef struct timer timer;

void timer_on(timer *T);
void timer_off(timer *T);
double timer_sec(timer *T);
void timer_reset(timer *T);
void timer_nano_add(timer *T,uint64_t nsec);
double timer_tot(timer *T);
double timer_avg(timer *T);
int timer_n(timer *T);
double timer_sec(timer *T);

#endif
