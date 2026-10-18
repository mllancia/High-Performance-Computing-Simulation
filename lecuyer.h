#ifndef LECUYER_H
#define LECUYER_H

#include <stdint.h>

struct prng
{
    uint32_t s[6];
};

double prng_next(struct prng *R);
unsigned prng_next_u(struct prng *R);
void prng_leap(struct prng *R);
void prng_init(struct prng *R,unsigned seed);
void prng_init_s(struct prng *R,char seedstring[]);
void prng_init_u(struct prng *R,unsigned seed);

#endif
