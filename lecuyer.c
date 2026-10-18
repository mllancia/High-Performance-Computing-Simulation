/* L'Ecuyer 's random number generator
 * http://www.iro.umontreal.ca/~lecuyer/myftp/papers/streams00.pdf
 */

#include "lecuyer.h"

#define prng_m1 4294967087UL
#define prng_m2 4294944443UL


double prng_next(struct prng *R)/*{{{*/
{
    /* returns a random number between 0 and 1, advances the prng */
    uint64_t xnew,ynew;
    uint32_t z;
    xnew=(prng_m1+1403580UL*(uint64_t)R->s[1]- 810728UL*(uint64_t)R->s[2]) % prng_m1;
    ynew=(prng_m2+ 527612UL*(uint64_t)R->s[3]-1370589UL*(uint64_t)R->s[5]) % prng_m2;
    R->s[2]=R->s[1];R->s[1]=R->s[0];R->s[0]=(uint32_t)xnew;
    R->s[5]=R->s[4];R->s[4]=R->s[3];R->s[3]=(uint32_t)ynew;
    z=(xnew-ynew+prng_m1)%prng_m1;
    return z*2.32830654929573e-10;
}/*}}}*/

unsigned prng_next_u(struct prng *R)/*{{{*/
{
    /* returns a random number between 0 and 1, advances the prng */
    uint64_t xnew,ynew;
    uint32_t z;
    xnew=(prng_m1+1403580UL*(uint64_t)R->s[1]- 810728UL*(uint64_t)R->s[2]) % prng_m1;
    ynew=(prng_m2+ 527612UL*(uint64_t)R->s[3]-1370589UL*(uint64_t)R->s[5]) % prng_m2;
    R->s[2]=R->s[1];R->s[1]=R->s[0];R->s[0]=(uint32_t)xnew;
    R->s[5]=R->s[4];R->s[4]=R->s[3];R->s[3]=(uint32_t)ynew;
    z=(xnew-ynew+prng_m1)%prng_m1;
    return z;
}/*}}}*/

static void prng_mulmod(uint32_t A[3][3],uint32_t v[3],uint32_t m)/*{{{*/
{
    uint64_t r[3];
    int a;
    for(a=0;a<3;a++) {
        r[a]=(
                +(uint64_t)A[a][0]*(uint64_t)v[0]
                +(uint64_t)A[a][1]*(uint64_t)v[1]
                +(uint64_t)A[a][2]*(uint64_t)v[2]
             ) % m;
    }
    v[0]=(uint32_t)r[0];
    v[1]=(uint32_t)r[1];
    v[2]=(uint32_t)r[2];
}/*}}}*/

static uint32_t A1p127[3][3] = {
          {    2427906178UL, 3580155704UL,  949770784UL }, 
          {     226153695UL, 1230515664UL, 3580155704UL },
          {    1988835001UL,  986791581UL, 1230515664UL }
          };

static uint32_t A2p127[3][3] = {
          {    1464411153UL,  277697599UL, 1610723613UL },
          {      32183930UL, 1464411153UL, 1022607788UL },
          {    2824425944UL,   32183930UL, 2093834863UL }
          };

void prng_leap(struct prng *R)
{
    prng_mulmod(A1p127,&R->s[0],prng_m1);
    prng_mulmod(A2p127,&R->s[3],prng_m2);
}

static unsigned prng_lc(unsigned old,int i)
{
    return (1812433253UL * (old ^ (old >> 30)) + i);
}

void prng_init_u(struct prng *R,unsigned seed)
{
    int i;
    R->s[0]=prng_lc(seed,1);
    for(i=1;i<6;i++) R->s[i]=prng_lc(R->s[i-1],i+1) % prng_m2;
}

void prng_init_s(struct prng *R,char seedstring[])
{
    int i;
    /* default seed (in case of an empty string) */
    R->s[0]=17;
    for(i=1;i<6;i++) R->s[i]=R->s[i-1]*17;
    for(i=0;seedstring[i];i++) {
        int j=i%6;
        R->s[j]=((uint64_t)R->s[j]*127+(unsigned char)seedstring[i]) % prng_m2;
    }
}
