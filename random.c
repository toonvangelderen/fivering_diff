/* A C-program for MT19937: Real number version  (1998/4/6)        */
/*   genrand() generates one pseudorandom real number (double)     */
/* which is uniformly distributed on [0,1]-interval, for each      */
/* call. sgenrand(seed) set initial values to the working area     */
/* of 624 words. Before genrand(), sgenrand(seed) must be          */
/* called once. (seed is any 32-bit integer except for 0).         */
/* Integer generator is obtained by modifying two lines.           */
/*   Coded by Takuji Nishimura, considering the suggestions by     */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.               */
/*                                                                 */
/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */ 
/* 02111-1307  USA                                                 */
/*                                                                 */
/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */
/*                                                                 */
/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */
/*                                                                 */
/* See also 1) http://random.mat.sbg.ac.at                         */
/*          2) http://www.math.keio.ac.jp/~matumoto/emt.html       */
/*                                                                 */
/* Slightly modified by Thijs J.H. Vlugt on 21/12/1998             */


#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "path.h"


#define N                     624
#define M                     397
#define MATRIX_A              0x9908b0df   
#define UPPER_MASK            0x80000000 
#define LOWER_MASK            0x7fffffff 
#define TEMPERING_MASK_B      0x9d2c5680
#define TEMPERING_MASK_C      0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)


#define halfsqrt12 1.7320508075


static unsigned long mt[N];
static int mti;

double RandomNumber(void)
{
  int kk;
  unsigned long y;
  static unsigned long mag01[2]={0x0,MATRIX_A};
  double zzz=2.0;

  while(zzz<0.0000000000001e0||zzz>0.9999999999999e0)
  {
    if (mti>= N) 
    { 
      for (kk=0;kk<N-M;kk++) 
      {
        y=(mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
        mt[kk]=mt[kk+M]^(y>>1)^mag01[y&0x1];
      }
      for (;kk<N-1;kk++) 
      {
        y=(mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
           mt[kk]=mt[kk+(M-N)]^(y>>1)^mag01[y&0x1];
      }
      y=(mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
         mt[N-1]=mt[M-1]^(y>>1)^mag01[y & 0x1];
      mti = 0;
    }
    y=mt[mti++];
    y^=TEMPERING_SHIFT_U(y);
    y^=TEMPERING_SHIFT_S(y)&TEMPERING_MASK_B;
    y^=TEMPERING_SHIFT_T(y)&TEMPERING_MASK_C;
    y^=TEMPERING_SHIFT_L(y);
    zzz=((double)y/(unsigned long)0xffffffff);
  }
  return(zzz);
}

int Random_urandom(void){
  int rng;
  int urnd = open("/dev/urandom", O_RDONLY);
  read(urnd, &rng, sizeof(int));
  close(urnd);
  // printf("the random number is  %d\n", abs(rng));
  return abs(rng);
}

void InitializeRandomNumberGenerator(double seed)
{
  unsigned long myint;
  int kk;
  double dummy;

  myint=(unsigned long)((seed)*429496729);
  myint=myint+1000;

  if(myint%2==0) 
    myint=myint+1;
  else
    myint=myint+2;

  mt[0]=myint&0xffffffff;

  for(mti=1;mti<N;mti++)
    mt[mti]=(69069*mt[mti-1])&0xffffffff;

  for(kk=1;kk<10000;kk++)
    dummy=RandomNumber();

  return;
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK 
#undef TEMPERING_MASK_B
#undef TEMPERING_MASK_C
#undef TEMPERING_SHIFT_U
#undef TEMPERING_SHIFT_S
#undef TEMPERING_SHIFT_T
#undef TEMPERING_SHIFT_L


double ran3()
{
  static int ma[60],mj,mk,mz,i,ii,k,inext,inextp,iff,mbig,mseed;
  
  
  if (iff ==0) {
    iff=1;
    mbig = 1000000000;
    mseed = 161803398;
    mz =0;
    mj = mseed ;
    mj = mj % mbig;
    ma[55] = mj;
    mk = 1;
    for (i=1;i<55;i++){
      ii = (21*i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk<mz) mk=mk+mbig;
      mj = ma[ii];
    }
    for(k=1;k<=4;k++){
      for(i=1;i<=55;i++){
	ma[i] = ma[i] - ma[1 + ((i+30)%55)];
	if (ma[i]<mz) ma[i] = ma[i] + mbig;
      }
    }
    inext=0;
    inextp=31;
   }
  if (++inext == 56)   inext  =1;
  if (++inextp ==56)  inextp =1;
  mj = ma[inext] - ma[inextp];
  if (mj<mz) mj=mj+mbig;
  ma[inext]=mj;
  return  (double)mj / mbig;
}
 


double RandomGaussianNumber() {

    double ran1,ran2,ransq,l;

    do {
        ran1=2.0*RandomNumber()-1.0;
        ran2=2.0*RandomNumber()-1.0;
        ransq=ran1*ran1+ran2*ran2;
    } while( ransq>=1.0 );

    l=ran1*sqrt(-2.0*log(ransq)/ransq);

    return l;
}





vector RandomBrownianVector(double variance) {

    vector u;

    u.x = variance*RandomGaussianNumber();
    u.y = variance*RandomGaussianNumber();
    u.z = variance*RandomGaussianNumber();

    return u;
}


vector RandomVector(double drmax) {

    vector u;

    u.x=drmax*(RandomNumber()-0.5);
    u.y=drmax*(RandomNumber()-0.5);
    u.z=drmax*(RandomNumber()-0.5);

    return u;
}

vector RandomVector3D(vector drmax) {
    // you can give three different (box) lengths for the random vector.

    vector u;

    u.x=drmax.x*(RandomNumber()-0.5);
    u.y=drmax.y*(RandomNumber()-0.5);
    u.z=drmax.z*(RandomNumber()-0.5);

    return u;
}



vector RandomUnitVector(void) {

    vector u;
    double ran1,ran2,ranh,ransq;

    do{
        ran1 = 1.0 - 2.0*RandomNumber();
        ran2 = 1.0 - 2.0*RandomNumber();
        ransq = ran1*ran1 + ran2*ran2;
    } while( ransq >= 1.0 );

    ranh = 2.0*sqrt(1.0-ransq);
    u.x = ran1*ranh;
    u.y = ran2*ranh;
    u.z = (1.0-2.0*ransq);

    return u;
}

// Generates A Random Velocity According To A Boltzmann Distribution
double RandomVelocity(double temperature) {

    return sqrt(temperature)*RandomGaussianNumber();

}

double RandomNumberRange(double boundary_min, double boundary_plus){
    double x, length=boundary_plus-boundary_min;
    // x is the random number between [boundary_min,boundary_plus>

    x= (length*RandomNumber()+boundary_min);

    return x;
}

int RandomIntegerRange(int boundary_min, int boundary_plus){
    int x, length=boundary_plus-boundary_min;
    // x is the random number between [boundary_min,boundary_plus]

    x= ((int)floor(length*RandomNumber())+boundary_min);

    return x;
}


