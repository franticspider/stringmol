/* A C-program for MT19937: Real number version([0,1)-interval) (1998/4/6) */
/*   genrand() generates one pseudorandom real number (double) */
/* which is uniformly distributed on [0,1)-interval, for each  */
/* call. sgenrand(seed) set initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be      */
/* called once. (seed is any 32-bit integer except for 0).     */
/* Integer generator is obtained by modifying two lines.       */
/*   Coded by Takuji Nishimura, considering the suggestions by */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.           */

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

/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */

#include<stdio.h>

#include "mt19937-2.h"

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializing the array with a NONZERO seed */
void
sgenrand(unsigned long seed)
//    unsigned long seed;
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}


double /* generating reals */
/* unsigned long */ /* for integer generation */
genrand()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand(4357); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }

    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return ( (double)y * 2.3283064365386963e-10 ); /* reals: [0,1)-interval */
    /* return y;   for integer generation */
}



//double /* generating reals */
/* unsigned long */ /* for integer generation */
unsigned long genrandint()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand(4357); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }

    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    //*return ( (double)y * 2.3283064365386963e-10 ); /* reals: [0,1)-interval */
     return y;  /* for integer generation */
}

/*UTILITY FUNCTION FOR RE-SEEDING ON RESTART*/
int mt_get_mti(){
	return mti;
}

void mt_set_mti(int val){
	mti = val;
}


/*RECORDING FUNCTION*/
void print_mt(FILE *fp){
	int i=0;
	fprintf(fp,"MTI %d\n",mti);
	for(i=0;i<N;i++){
		fprintf(fp,"%lu\n",mt[i]);
	}
}



int report_load_error(enum load_mt_errcode ec, int posn, const char *fn){
	char msg[128];
	switch(ec){
	case load_mt_nofile:
		sprintf(msg,"ERROR: can't open mt file %s\n",fn);
		break;
	case load_mt_bad_mti:
		sprintf(msg,"ERROR: bad mti value in file %s\n",fn);
		break;
	case load_mt_bad_mt:
		sprintf(msg,"ERROR: bad mt value at posn %d in file %s\n",posn,fn);
		break;
	default:
		sprintf(msg,"UNKNOWN ERROR LOADING MT RNG\n");
		break;
	}

	if(ec != load_mt_success){
		printf("%s",msg);
		fflush(stdout);
	}

	return ec;
}



/*LOADING FUNCTION */
int load_mt(const char *fn){

	FILE *fp;
	const int maxl = 128;
	char line[maxl];
	int mtival,args_read=0,ii;
	enum load_mt_errcode errcode = load_mt_success;
	unsigned long mtval;

	if((fp = fopen(fn,"r"))!=NULL){

		//First line gets the position of mti
		if((fgets(line,maxl,fp))!=NULL){
			args_read = sscanf(line,"MTI %d",&mtival);
			if(args_read == 1){
				mti = mtival;
				//Now we can read the rest of the data;
				for(ii=0;ii<N;ii++){
					if((fgets(line,maxl,fp))!=NULL){

						args_read = sscanf(line,"%lu",&mtval);
						if(args_read == 1){
							mt[ii] = mtval;
						}
						else{
							errcode = load_mt_bad_mt;
							break;
						}
					}
					else{
						errcode = load_mt_bad_mt;
						break;
					}
				}
			}
			else{
				errcode = load_mt_bad_mti;
			}
		}
		else{
			errcode = load_mt_bad_mti;
		}

		fclose(fp);

	}
	else{

		errcode = load_mt_nofile;
	}

	return report_load_error(errcode,ii,fn);
}



// this main() outputs first 1000 generated numbers
/*
main()
{
    int j;

    sgenrand(4357); //any nonzero integer can be used as a seed
    for (j=0; j<1000; j++) {
        printf("%10.8f ", genrand());
        if (j%8==7) printf("\n");
    }
    printf("\n");
}
*/
