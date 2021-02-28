/* Copyright (C) 2009-2012 Simon Hickinbotham                           */
/* When you use this, send an email to: sjh436@gmail.com                */
/* with an appropriate reference to your work.                          */

/* This file is part of STRINGMOL										*/

/* STRINGMOL is free software: you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* This program is distributed in the hope that it will be useful,      */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License    */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>

#include "randutil.h"

#define USING_MT 			/* This means we are using the Mersenne twister      */
#define USING_SEED_DEVRAND 	/* This means we are using dev/random to get a seed. */

#ifdef USING_MT
#include "mt19937-2.h"
#endif

/*
 * Obtain a seed from /dev/random - better than using clock, especially for array jobs
 * NB: this will only work if /dev/random is set up!
 */
int devrandomseed(){

	int randomData = open("/dev/random", O_RDONLY);
	int sjhRandomInteger;
  ssize_t rs;
	rs=read(randomData, &sjhRandomInteger, sizeof sjhRandomInteger);
  if(rs < 0) printf("Error reading randomData in devrandomseed\n");
	// you now have a random integer!
	close(randomData);

	printf("in devrandomseed, seed is %d (%u)\n",sjhRandomInteger,sjhRandomInteger);
	return sjhRandomInteger;
}





int initmyrand(int seed){

	if(seed<0){
#ifdef USING_SEED_DEVRAND
		seed = devrandomseed();
#else
		seed = time(NULL);
#endif
	}

#ifdef VERBOSE
	printf("in initmyrand, seed is %d (%u)\n",seed,seed);
#endif

#ifdef USING_MT
	sgenrand(seed);
#else
	srand(seed);
#endif
	return seed;
}




unsigned long longinitmyrand(unsigned long *inseed){

	unsigned long seed;
	if(inseed==NULL){
#ifdef USING_SEED_DEVRAND
		seed = devrandomseed();
		printf("in initmyrand, seed is %ld", (long int) seed  );
		printf("(%lu)\n",(unsigned long int) seed);
#else
		seed = time(NULL);
#endif
	}
	else{
		seed = *inseed;
	}
#ifdef USING_MT
	sgenrand(seed);
#else
	srand(seed);
#endif

	return seed;
}





double rand0to1(){
	double x;
#ifdef USING_MT
	x = genrand();
#else
	x = (double) rand() / (double) RAND_MAX;
#endif
	return x;
}

unsigned long randint(){

	unsigned long x;

#ifdef USING_MT
	x = genrandint();
#else
	x = rand();
#endif
	return x;
}




/* Create an array of random integers between the range [min,max) */
int * randintarray(const int size,const int Min,const int max){
	int i, * array;
	array = (int *) malloc(size*sizeof(int));
	for(i=0;i<size;i++)
		array[i] = Min + floor((double)max*rand0to1());
	return array;
}




/* Create an array of random integers between the range [min,max) */
int * randboolarray(const int size){
	int i, * array;
	float v;
	array = (int *) malloc(size*sizeof(int));
	for(i=0;i<size;i++){
		v=rand0to1();
		if(v<0.5)
			array[i] = 0;
		else
			array[i] = 1;
	}
	return array;
}


/*UTILITY FUNCTION FOR RE-SEEDING ON RESTART*/
int get_mti(){
#ifdef USING_MT
	return mt_get_mti();
#else
	printf("NOT USING MERSENNE TWISTER - CAN'T GET MTI!!\n");
	return 0;
#endif
}

/* Only used in tests...see load_mt() in mt19937-2.c
void set_mti(int val){
#ifdef USING_MT
	mt_set_mti(val);
#else
	printf("NOT USING MERSENNE TWISTER - CAN'T SET MTI!!\n");
#endif
}*/


