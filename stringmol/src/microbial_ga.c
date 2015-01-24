/* Copyright (C) 2009-2012 Verena Fischer and Simon Hickinbotham        */
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
#include <stdio.h>
#include <math.h>

#include "microbial_ga.h"
#include "randutil.h"

/* GA parameters */
//const int POPSIZE = 100;
//const int RUNS = 100;
//const int PARAMETERS = 25; //25 rules in the metabolism

const double REC = 0.5;
const double MUT = 0.1;
double MDIST = 0.5; //0.05;
/* GA parameters END*/


/*	 COMMENTED THESE OUT - NOW PASSED INTO THE PROGRAM
double pop[POPSIZE][PARAMETERS];
double eval[POPSIZE];
*/

/*
 * "Mutates" a double value, i.e. adds or subtracts MDIST to/from the value.
 */
double mutate (double f) {
	if (randint()%2 == 0){
		f = (f + MDIST);
	}
	else{
		f = (f - MDIST);
	}
	return f;
}

/*
 * "Mutates" an integer value, by doing the equivalent of a bitflip...
 * NB: not very efficient at the moment!
 */
int mutate_int (int val, int min, int max) {

	int range = max-min;
	int brange = 1,rpow = 0;
	int mpow;

	//get the powers of 2 that span the range:
	while(brange<range){
		brange *= 2;
		rpow++;
	}

	int outofrange=1;
	while(outofrange){
		mpow = (float) rpow * rand0to1();

		if (randint()%2 == 0){
			val = (val + pow(2,mpow));
		}
		else{
			val = (val - pow(2,mpow));
		}
		if(val>=min && val<max)
			outofrange=0;
	}
	return val;
}

/**
 * Microbial GA step.
 */
//void step() {
int ga_step(double **pop, double *eval, const int POPSIZE, const int PARAMETERS){
	int w;
	int l;
	int i;
	int a=randint()%POPSIZE;
	int b=randint()%POPSIZE;
	while (a==b) {
		b=randint()%POPSIZE;
	}
	if (eval[a]<eval[b]) {//note the lower the value the fitter the individual!
		w=a;
		l=b;
	}
	else {
		w=b;
		l=a;
	}
	for (i=0; i<PARAMETERS; i++) {
		//recombine
		if (((randint()%1000)/1000.0)<REC) {
			pop[l][i]=pop[w][i];
		}
		//mutate
		if ((randint()%1000)/1000.0<MUT) {
			pop[l][i]=mutate(pop[l][i]);
		}
	}
	//eval[l] = evaluate(l);
	return l;
}




//void step() {
int ga_step_int(int **pop, double *eval, const int POPSIZE, const int PARAMETERS, const int minval, const int maxval, int *wn){
	int winner;
	int loser;
	int parameter;
	int a=randint()%POPSIZE;
	int b=randint()%POPSIZE;
	while (a==b) {
		b=randint()%POPSIZE;
	}
	if (eval[a]<eval[b]) {//note the lower the value the fitter the individual!
		winner=a;
		loser=b;
	}
	else {
		winner=b;
		loser=a;
	}
	for (parameter=0; parameter<PARAMETERS; parameter++) {
		//recombine
		if (((randint()%1000)/1000.0)<REC) {
			pop[loser][parameter]=pop[winner][parameter];
		}
		//mutate
		if ((randint()%1000)/1000.0<MUT) {
			pop[loser][parameter]=mutate_int(pop[loser][parameter],minval,maxval);
		}
	}
	//sometimes we want to know the winner too!
	if(wn!=NULL)
		*wn = winner;
	return loser;
}

void recmut_bool(int **pop, const int winner, const int loser, const int L){
	int p;
	for(p=0;p<L;p++){
		//recombine
		if (((randint()%1000)/1000.0)<REC) {
			pop[loser][p]=pop[winner][p];
		}
		//mutate
		if ((randint()%1000)/1000.0<MUT) {
			pop[loser][p]=1-pop[loser][p];
		}
	}
}










/**
 * (returns euclidian distance to goal point)
 */
double evaluate(int i) {
	//do things with Player Stage
	return 0.;
}


int * randinitint(const int N,const int min,const int upto){

	int i,*A;
	A = (int *) malloc(N*sizeof(int));

	for(i=0;i<N;i++){
		A[i] = min + (upto*rand0to1());
	}

	return A;

}




















