/* Copyright (C) 2009-2015 Simon Hickinbotham                           */
/* When you use this, send an email to: sjh436@gmail.com                */
/* with an appropriate reference to your work.                          */

/* This file is part of STRINGMOL										*/

/* STRINGMOL is free software: you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* STRINGMOL is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License    */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.*/
/************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

/*TODO: To many interdependencies here - stringPM.h requires agests_base.h requires rules.h ...*/
//utilities
#include "mt19937-2.h"
#include "randutil.h"

//abm stuff
#include "rules.h"
#include "agents_base.h"

//stringmol
#include "params.h"
#include "alignment.h"
#include "SMspp.h"
#include "stringPM.h"


// Writing PNGs
#include "lodepng.h"
#include <iostream>

#include "setupSM.h"

#include "tests.h"



/*Test the random number generator. It must:
 *
 * 1: give the same sequence from a seed
 * 2: give the current position so that we can re-seed
 */
int test_rand(int verbose){

	int failed = 0;
	int rin = 436;

	if(verbose)printf("Testing seeding using dev/rand... \n");fflush(stdout);
	int rout = initmyrand(-1);
	if(rout == rin){
		printf("FAILED - requested seed from dev/random, but got -1\n");
	}
	else
		if(verbose)printf("PASSED - init set seed as %d (%u)\n",rout,(unsigned int) rout);

	if(verbose)printf("Testing seeding using dev/rand again... \n");fflush(stdout);
	rout = rout - initmyrand(-1);
	if(rout == 0 ){
		printf("FAILED - requested new seed from dev/random, but got same one\n");
	}
	else
		if(verbose)printf("PASSED - init set seed as %d\n",rout);

	if(verbose)printf("Testing seeding using ingeter %d... \n",rin);fflush(stdout);
	rout = initmyrand(rin);
	if(rout != rin){
		printf("FAILED - seed not set - different seed used\n");
	}
	else
		if(verbose)printf("PASSED - init set seed as %d\n",rout);


	if(verbose)printf("Testing re-setting mt index... \n");fflush(stdout);
	int pos = 22;
	set_mti(pos);
	pos = get_mti();




	if(!failed)
		printf("ALL RNG TESTS PASSED\n\n");
	return failed;
}



//stringPM * test_config_settings( int argc, char *argv[], int return_SM){
stringPM * test_config_settings( int argc, char *argv[], int return_SM){
	/** The idea here is to report the default values of the parameters, then parse the config and report them again. */

	stringPM *A;

	A = new stringPM(NULL);
	unsigned int ntrials = 0;
	unsigned int nsteps = 0;

	printf("\nBEFORE loading the config, params are:\n");
	print_params(A,ntrials,nsteps);
	printf("..c'est ca!\n\n");

	//int readordef_param_int(char *fn, const char *label, int *val, const int defaultvalue, const int verbose)
	readordef_param_int(argv[2], "NTRIALS", &ntrials, 1, 1);
	int nns = readordef_param_int(argv[2], "NSTEPS", &nsteps, -1, 1);

	A->load(argv[2],NULL,0,1);
	//if(!arg_load(A, argc, argv, 0))
	//	return NULL;

	printf("\n\nAFTER loading the config, params are:\n");
	print_params(A,ntrials,nsteps);
	if(nns==1)
		printf("NSTEPS was not specified. Simulations will run indefinitely");
	printf("..c'est ca!\n\n");

	if(return_SM)
		return A;
	else{
		A->clearout();
		delete A;
		return NULL;
	}
}

int compare_config(stringPM *A, stringPM *B){

	/*TODO: compare non-stringPM variables
	printf("Non-stringPM variables:\n");
	if(ntrials<0)
		printf("NTRIALS     not set - the default value would be used if needed\n");
	else
		printf("NTRIALS     %d\n",ntrials);
	*/

	//load params:
	if(A->cellrad-B->cellrad>FLT_MIN){
		printf("ERROR - cellrad not saved properly\n");
	}
	if(A->vcellrad-B->vcellrad>FLT_MIN){
		printf("ERROR - vcellrad not saved properly\n");
	}
	if(A->energy!=B->energy){
		printf("ERROR - energy not saved properly\n");
	}
	if(A->nsteps!=B->nsteps){
		printf("ERROR - nsteps not saved properly\n");
	}


	/** Although agents_base contains a call to load_influx, Stringmol never uses it */
	//TODO: test blosum - probably needs its own function
	//load agents: load_table(_mtx); load_mut; load_decay
	//if(A->blosum == NULL)
	//	printf("BLOSUM      not set - needs to be loaded explicitly\n");
	//else
	//	printf("BLOSUM      %d size table loaded\n",A->blosum->N);


	if(A->indelrate-B->indelrate>FLT_MIN){
		printf("ERROR - indelrate not saved properly\n");
	}
	if(A->subrate-B->subrate>FLT_MIN){
		printf("ERROR - subrate not saved properly\n");
	}
	if(A->decayrate-B->decayrate>FLT_MIN){
		printf("ERROR - decayrate not saved properly\n");
	}
	if(A->maxl0!=B->maxl0){
		printf("ERROR - maxl0 not saved properly\n");
	}
	if(A->estep!=B->estep){
		printf("ERROR - estep not saved properly\n");
	}

	return 0;

}


/* Test loading and saving of configs..
 * STRATEGY:
 * 		1: Load a file with known settings - see if we've got the right number.
 * 		1a: Save the file and see if the settings are the same
 * 		2: Iterate 10,000 time steps.
 * 		3: Save state.
 * 		4: Load state in a new object
 * 		5: Compare states.
 */
int test_loadsave(int argc, char *argv[]){

	/*TODO: test arguments */
	stringPM *A;
	stringPM *B;
	stringPM *C;
	const int fnlen =200;
	FILE *fp;
	char fn[] = "test_output.cfg";
	char fn1000[] = "test_output_1000.cfg";
	char **argv2;

	//Load the simulation and test that the config settings are correct
	A = test_config_settings(argc,argv,1);

	//Write the resulting config to file
	fp = fopen(fn,"w");
	A->print_conf(fp);
	fclose(fp);

	argv2 = (char **) malloc(argc*sizeof(char *));
	for(int c=0;c<argc;c++){
		argv2[c] = (char *)malloc(fnlen*sizeof(char));
		memset(argv2[c],0,fnlen*sizeof(char));
		sprintf(argv2[c],"%s",argv[c]);
	}
	sprintf(argv2[2],"%s",fn);

	//Load the simulation and test that the config settings are correct
	B = test_config_settings(argc,argv2,1);

	int csc = compare_config(A,B);

	//Run the Trial forward

	A->print_agents(stdout,"NOW",0);


	run_one_AlifeXII_trial(A);

	//Write the resulting config to file
	fp = fopen(fn1000,"w");
	A->print_conf(fp);
	fclose(fp);

	sprintf(argv2[2],"%s",fn1000);


	C = test_config_settings(argc,argv2,1);


	csc = compare_config(A,C);


	return csc;
}




/* Refactoring practice demands that tests are written whilst writing code!
 * This function tests that each stringmol configuration works
 */
int test_all(int argc, char *argv[]){

	int failed = 0;

	printf("Testing rng\n");
	failed = test_rand(0);


	failed = test_loadsave(argc,argv);



	printf("Check config test\n");



	return failed;
}


void test_rand_config(int argc, char *argv[]){
	int j=0;

	//Seed the rng
	//TODO: Get the Seed from the config file...
	sgenrand(1); //any nonzero integer can be used as a seed

	/* The size of the MT array is 624, but mti is incremented
	 * *after* being used, so the range of values of mti is
	 * between 1 and 624...
	 */
	int lim = 625;
	for (j=0; j<lim; j++) {
		printf("%0.8f ", rand0to1());//genrand());
		if (j%8==7) printf("\n");
	}
	printf("\nMTSTATE IS:\n");

	//Save the state:
	print_mt(stdout);
	FILE *fp,*tfp1,*tfp2;
	fp = fopen("TestMT1.dat","w");
	print_mt(fp);
	fclose(fp);

	//Run the rng on another 1000, save output to file
	fp = fopen("testdata1.dat","w");
	for (j=0; j<lim; j++) {
		fprintf(fp, "%0.8f", rand0to1());//genrand());
		if (j%8==7)fprintf(fp,"\n");else fprintf(fp,"\t");
	}
	fclose(fp);

	load_mt("TestMT1.dat");


	fp = fopen("testdata1_again.dat","w");
	for (j=0; j<lim; j++) {
		fprintf(fp, "%0.8f", rand0to1());//genrand());
		if (j%8==7)fprintf(fp,"\n");else fprintf(fp,"\t");
	}
	fclose(fp);

	//Just to be sure, let's run this a few hundred thousand times and add a prime...
	for (j=0; j<((lim*1000)+13); j++) {
		rand0to1();
	}

	//Now see if the two sets of values are the same:
	tfp1 = fopen ("testdata1.dat","r");
	tfp2 = fopen ("testdata1_again.dat","r");
	//TODO: write this test (doing it manually atm)

	fclose(tfp1);
	fclose(tfp2);




}

