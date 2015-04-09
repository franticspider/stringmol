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

#ifndef SETUPSM_H_
#define SETUPSM_H_

	struct runparams{
		int indefinite;
		int maxnsteps;
	};

    void clearfiles( char *argv[]);
	void setupSMol(struct runparams &R, int argc, char *argv[]);
	void record_spp(stringPM *A);
	void printsppct(stringPM *A, int t);

	void setmaxcode(stringPM *A, int *maxcode);
	int run_one_comass_trial(const int rr, stringPM *A, int * params, struct runparams *R);

	void setmutnet(int * mutnet, swt *blosum);

	//count species in a containers nowhead
	float ctspp(stringPM *A, const int spp);

	//get stats for evolution cf seed community
	//float * evostats(char * Afn, stringPM *B,s_sw **spp_matches, float *class_score,float *self);
	void evostats(char * Afn, stringPM *B,s_sw **spp_matches, float *self, float *gvm);

	//Different ways of loading, dependant upon the no. of arguments
	int arg_load(stringPM *A, int argc, char *argv[], int verbose=0);


	float * self_stats(char * Afn);

	/* Premiminary check that a config is sane */
	void check_config( int argc, char *argv[]);

	/* Standard initialisation of random number seed */
	void init_randseed_config(int argc, char *argv[]);


#endif /* SETUPSM_H_ */
