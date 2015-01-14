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


#ifdef __cplusplus
extern "C" {
#endif

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_


typedef struct td_s_sw{
	//First, equivalences to s_align
	int match;		// the number of matching characters.
	float score; 	// the score of the match
	float prob;		// the probability of the match - used for determining events based on the score/match
	int s1;			// start of the match in string 1
	int e1;			// end of the match in string 1
	int s2;			// start of the match in string 2
	int e2;			// end of the match in string 2

	//now the bits we need to store in a list:
	int sp1;		// First chosen species
	int sp2;		// Second chosen species

	//Now the bits we need when we are measuring distances between species

	struct td_s_sw *next;	//next struct.
} s_sw;

typedef struct s_align{
	int match;		// the number of matching characters.
	float score; 	// the score of the match
	float prob;		// the probability of the match - used for determining events based on the score/match
	int s1;			// start of the match in string 1
	int e1;			// end of the match in string 1
	int s2;			// start of the match in string 2
	int e2;			// end of the match in string 2
} align;

typedef struct s_pset{
	char *i[2];	//instruction pointer
	char *r[2]; //read pointer
	char *w[2]; //write pointer
	char *f[2];	//flow (loop) pointer
	int	it;
	int rt;
	int wt;
	int ft;
} pset;

typedef struct s_swt{
	int 	**adj;	//adjacency matrix values
	float  	**T;	//substitution matrix values
	char  	*key;	//Instruction set list
	int		N;		//Number of instructions.
} swt;


//Linked list stuff

s_sw * 	read_sw(s_sw *swlist, int sp1, int sp2);
int 	store_sw(s_sw **swlist, align * sw,int sp1, int sp2);
int		load_sw(s_sw *b, align *sw);
void 	free_swlist(s_sw **head);


//Important to have nonzero values for this enumerartion!
enum sw_subs{swMATCH=1,swDEL=2,swINS=3};


int LongestCommonSubsequence(char *s1, char *s2);

int SmithWaterman(char *s1, char *s2, align *A, swt *T, int verbose);
int SmithWatermanV2(char *s1, char *s2, align *A, swt *swT, int verbose);

void align_prob(align *A);
int align_event(align *A, int len);


void print_swt(FILE *fp, swt *sss);

int load_table(char *fn,swt *T);
void table_from_string(float **T, char *key, const int N);
int tab_idx(char X, swt*T);

//Get an adjacent symbol in the mutation space
char sym_from_adj(char X, swt *swt);


/* Get score only and don't worry about anything else! */

float score_sw(char *s1, char *s2, swt *swT);

/*TESTING:*/
void test_adj(swt *swt);

#endif /* ALIGNMENT_H_ */
#ifdef __cplusplus
}
#endif
