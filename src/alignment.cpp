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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "randutil.h"

#include "alignment.h"

/*
 * LINKED LIST STUFF
 * ***************************************************************/

//swa = read_sw(swlist,a1->spp,a2->spp);

s_sw * read_sw(s_sw *swlist, int sp1, int sp2){

	s_sw *p;

	for(p=swlist;p!=NULL;p=p->next){
		if(p->sp1 == sp1)
			if(p->sp2 == sp2)
				return p;
	}
	return NULL;
}



//store_sw(&swlist,sw);

int store_sw(s_sw **head, align *sw, int sp1, int sp2){

	//idea is to bung new alignments at the front of the list, since they are more likely to be used..

	s_sw *p;
	s_sw *old;


	old = *head;

	p=(s_sw *) malloc(sizeof(s_sw));

	p->next = old;
	p->sp1 = sp1;
	p->sp2 = sp2;

	p->match = sw->match;		// the number of matching characters.
	p->score = sw->score; 		// the score of the match
	p->prob =  sw->prob;		// the probability of the match - used for determining events based on the score/match
	p->s1 =    sw->s1;			// start of the match in string 1
	p->e1 =    sw->e1;			// end of the match in string 1
	p->s2 =    sw->s2;			// start of the match in string 2
	p->e2 =    sw->e2;			// end of the match in string 2

	*head = p;




	//if(*head!=NULL)
	return 0;
}


//load_sw(swa,sw);
int load_sw(s_sw *b, align *sw){
	sw->match = b->match;		// the number of matching characters.
	sw->score = b->score; 		// the score of the match
	sw->prob =  b->prob;		// the probability of the match - used for determining events based on the score/match
	sw->s1 =    b->s1;			// start of the match in string 1
	sw->e1 =    b->e1;			// end of the match in string 1
	sw->s2 =    b->s2;			// start of the match in string 2
	sw->e2 =    b->e2;			// end of the match in string 2
	return 0;
}



void free_swlist(s_sw **head){

	s_sw *p,*pold;

	pold=*head;
	while(pold!=NULL){
		p = pold->next;
		free(pold);
		pold=p;
	}
	*head = NULL;
}


////////////////////////////////////////////////////////////////////////////////////


//public
int LongestCommonSubsequence(char *s1, char *s2)
{
	int l1,l2;
	int i,j;

	//if either string is empty, the length must be 0
	if ((l1=strlen(s1))==0 || (l2=strlen(s2))==0)
		return 0;

	int **num;
	num = (int **) malloc(l1*sizeof(int *));// new int[s1.Length, s2.Length];  //2D array
	for(i=0;i<l1;i++){
		num[i]=(int *) malloc(l2*sizeof(int));
		memset(num[i],0,l2*sizeof(int));
	}


	char letter1;
	char letter2;

	//Actual algorithm
	for(i=0;i<l1;i++){
		letter1 = s1[i];
		for(j=0;j<l2;j++){
			letter2 = s2[j];

			if(letter1 == letter2)
			{
				if((i == 0) || (j == 0))
					num[i][j] = 1;
				else
					num[i][j] = 1 + num[i-1][j-1];
			}
			else
			{
				if ((i == 0) && (j == 0))
					num[i][j] = 0;
				else if ((i == 0) && !(j == 0))   //First ith element
					num[i][j] = 0>num[i][j - 1]?0:num[i][j - 1];
				else if (!(i == 0) && (j == 0))   //First jth element
					num[i][j] = 0>num[i-1][j]?0:num[i-1][j];
				else // if (!(i == 0) && !(j == 0))
					num[i][j] = num[i-1][j]>num[i][j-1]?num[i-1][j]:num[i][j-1];
			}
		}//end j
	}//end i

	//let's printout the sequences and look at the numbers...

	void align_prob(align *A);

	//printf("\t");
	printf(" ");
	for(i=0;i<l1;i++)
		printf("%c",s1[i]);
	printf("\n");
	for(j=0;j<l2;j++){
		printf("%c",s2[j]);
		for(i=0;i<l1;i++){
			printf("%d",num[i][j]);
		}
		printf("\n");
	}

	//ok - now let's try and find where the best match is now...


	//We know where the first matching character is in the sequence, but is is possible that there may be a match closer down...


	return num[l1-1][l2-1];
} //end LongestCommonSubsequence
//Usage: LongestCommonSubsequence("computer", "boathouse")




int swt_index(int c,swt *T){
	int i = 0;

	for(i=0;i<T->N;i++){
		if(T->key[i]==c)
			return i;
	}

	printf("Code %d (%c) not found in SWT key\n",c,c);
	fflush(stdout);
	return -1;

}




float wT(int a, int b,swt *T){

	int ai,bi;

	if(a)
		ai = swt_index(a,T);
	if(b)
		bi = swt_index(b,T);

	if(a<0 || b<0){
		printf("swt_index returned an error w()\n");fflush(stdout);
		return -10000;
	}

	if(!a && !b){// ERROR
		printf("BAD CALL TO w()\n");fflush(stdout);
		return -10000;
	}
	if(!b){//Indel with a
		return T->T[T->N][ai];
	}
	if(!a){//Indel with b
		return T->T[T->N][bi];
	}

	//
	return T->T[ai][bi];

}


float w(int a, int b){

	float mismatch = -1.0/3.0;

	if(!a && !b){// ERROR
		printf("BAD CALL TO w()\n");fflush(stdout);
		return -10000;
	}

	if(!b){//DELETION
		return mismatch;
	}

	if(!a){//INSERTION
		return mismatch;
	}

	if(a==b)
		return 1.0;
	else
		return mismatch;
}



int SmithWaterman(char *s1, char *s2, align *A, swt *T, int verbose){

	//align *A;
	int l1,l2;
	int i,j;
	int si,sj;
	int
			match,
		endi=0,
		endj=0;
	float
		m,
		**H,
		emax = 0;
	float l,d,t;
	int action;

	//Set up the structures and test for errors.

	//if either string is empty, the length must be 0
	if ((l1=strlen(s1))==0 || (l2=strlen(s2))==0)
		return 0;

	H = (float **) malloc((l1+1)*sizeof(float *));// new int[s1.Length, s2.Length];  //2D array
	for(i=0;i<=l1;i++){
		H[i]=(float *) malloc((l2+1)*sizeof(float));
		memset(H[i],0,(l2+1)*sizeof(int));
	}


	//Populate H
	for(i=1;i<=l1;i++){
		si=i-1;
		for(j=1;j<=l2;j++){
			sj=j-1;
			//Calculate the three values:
			H[i][j]=0;

			//Match/Mismatch
			//m = H[i-1][j-1] + w(s1[si],s2[sj]);
			m = H[i-1][j-1] + wT(s1[si],s2[sj],T);
			H[i][j] = m>H[i][j]?m:H[i][j];

			//Deletion

			//m = H[i-1][j] + w(s1[si],0);
			m = H[i-1][j] + wT(s1[si],0,T);
			H[i][j] = m>H[i][j]?m:H[i][j];

			//Insertion
			//m = H[i][j-1] + w(0,s2[sj]);
			m = H[i][j-1] + wT(0,s2[sj],T);
			H[i][j] = m>H[i][j]?m:H[i][j];

			//record the posn of the highest score:
			if(H[i][j]>emax){
				emax = H[i][j];
				endi = i;
				endj = j;
			}
		}
	}

	//OK - now we know the max score and the end posn. Need to track back from the end...

	i=endi;
	j=endj;
	match = 0;

	if(verbose){
		printf("SW alignment table\n");
		printf("\t");
		for(i=0;i<=l1;i++){
			printf("\t%c",s1[i]);
		}
		printf("\n");
		fflush(stdout);
		for(j=0;j<=l2;j++){
			printf("%c",s2[j]);
			for(i=0;i<=l1;i++){
				printf("\t%0.2f", H[i][j]);
			}
			printf("\n");
		}
	}



	while(H[i][j]>0){
		l=H[i-1][j];
		d=H[i-1][j-1];
		t =H[i][j-1];

		action = -1;
		if(l>=d&&l>=t){
			action = 0; //LEFT 	(INSERTION)
		}
		if(t>=d&&t>=l){
			action = 2; //TOP (DELETION)
		}
		if(d>=l&&d>=t){
			action = 1; //DIAGONAL (MATCH)
		}
		switch(action){
			case -1:
				printf("no action found??\n");
				break;
			case 0:
				i--;
				break;
			case 1:
				match++;
				i--;
				j--;
				break;
			case 2:
				j--;
				break;
		}
	}

	A->s1=i;
	A->s2=j;
	A->e1=endi;
	A->e2=endj;
	A->score = emax;
	A->match = match;

	for(i=0;i<=l1;i++){
		free(H[i]);
	}
	free(H);

	return 0;
}





/*
 * THIS VERSION USES A PROPER TRACE BACK MATRIX
 */
int SmithWatermanV2(char *s1, char *s2, align *A, swt *swT, int verbose){

	//align *A;
	int l1,l2;
	int i,j;
	int si,sj;
	int
		match,
		endi=0,
		endj=0,
		**T;
	float
		m,
		**H,
		emax = 0;

	//Set up the structures and test for errors.

	//if either string is empty, the length must be 0
	if ((l1=strlen(s1))==0 || (l2=strlen(s2))==0)
		return 0;

	H = (float **) 	malloc((l1+1)*sizeof(float *));// new int[s1.Length, s2.Length];  //2D array
	T = (int **) 	malloc((l1+1)*sizeof(int *  ));
	for(i=0;i<=l1;i++){
		H[i]=(float *) malloc((l2+1)*sizeof(float));
		T[i]=(  int *) malloc((l2+1)*sizeof(  int));

		memset(H[i],0,(l2+1)*sizeof(float));
		memset(T[i],0,(l2+1)*sizeof(  int));
	}


	//Populate H
	for(i=1;i<=l1;i++){
		si=i-1;
		for(j=1;j<=l2;j++){
			sj=j-1;
			//Calculate the three values:
			H[i][j]=0;

			//Match/Mismatch
			//m = H[i-1][j-1] + w(s1[si],s2[sj]);
			m = H[i-1][j-1] + wT(s1[si],s2[sj],swT);
			if(m>H[i][j]){
				H[i][j] = m;
				T[i][j] = swMATCH;
			}

			//Deletion
			m = H[i-1][j] + wT(s1[si],0,swT);
			if(m>H[i][j]){
				H[i][j] = m;
				T[i][j] = swDEL;
			}

			//Insertion
			m = H[i][j-1] + wT(0,s2[sj],swT);
			if(m>H[i][j]){
				H[i][j] = m;
				T[i][j] = swINS;
			}

			//record the posn of the highest score:
			//NOTE this is the topleft high value - we may want bottomright
			if(H[i][j]>emax){
				emax = H[i][j];
				endi = i;
				endj = j;
			}
		}
	}

	//OK - now we know the max score and the end posn. Need to track back from the end...

	i=endi;
	j=endj;
	match = 0;

	if(verbose){
		printf("SW alignment table\n");
		printf("\t");
		for(i=0;i<=l1;i++){
			printf("\t%c",s1[i]);
		}
		printf("\n");
		fflush(stdout);
		for(j=0;j<=l2;j++){
			printf("%c",s2[j]);
			for(i=0;i<=l1;i++){
				printf("\t%0.2f", H[i][j]);
			}
			printf("\n");
		}

		printf("SW trace back\n");
		printf("\t");
		for(i=0;i<=l1;i++){
			printf("\t%c",s1[i]);
		}
		printf("\n");
		fflush(stdout);
		for(j=0;j<=l2;j++){
			printf("%c",s2[j]);
			for(i=0;i<=l1;i++){
				printf("\t%d", T[i][j]);
			}
			printf("\n");
		}
	}

	//printf("\nfinished\n");

	i=endi;
	j=endj;
	match = 0;

	while(T[i][j]){
		switch(T[i][j]){
			case swMATCH:
				match++;
				i--;
				j--;
				break;
			case swDEL:
				i--;
				break;
			case swINS:
				j--;
				break;
			case 0:
				if(verbose)
					printf("Trace runs from %d,%d to %d,%d\n",endi,endj,i,j);
				break;
			default:
				printf("Bad value %d for trace back!\n",T[i][j]);

		}
	}

	if(i<0||i>endi)
		printf("ERROR in traceback for i = %d, endi = %d, l1 = %d\n",i,endi,l1);
	if(j<0||j>endj)
		printf("ERROR in traceback for j = %d, endj = %d, l2 = %d\n",j,endj,l2);

	A->s1=i;
	A->s2=j;
	A->e1=endi;
	A->e2=endj;
	A->score = emax;
	A->match = match;

	for(i=0;i<=l1;i++){
		free(H[i]);
		free(T[i]);
	}
	free(H);
	free(T);

	return 0;
}















void align_prob(align *A){

	if(A->match)
		A->prob = A->score/A->match;
	else
		A->prob = 0.;

}


int align_event(align *A,int len){

	float rand = (float) rand0to1();
	if(len<0){
		printf("Align with length -1\n");
		fflush(stdout);
		align_prob(A);
	}
	else
		A->prob = pow((float)A->score/len,len);//factorial(len);//A->score/len;
	if(rand<A->prob)
		return 1;
	return 0;
}


float instr_wt(char C){

	if(C>='A' && C<='Z')
		return 1.0;
	return 0.5;
}


void print_swt(FILE *fp, swt *sss){
	int i,j;
	printf("   ");
	for(i=0;i<sss->N;i++)
		printf("%c   ",sss->key[i]);
	printf("\n");
	for(i=0;i<sss->N;i++){
		printf("%c  ",sss->key[i]);
		for(j=0;j<sss->N;j++){
			printf("%0.2f ",sss->T[i][j]);
		}
		printf("\n");
	}
	//Print the indel value
	printf("-   ");
	for(i=0;i<sss->N;i++){
		printf("%0.2f ", sss->T[sss->N][i]);
	}
	printf("\n");
}


char sym_from_adj(char X, swt *swt){

	int idx = tab_idx(X,swt);
	int count=0,i;
	float rno=rand0to1();

	if(idx>-1){
		for(i=0;i<swt->N;i++)
			if(swt->adj[idx][i]==1)
				count++;

		count = ceil((float)count*rno);
		i=0;
		while(count){
			if(swt->adj[idx][i++]==1)
				count--;
		}
		return swt->key[i-1];
	}
	else{
		printf("ERRROR reading index for %c\n",X);
	}


	return 0;
}



void table_from_string(float **T, char *key, const int N){

	int i,j,d1;
	int sum;

//	i=0;
	//Step 1: fill the mutations.



	i=0;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			d1 = abs(i-j);
			T[i][j] = -1. * (d1<N-d1?d1:N-d1);
		}
	}


	i=0;
	//step 2: Calculate the diagonal
	sum=0;
	for(i=0;i<N;i++)
		sum += T[i][0];
	sum/=N;


	i=0;
	float indval = -4./3.;
	//step 3: Calculate the indel
	for(i=0;i<N;i++){
		T[i][i]= -sum * instr_wt(key[i]);
		T[N][i]= indval;
	}



	i=0;
	//step 4: Normalize over the diagonal.
	//NB: Indel row is already done in step 3!
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			T[i][j] /= (float) abs(sum);


}


int load_table(char *fn,swt *T){

	const int maxl=256;
	FILE *fp,*fp2;
	char fn2[maxl];
	char line[maxl];
	char label[maxl];
	int i,j,found;

	T->N=0;
	T->T=NULL;
	T->key=NULL;

	if((fp=fopen(fn,"r"))!=NULL){
		found = 0;
		while((fgets(line,maxl,fp))!=NULL){
			memset(label,0,maxl);
			sscanf(line,"%s",label);
			//printf("line = %s",line);
			if(!strncmp(line,"USING",5)){
				sscanf(line,"%*s %s",fn2);
				found = 1;
				fclose(fp);
				break;
			}
		}
		if(!found){
			fclose(fp);
			printf("No MIS file name found\n");
			return 1;
		}



		if((fp2=fopen(fn2,"r"))!=NULL){
			found = 0;
			while((fgets(line,maxl,fp))!=NULL){
				memset(label,0,maxl);
				sscanf(line,"%s",label);
				//printf("line = %s",line);
				if(!strncmp(line,"SET",3)){
					sscanf(line,"%*s %s",fn2); //using fn2 to temporarily hold the alphabet...
					T->N = strlen(fn2);
					T->key = (char *)malloc((T->N+1) * sizeof(char));
					memset(T->key,0,(T->N+1)*sizeof(char));
					strncpy(T->key,fn2,T->N);

					T->T = (float **) malloc((T->N +1) * sizeof(float *));
					for(i=0;i<T->N+1;i++)
						T->T[i] = (float *) malloc(T->N * sizeof(float));

					table_from_string(T->T,T->key,T->N);

					//print for debug.
					printf("   ");
					for(i=0;i<T->N;i++)
						printf("%c   ",T->key[i]);
					printf("\n");
					for(i=0;i<T->N;i++){
						printf("%c  ",T->key[i]);
						for(j=0;j<T->N;j++){
							printf("%0.2f ",T->T[i][j]);
						}
						printf("\n");
					}
					//Print the indel value
					printf("-   ");
					for(i=0;i<T->N;i++){
						printf("%02d ",(int) T->T[i][T->N]);
					}
					printf("\n");


					found = 1;
					fclose(fp2);

					//SUCCESS
					return 0;
				}
			}
			if(!found){
				fclose(fp2);
				printf("No MIS string found in MIS file\n");
				return 2;
			}

		}
		else{
			printf("Unable to open MIS file %s",fn2);
			return 3;
		}
	}
	else
		return 4;

	return 5;
}



int tab_idx(char X, swt *T){

	int i;
	for(i=0;i<T->N;i++)
		if(T->key[i]==X)
			return i;

	printf("ERROR: Code %c (%d) not found in instruction set!\n",X,X);fflush(stdout);
	return -1;
}

/* FUNCTIONS FOR ADAM NELLIS */

float score_sw(char *s1, char *s2, swt *swT){
	align A;
	//SmithWatermanV2(char *s1, char *s2, align *A, swt *swT, int verbose)
	SmithWatermanV2(s1,s2,&A,swT,0);
	return A.score;
}

/////////TESTING

void test_adj(swt *swt){
	int i,j,ntrials=10;
	printf("Adjacency matrix:\n");
	printf("  %s\n",swt->key);
	for(i=0;i<swt->N;i++){
		printf("%c ",swt->key[i]);
		for(j=0;j<swt->N;j++){
			if(swt->adj[i][j])
				printf("*");
			else
				printf(" ");
		}
		printf("\n");
	}
	printf("Testing adjacent recall\n");
	for(i=0;i<swt->N;i++){
		printf("%c; ",swt->key[i]);
		for(j=0;j<ntrials;j++){
			printf("%c,",sym_from_adj(swt->key[i],swt));
		}
		printf("\n");
	}
}



/* We want the default alignment to be:
 *


ABC$DEF%GH^IJK?LMN}OPQ>RST=UVWXYZ


1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121
-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243
-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364
-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485
-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607
-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728
-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849
-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971
-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092
-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213
-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335
-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456
-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577
-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699
-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820
-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941
-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941
-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820
-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699
-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577
-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456
-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335
-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213
-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092
-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971
-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849
-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728
-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607
-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485
-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364
-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243
-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121
-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000
-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333
0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0

*/



swt * default_table(){

	swt * table;
	char ** TData;
	const int tlen = 500;
	int i,j;
	char *p;

	table = (swt *) malloc(sizeof(swt));

	//Set N and key:
	table->N = 33;
	table->key = (char *)malloc((table->N+1) * sizeof(char));
	memset(table->key,0,(table->N+1)*sizeof(char));

	strncpy(table->key,"ABC$DEF%GH^IJK?LMN}OPQ>RST=UVWXYZ",table->N);

	//Set up the data structure
	table->T = (float **) malloc((table->N +1) * sizeof(float *));
	TData = (char **)malloc((table->N +1) * sizeof(char *));
	for(i=0;i<table->N+1;i++){
		table->T[i] = (float *) malloc(table->N * sizeof(float));
		TData[i] = (char *) malloc(tlen *sizeof(char));
		memset(TData[i],0,tlen *sizeof(char));
	}

	strcpy(TData[ 0], "1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121");
	strcpy(TData[ 1],"-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243");
	strcpy(TData[ 2],"-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364");
	strcpy(TData[ 3],"-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485");
	strcpy(TData[ 4],"-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607");
	strcpy(TData[ 5],"-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728");
	strcpy(TData[ 6],"-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849");
	strcpy(TData[ 7],"-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971");
	strcpy(TData[ 8],"-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092");
	strcpy(TData[ 9],"-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213");
	strcpy(TData[10],"-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335");
	strcpy(TData[11],"-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456");
	strcpy(TData[12],"-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577");
	strcpy(TData[13],"-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699");
	strcpy(TData[14],"-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820");
	strcpy(TData[15],"-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941");
	strcpy(TData[16],"-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941");
	strcpy(TData[17],"-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820");
	strcpy(TData[18],"-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699");
	strcpy(TData[19],"-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577");
	strcpy(TData[20],"-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456");
	strcpy(TData[21],"-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335");
	strcpy(TData[22],"-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213");
	strcpy(TData[23],"-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092");
	strcpy(TData[24],"-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971");
	strcpy(TData[25],"-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849");
	strcpy(TData[26],"-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	0.500	-0.121	-0.243	-0.364	-0.485	-0.607	-0.728");
	strcpy(TData[27],"-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485	-0.607");
	strcpy(TData[28],"-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364	-0.485");
	strcpy(TData[29],"-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243	-0.364");
	strcpy(TData[30],"-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121	-0.243");
	strcpy(TData[31],"-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000	-0.121");
	strcpy(TData[32],"-0.121	-0.243	-0.364	-0.485	-0.607	-0.728	-0.849	-0.971	-1.092	-1.213	-1.335	-1.456	-1.577	-1.699	-1.820	-1.941	-1.941	-1.820	-1.699	-1.577	-1.456	-1.335	-1.213	-1.092	-0.971	-0.849	-0.728	-0.607	-0.485	-0.364	-0.243	-0.121	1.000");
	strcpy(TData[33],"-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333	-1.333");

	//Now read the data into the table:
	for(i=0;i<table->N+1;i++){
		//fgets(line,maxl,fp);
		p=strtok(TData[i],", \t");
		for(j=0;j<table->N;j++){
			sscanf(p,"%f",&(table->T[i][j]));
			p = strtok(NULL,", \t");
		}
	}

	table->adj = (int **) malloc((table->N) * sizeof(int *));
	for(i=0;i<table->N;i++){
		table->adj[i] = (int *) malloc(table->N * sizeof(int));
		memset(TData[i],0,tlen *sizeof(char));
	}

	strcpy(TData[ 0],  "0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1");
	strcpy(TData[ 1],  "1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[ 2],  "0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[ 3],  "0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[ 4],  "0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[ 5],  "0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[ 6],  "0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[ 7],  "0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[ 8],  "0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[ 9],  "0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[10],  "0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[11],  "0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[12],  "0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[13],  "0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[14],  "0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[15],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[16],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[17],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[18],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[19],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[20],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[21],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0	0");
	strcpy(TData[22],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0	0");
	strcpy(TData[23],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0	0");
	strcpy(TData[24],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0");
	strcpy(TData[25],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0	0");
	strcpy(TData[26],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0	0");
	strcpy(TData[27],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0	0");
	strcpy(TData[28],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0	0");
	strcpy(TData[29],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0	0");
	strcpy(TData[30],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1	0");
	strcpy(TData[31],  "0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	1");
	strcpy(TData[32],  "1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0");


	for(i=0;i<table->N;i++){
		//fgets(line,maxl,fp);
		p=strtok(TData[i],", \t");
		for(j=0;j<table->N;j++){
			sscanf(p,"%d",&(table->adj[i][j]));
			p = strtok(NULL,", \t");
		}
		free(TData[i]);
	}

	free(TData);

	return table;
}












