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


/* Description
 *
 * This file contains helper functions to set up stringmol runs
 *
 * */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "randutil.h"
#include "params.h"
//}

//stringmol
#include "alignment.h"

//metabolism
#include "rules.h"
#include "agents_base.h"
#include "SMspp.h"
#include "stringPM.h"

#include "setupSM.h"

/* Used in comass_ga and comass_ga_boostwinners
 *
 */
void clearfiles(char *argv[]){

	char fn[256];
	FILE *ftmp;

	sprintf(fn,"%s.spatial.summary.dat",argv[1]);
	ftmp = fopen(fn,"w");
	fclose(ftmp);
	ftmp = fopen("epochs.dat","w");
	fclose(ftmp);
}



/* Used in comass_ga and comass_ga_boostwinners
 *
 */
void setupSMol(struct runparams &R, int argc, char *argv[]){

	FILE *fp;

	int rerr=1,rlim=20;
	if((fp=fopen(argv[2],"r"))!=NULL){
		rerr = read_param_int(fp,"NTRIALS",&rlim,1);
		switch(rerr){
		case 2:
			printf("Multiple NTRIALS specified. Check config file\n");
			getchar();
			exit(0);
		case 0:
			printf("Setting NTRIALS to %d\n",rlim);
			break;
		default:
			printf("NTRIALS not sepcifid;\nSetting NTRIALS to %d\n",rlim);
			break;
		}
		fclose(fp);
	}

	//Read nsteps:
	if((fp=fopen(argv[2],"r"))!=NULL){
		rerr = read_param_int(fp,"NSTEPS",&(R.maxnsteps),1);
		switch(rerr){
		case 2:
			printf("Multiple NSTEPS specified. Check config file\n");
			getchar();
			exit(0);
		case 0:
			printf("Setting NSTEPS to %d\n",R.maxnsteps);
			R.indefinite=0; //TODO: fix the indefinite thing if NSTEPS is not specified..
			break;
		default:
			printf("NSTEPS not sepcified;\nEach trial will run to extinction\n");
			break;
		}
		fclose(fp);
	}
	fflush(stdout);

	return;

}



/* Used in comass_ga and comass_ga_boostwinners
 *
 */
void record_spp(stringPM *A){
	char fn[256];
	FILE *fp;

	A->print_spp_count(stdout,0,-1);

	//printf("Printing species list\n");
	sprintf(fn,"splist%03d.dat",A->r);
	if((fp = fopen(fn,"w"))!=NULL){
		A->spl->print_spp_list(fp);
		fclose(fp);
	}
	else{
		printf("Unable to write to %s\n",fn);
	}
}

/* used in SmPm_conpop
 *
 */
float ctspp(stringPM *A, const int spp){

	float count = 0;

	s_ag *pa;


	for(pa = A->nowhead; pa!=NULL; pa = pa->next){
		if(pa->spp->spp==spp)
			count++;
	}

	return count;

}


/* used in comass_AlifeXII, energetic_AlifeXII, origlife,
 * comass_GA, comass_GA_boostwinners, SmPm_AlifeXII, SmPm_conpop and speigmonst
 *
 */
void printsppct(stringPM *A, int t){

	char fn[128];
	FILE *fp;
	s_ag *pa;
	int spc,count;
	int finished = 0;
	int nag,*done;
	int i,found;

	nag = A->nagents(A->nowhead,-1);

	done = (int *) malloc(nag*sizeof(int));
	memset(done,0,nag*sizeof(int));

	memset(fn,0,128*sizeof(char));
	sprintf(fn,"popdy%03d.dat",A->r);
	fp = fopen(fn,"a");



	do{
		i = 0;
		found=0;
		finished = 1;
		for(i=0,pa=A->nowhead;i<nag;i++,pa=pa->next){
			if(!done[i]){
				if(!found){
					done[i]=1;
					count=1;
					finished=0;
					found=1;
					spc = pa->spp->spp;
				}
				else{
					if(pa->spp->spp==spc){
						done[i]=1;
						count++;
					}
				}
			}
		}

		//Write to file
		if(!finished)
			fprintf(fp,"%d,%d,%d\n",t,spc,count);

	}while(!finished);

	fflush(fp);
	fclose(fp);
	free(done);
}




int run_one_comass_trial(const int rr, stringPM *A,  int * params, struct runparams *R){


	FILE *mc;char fn[256],pfn[256];FILE *ftmp;
	//sprintf(fn,"maxcodes%03d.txt",rr);
	//mc = fopen(fn,"w");
	//fclose(mc);


	int *maxcode;
	//todo: what is the relationship between maxcode and params...?
	//if(!rr)
	maxcode = (int *) malloc(A->blosum->N * sizeof(int));
	memset(maxcode,0,A->blosum->N*sizeof(int));


	A->r=rr;
	sprintf(pfn,"popdy%03d.dat",A->r);
	ftmp = fopen(pfn,"w");
	fclose(ftmp);

	//todo DELETE if we don't need epochs any more
	//int lastepoch=A->get_ecosystem(),thisepoch,nepochs=1;

	A->domut=0;
	int nsteps=0;
	int i;
	for(i=0;R->indefinite || nsteps <= R->maxnsteps;i++){

		A->extit = i;

		A->comass_make_next();
		A->update();


		if(!(i%1000)){
			record_spp(A);

			//todo: put the below in a function
			printf("%03d At  time %d e=%d\n",rr,i,(int)A->energy);
			printsppct(A,i);

			setmaxcode(A,maxcode);

			/*
			printf("CODE\tMAX\tMASS\tMAX_USED\n");
			for(int k=0;k<A->blosum->N;k++){
				printf("%c:\t%d\t%d\t%d",A->blosum->key[k],params[k],A->mass[k],maxcode[k]);
				if(A->mass[k]<0 || A->mass[k]>params[k])
					printf(" ERROR\n");
				else
					printf("\n");

			}
			printf("\n");
			*/
		}

		if(!A->nagents(A->nowhead,-1)){// || (!rr && nsteps >1500000) || nsteps >15000000){
			printf("DEATH\n");
			printf("At  time %d e=%d\t",i,(int)A->energy);
			A->print_spp_count(stdout,0,-1);
			nsteps=i;
			break;
		}

		nsteps++;
		A->energy += 20;
	}

	free(maxcode);
	return i;
}


void setmaxcode(stringPM *A, int *maxcode){

	s_ag *pag;


	for(int i=0;i<A->blosum->N;i++){
		int count=0;
		for(pag=A->nowhead;pag!=NULL;pag=pag->next){
			int C = strlen(pag->S);
			for(int c=0; c<C; c++){
				if(pag->S[c]==A->blosum->key[i]){
					count++;
				}
			}
		}
		maxcode[i]=count<maxcode[i]?maxcode[i]:count;
	}
}

//These are the functions we need to manipulate mutation networks:
//int ** random_mtx(const int N){
//	int ** matrix;
//
//
//}


int make_mtx_file(char *fn, char *basisfn, int **mut, const int N){
	FILE *in,*out;
	char line[2048];
	int ncodes;
	//int **mut;
	in = fopen(basisfn,"r");
	out = fopen(fn,"w");

	//Get the 1st line - it contains the codes:
	fscanf(in,line);
	ncodes = strlen(line);
	fprintf(out,"%s\n",line);

	//Copy the sw table
	for(int i=0;i<ncodes+1;i++){
		fscanf(in,line);
		fprintf(out,"%s\n",line);
	}

	//print the mutation network

}

void setmutnet(int * mutnet, swt *blosum){

	int i,j;
	for(i=0;i<blosum->N;i++)
		for(j=0;j<blosum->N;j++){
			blosum->adj[i][j] = mutnet[(i*blosum->N)+j];

		}
}



float **swdt(stringPM * A, stringPM * B,s_sw **spp_matches, int *mol_class, float *class_score, float *self){

	int a,b,
		nA,//=A->spp_count(),
		nB;//=B->spp_count();


	A->get_spp_count(-1);
	A->count_spp();
	nA = A->spl->spp_count - 1;

	B->get_spp_count(-1);
	B->count_spp();
	nB = B->spl->spp_count - 1;

	float **dt;
	l_spp *pA,*pB;
	align sw;
	s_sw * swa;
	int L,lA,lB;

	dt = (float **) malloc(nA * sizeof(float *));

	for(a=0,pA=A->spl->species;a<nA;a++,pA=pA->next){
		dt[a]= (float *) malloc(nB * sizeof(float));
		pB = B->spl->species;
		lA = strlen(pA->S);
		for(b=0,pB=B->spl->species;b<nB;b++,pB=pB->next){
			if(pA->count && pB->count){

				int foundsw=0;
				swa = read_sw(*spp_matches,pA->spp,pB->spp);
				if(swa==NULL){
					swa = read_sw(*spp_matches,pB->spp,pA->spp);
					if(swa==NULL){
						float	bprob = SmithWatermanV2(pA->S,pB->S,&sw,A->blosum,0);

						lB = strlen(pB->S);
						L = lA>lB?lA:lB;
						sw.score = sw.score/L;
						dt[a][b]=sw.score/self[a];

						store_sw(spp_matches,&sw,pA->spp,pB->spp);
					}
					else foundsw=1;
				}
				else foundsw=1;
				if(foundsw){
					load_sw(swa,&sw);
					dt[a][b]=sw.score/self[a];
				}
			}
			else{
				dt[a][b]=-1.0;
			}
		}
	}

	//Now we can find the gene/machine ratio - assuming 1st 3 spp are
	//int mol_class{0,0,0,1,1,1};
	//float class_score{0.0,0.0};


	//For each species present:
	for(b=0,pB=B->spl->species;b<nB;b++,pB=pB->next){
		if(pB->count){
			for(a=0,pA=A->spl->species;a<nA;a++,pA=pA->next){
				if(dt[a][b] > -0.5){
					class_score[mol_class[a]] += pB->count * dt[a][b];
				}
			}
		}
	}

	return dt;
}



void print_class_scores(FILE *fp, float *scores,  int N){

	float sum=0.;
	int i;

	for(i=0;i<N;i++){
		sum += scores[i];
	}

	fprintf(fp,"Ratio: ");
	for(i=0;i<N;i++){
		if(i)
			fprintf(fp,":");
		fprintf(fp,"%f",scores[i]/sum);
	}

	fprintf(fp,",\tTotals: ");
	for(i=0;i<N;i++){
		if(i)
			fprintf(fp,":");
		fprintf(fp,"%f",scores[i]);
	}

	fprintf(fp,"\n");
}

void print_swdt(FILE *fp, float **dt,  stringPM * A, stringPM * B){

	int a,b,
		nA,//=A->spp_count(),
		nB;//=B->spp_count();

	A->count_spp();
	nA = A->spl->spp_count - 1;

	B->count_spp();
	nB = B->spl->spp_count - 1;

	l_spp *pA,*pB;

	int *found;
	found = (int *) malloc(nB*sizeof(int));
	memset(found,0,nB*sizeof(int));

	//Print the column names
	printf("Species distance table:\n");
	//fprintf(fp,"\t");
	//for(b=0,pB=B->spl->species;b<nB;b++,pB=pB->next){
	for(b=0,pB=B->spl->species;b<nB;b++,pB=pB->next){
		if(pB->count){
			found[b]=1;
			fprintf(fp,"\t%03d",pB->spp);
		}
	}
	fprintf(fp,"\n");
	for(a=0,pA=A->spl->species;a<nA;a++,pA=pA->next){
		fprintf(fp,"%03d",pA->spp);
		for(b=0;b<nB;b++){

			if(found[b]){
				fprintf(fp,"\t%0.3f",dt[a][b]);
			}
		}
		fprintf(fp,"\n");
	}
}



//get stats for evolution cf seed community
void evostats(char * Afn, stringPM *B,s_sw **spp_matches, float *self, float *gvm){

	//Create a copy of the seed replicase population

	SMspp		SP;
	stringPM A(&SP);
	A.load(Afn,NULL,0,0);

	//TODO: read this data from the config file somehow...
	const int nseedspp = 6;
	int *mol_class;//{0,0,0,1,1,1};
	float *class_score;//{0.0,0.0};
	mol_class = (int *)malloc(nseedspp * sizeof(int));
	class_score = (float *)malloc(2 *sizeof(float));

	memset(class_score,0,2 *sizeof(float));
	for(int i=0;i<nseedspp;i++){
		if(i<3)
			mol_class[i]=0;
		else
			mol_class[i]=1;
	}

	//Get min dist from spp to spp
	float ** dt;
	dt = swdt(&A,B,spp_matches,mol_class,class_score,self);

	//print_swdt(stdout,dt,&A,B);

	//print_class_scores(stdout,class_score,2);

	//printf("Species mapping in B is:\n");

	//TODO: cleanup properly!
	free(dt);
	free(mol_class);
	free(class_score);


}



float * self_stats(char * Afn){

	//Create a copy of the seed replicase population

	SMspp		SP;
	stringPM A(&SP);
	A.load(Afn,NULL,0,0);

	//Get min dist from spp to spp

	A.get_spp_count(-1);
	A.count_spp();
	int nA = A.spl->spp_count - 1;
	l_spp *pA;

	float *self;
	self = (float *)malloc(nA*sizeof(float));
	pA=A.spl->species;
	for(int a=0;a<nA;a++,pA=pA->next){

		align sw;
		SmithWatermanV2(pA->S,pA->S,&sw,A.blosum,0);

		int L = strlen(pA->S);
		sw.score = sw.score/L;
		self[a]=sw.score;

	}

	return self;
}


int arg_load(stringPM *A, int argc, char *argv[], int verbose ){


	switch(argc){
	case 3:
		if(verbose)printf("Traditional config\n");
		A->load(argv[2],NULL,0,1);
		return 1;
	case 4:
		if(verbose)printf("youShare-compatible config\n");
		A->load(argv[2],argv[3],0,1);
		return 1;
	default:
		if(verbose)printf("Error: wrong number of arguments - try 2 or 3\n");
		return 0;
	}


}

void print_params(stringPM *A, int ntrials, int nsteps){

	printf("Non-stringPM variables:\n");
	if(ntrials<0)
		printf("NTRIALS     not set - the default value would be used if needed\n");
	else
		printf("NTRIALS     %d\n",ntrials);

	if(nsteps<0)
		printf("NSTEPS      not set - the default value would be used if needed\n");
	else
		printf("NSTEPS      %d\n",nsteps);

	//load params:
	printf("CELLRAD     %f (vcellrad = %f)\n",A->cellrad,A->vcellrad);
	printf("AGRAD       %f\n",A->agrad);
	printf("ENERGY      %d\n",(int) A->energy);
	printf("NSTEPS      %f\n",A->nsteps);

	/** Although agents_base contains a call to load_influx, Stringmol never uses it */

	//load agents: load_table(_mtx); load_mut; load_decay
	if(A->blosum == NULL)
		printf("BLOSUM      not set - needs to be loaded explicitly\n");
	else
		printf("BLOSUM      %d size table loaded\n",A->blosum->N);
	printf("MUTATE      indelrate = %f; subrate = %f\n",A->indelrate,A->subrate);
	printf("DECAY       %f\n",A->decayrate);
	printf("MAXLEN      %d, (maxl0 = %d)\n",A->maxl, A->maxl0);
	printf("ESTEP       %d\n",A->estep);


}

void check_config( int argc, char *argv[]){
	/** The idea here is to report the default values of the parameters, then parse the config and report them again. */
	SMspp		SP;
	stringPM A(&SP);
	int ntrials = -1;
	int nsteps = -1;

	printf("\nBEFORE loading the config, params are:\n");
	print_params(&A,ntrials,nsteps);
	printf("..c'est ca!\n\n");

	//int readordef_param_int(char *fn, const char *label, int *val, const int defaultvalue, const int verbose)
	readordef_param_int(argv[2], "NTRIALS", &ntrials, 1, 1);
	int nns = readordef_param_int(argv[2], "NSTEPS", &nsteps, -1, 1);

	if(!arg_load(&A, argc, argv, 0))
		return;

	printf("\n\nAFTER loading the config, params are:\n");
	print_params(&A,ntrials,nsteps);
	if(nns==1)
		printf("NSTEPS was not specified. Simulations will run indefinitely");
	printf("..c'est ca!\n\n");
}

/* This is used in the comass GA - see test.cpp for how to set up randseed properly
 *
 */
void init_randseed_config(int argc, char *argv[]){

	/*TODO: set up random seed properly - the best result was *without* a random seed, and has been lost */
	/* Funny business with unsigned longs... */
	//printf("%u\n", -873302838);//-2041524348);
	char test[128];
	memset(test,0,128*sizeof(char));
	sprintf(test,"-873302838");
	long dummy;
	sscanf(test,"%d",&dummy);
	//printf("Random seed is %d\n",dummy);
	//printf("Random seed is %lu (unsigned)\n",dummy);


	unsigned long seedin =  2846144656u;

	int qnnscoring = 1;

	FILE *fpr;
	if((fpr=fopen(argv[2],"r"))!=NULL){
		int stmp;
		int rerr = read_param_int(fpr,"RANDSEED",&stmp,1);
		rerr = read_param_int(fpr,"GAQNN",&qnnscoring,1);
		if(rerr)
			qnnscoring=1;
		seedin = stmp;
		fclose(fpr);
	}

	unsigned long rseed = longinitmyrand(&seedin);//437);//-1);//437);
	//unsigned long rseed = longinitmyrand(NULL);//437);//-1);//437);
	FILE *frs;
	if((frs=fopen("randseed.txt","w"))==NULL){
		printf("Coundln't open randseed.txt\n");
		getchar();
	}
	fprintf(frs,"(unsigned) random seed is %lu \n",rseed);
	fflush(frs);
	fclose(frs);

}








