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
#include "instructions.h"

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
void setupSMol(struct runparams &RunPar, int argc, char *argv[]){

	FILE *fp;

	unsigned int rerr=1,rlim=20;
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
			printf("NTRIALS not specified;\nSetting NTRIALS to %d\n",rlim);
			break;
		}
		fclose(fp);
	}

	//Read nsteps:
	if((fp=fopen(argv[2],"r"))!=NULL){
		rerr = read_param_int(fp,"NSTEPS",&(RunPar.maxnsteps),1);
		switch(rerr){
		case 2:
			printf("Multiple NSTEPS specified. Check config file\n");
			getchar();
			exit(0);
		case 0:
			printf("Setting NSTEPS to %d\n",RunPar.maxnsteps);
			RunPar.indefinite=0; //TODO: fix the indefinite thing if NSTEPS is not specified..
			break;
		default:
			printf("NSTEPS not specified;\nEach trial will run to extinction\n");
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

	char pfn[256];FILE *ftmp;

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
	unsigned int nsteps=0;
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
						//TODO: this returns a float - we don't use it...
						SmithWatermanV2(pA->S,pB->S,&sw,A->blosum,0);

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

	//TODO: This can be a bit fragile for runs with >2 arguments... be careful!
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
	sscanf(test,"%ld",&dummy);
	//printf("Random seed is %d\n",dummy);
	//printf("Random seed is %lu (unsigned)\n",dummy);


	unsigned long seedin =  2846144656u;

	unsigned int qnnscoring = 1;

	FILE *fpr;
	if((fpr=fopen(argv[2],"r"))!=NULL){
		unsigned int stmp;
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








/*
 * NB: To get identical trials to those run for ALifeXII, do the following:
 * 1: Use the file "replicase.conf" as the input
 * 2: Fix the random number seed to 437
 * 3: #define DO_ANCESTRY to get ancestry files out...
 * 4: Run on a 32-bit linux slackware system, circa 2008 vintage...:)
 */
//todo: SmPm_AlifeXII() should call this
int run_one_AlifeXII_trial(stringPM *A){

	int i;



	A->print_agents(stdout,"NOW",0);

	A->r=0;
	//sprintf(pfn,"popdy%d%02d.dat",proc,A.r);
	//ftmp = fopen(pfn,"w");
	//fclose(ftmp);

#ifdef DO_ANCESTRY
	int lastepoch=A.get_ecosystem(),thisepoch,nepochs=1;
#endif

	int nsteps=0;
	//TODO: Accommodate indefinitre runs, like this:
	//for(i=0;indefinite || nsteps <= maxnsteps;i++){

	for(i=0;nsteps <= A->nsteps;i++){

		//TODO: find out what this does - rename the variable to  make it clear.
		A->extit = i;

		A->make_next();
		A->update();

		if(!(i%1000)){
			A->print_spp_count(stdout,0,-1);
		}

		if(!(i%1000)){
			printf("At  time %d e=%d, mutrate = %0.9f & %0.9f\n",i,(int)A->energy,A->subrate,A->indelrate);
			printsppct(A,i);
		}

//TODO: See equivalent line in SmPm_AlifeXII() function for what should be in the following #ifdef...
//#ifdef DO_ANCESTRY
//#endif

		if(!A->nagents(A->nowhead,-1)){
			printf("DEATH\n");
			printf("At  time %d e=%d, mutrate = %0.9f & %0.9f\t",i,(int)A->energy,A->indelrate,A->subrate);
			A->print_spp_count(stdout,0,-1);
			nsteps=i;
			break;
		}
		nsteps++;

		A->energy += A->estep;
	}

	printf("Finished - alls well!\nclear out memory now:\n");
	fflush(stdout);

	//TODO: need to do this outside the function!
	//A.clearout();

	return 0;
}

//////////////////////////////////////////////////////////////////////////
// SPATIAL STRINGMOL CODE - CALLED FROM stringmol.cpp AND smspatial.cpp //
//////////////////////////////////////////////////////////////////////////



int randy_Moore(const int X, const int Y, const int Xlim, const int Ylim, int *xout, int *yout){

	int pos = 8. * rand0to1();

	/*    012
	 *    3*4
	 *    567
	 */

	int xoff = 0;
	int yoff = 0;

	switch(pos){
	case 0:
		xoff = -1;
		yoff = -1;
		break;
	case 1:
		yoff = -1;
		break;
	case 2:
		xoff = 1;
		yoff = -1;
		break;
	case 3:
		xoff = -1;
		break;
	case 4:
		xoff = 1;
		break;
	case 5:
		xoff = -1;
		yoff = 1;
		break;
	case 6:
		yoff = 1;
		break;
	case 7:
		xoff = 1;
		yoff = 1;
		break;
	}

	*xout = (X+Xlim+xoff)%Xlim;
	*yout = (Y+Ylim+yoff)%Ylim;
	return 0;
}





smsprun * init_smprun(const int gridx, const int gridy){

	smsprun *run;
	run = (smsprun *) malloc(sizeof(smsprun));

	//TODO: make grid size changeable via config...
	run->gridx=gridx;
	run->gridy=gridy;

	run->grid=(s_ag ***) malloc(run->gridx*sizeof(s_ag **));
	run->status=(s_gstatus **) malloc(run->gridx*sizeof(s_gstatus *));

	for(int i=0;i<run->gridx;i++){
		run->grid[i] = (s_ag **) malloc(run->gridy*sizeof(s_ag *));
		run->status[i] = (s_gstatus *) malloc(run->gridy*sizeof(s_gstatus));

		for(int j=0;j<run->gridy;j++){
			run->grid[i][j]=NULL;
			run->status[i][j]=G_EMPTY;
		}
	}

	return run;
}

void obsolete_find_ag_gridpos(s_ag *pag,smsprun *run, int *x, int *y){
	for(int i=0;i<run->gridx;i++){
		for(int j=0;j<run->gridy;j++){
			if( run->grid[i][j] == pag ){
				*x=i;
				*y=j;
				return;
			}
		}
	}
	printf("Unable to find agent %d\n",pag);
}


s_ag * pick_partner(stringPM *A, smsprun *run,int x, int y){

	int i,j,xx,yy;

	//first, let's count the agents
	int count = 0;
	for(i=-1;i<2;i++){
		for(j=-1;j<2;j++){
			if( !(i == 0 && j ==0) ){
				xx = (x + i + run->gridx)%run->gridx;
				yy = (y + j + run->gridy)%run->gridy;
				if(run->grid[xx][yy]!=NULL){
					if(run->grid[xx][yy]->status == B_UNBOUND){
						if(run->status[xx][yy] == G_NOW){
							count++;
						}
					}
				}
			}
		}
	}
	if(!count)
		return NULL;

	int it = count * rand0to1();


	//now, let's choose the agents
	count = 0;
	for(i=-1;i<2;i++){
		for(j=-1;j<2;j++){
			if( !(i == 0 && j ==0) ){
				xx = (x + i + run->gridx)%run->gridx;
				yy = (y + j + run->gridy)%run->gridy;
				if(run->grid[xx][yy]!=NULL){
					if(run->grid[xx][yy]->status == B_UNBOUND){
						if(run->status[xx][yy] == G_NOW){
							if(count==it)
								run->status[xx][yy] = G_NEXT;
								return run->grid[xx][yy];
							count++;
						}
					}
				}
			}
		}
	}

	//We should never get to here
	printf("Something's wrong - neighbour detected in first pass but none selected\n");
	return NULL;
}


void update_grid(smsprun *run){
	for(int i=0;i<run->gridx;i++){
		for(int j = 0;j<run->gridy;j++){
			if(run->grid[i][j]!=NULL){
				if(run->status[i][j]==G_NEXT)
					run->status[i][j] = G_NOW;
			}
			else{
				run->status[i][j] = G_EMPTY;
			}
		}
	}
}

void place_mol(s_ag *ag,smsprun *run,int x, int y){
	//TODO: error checking!
	run->grid[x][y] = ag;
	//Set status to next - placing the molecule has used up the 'now'
	run->status[x][y] = G_NEXT;
	ag->set=true;
	ag->x=x;
	ag->y=y;
}


s_ag * place_neighbor(stringPM *A, smsprun *run,s_ag *c,int x,int y){


	int xx,yy;
	int nvacant=0;
	for(int i=-1;i<2;i++){
		for(int j=-1;j<2;j++){
			xx = (x+i+run->gridx)%run->gridx;
			yy = (y+j+run->gridy)%run->gridy;
			//No need to ingnore x,y because it is occupied by the parent!
			if(run->grid[xx][yy]==NULL)
				nvacant++;
		}
	}
	if(nvacant){
		//Decide where to put it:
		int pos = nvacant * rand0to1();
		int here=0;
		for(int i=-1;i<2;i++){
			for(int j=-1;j<2;j++){
				xx = (x+i+run->gridx)%run->gridx;
				yy = (y+j+run->gridy)%run->gridy;
				//No need to ingnore x,y because it is occupied by the parent!
				if(run->grid[xx][yy]==NULL){
					if(here == pos){
						place_mol(c,run,xx,yy);
						return c;
					}
					here++;
				}
			}
		}
	}
	else{
		return NULL;
	}
	/*TODO: Need to decide whether to replace or not
	else{
		int noccupied = 8-nvacant;
		if(run->grid[x][y]==NULL){
			//This should never happen....
			printf("Alert! empty parent cell!\n");
			noccupied++;
		}

		int pos = nvacant * rand0to1();
		int here=0;
		for(int i=-1;i<2;i++){
			for(int j=-1;j<2;j++){
				xx = (x+i+run->gridx)%run->gridx;
				yy = (y+j+run->gridy)%run->gridy;
				//No need to ingnore x,y because it is occupied by the parent!
				if(i!=0 && j!=0){
					if(run->grid[xx][yy]!=NULL){
						if(here == pos){
							//Remove the incumbent


							run->grid[xx][yy]=c;
							return;
						}
						here++;
					}
				}
			}
		}
	}
	*/
}



int spatial_cleave(stringPM *A, smsprun *run, s_ag *act){//, int x, int y){

	int dac = 0,cpy;
	s_ag *c,*pass,*csite;
	c = NULL;
	pass = act->pass;

	//pick the mol containing the cleave site:
	csite = act->ft?act:pass;

	if(act->f[act->ft]-csite->S < csite->len){

		//1: MAKE THE NEW MOLECULE FROM THE CLEAVE POINT
		c = A->make_ag(pass->label);//,1);

		//Copy the cleaved string to the agent
		char *cs;
		c->S =(char *) malloc(A->maxl0*sizeof(char));
		memset(c->S,0,A->maxl0*sizeof(char));
		cs = csite->S;
		cpy = strlen(cs);

		//Check that we aren't creating a zero-length molecule:
		if(!cpy){
			printf("WARNING: Zero length molecule being created!\n");
		}

		cpy -= act->f[act->ft]-cs;

		if(!cpy){
			printf("ERROR: Zero length molecule definitely being created!\nbail..\n");
			A->free_ag(c);
		}
		else{

			strncpy(c->S,act->f[act->ft],cpy);
			c->len = strlen(c->S);

#ifdef VERBOSE
		printf("String %d created:\n%s\n",c->idx,c->S);
#endif

			//Check the lineage
			A->update_lineage(c,'C',1,act->spp,pass->spp,act->biomass);
			act->biomass=0; //reset this; we might continue to make stuff!

			//TODO: place the new agent on the grid
			if((place_neighbor(A,run,c,act->x,act->y))!=NULL){//,x,y))!=NULL){
				//append the agent to nexthead
				A->append_ag(&(A->nexthead),c);
			}
			else{
				A->free_ag(c);
			}
		}
		//TODO: check string lens of act and pass?


		//2: HEAL THE PARENT
		memset(act->f[act->ft],0,cpy*sizeof(char));

		csite->len = strlen(csite->S);
		if(csite->len==0){
			printf("zero length parent string!\n");
		}


		//Get rid of zero-length strings...
		//NB - grid status will be updated at the end of the timestep - simpler.
		if((dac = A->check_ptrs(act))){
			int x,y;
			switch(dac){
			case 1://Destroy active - only append passive
				A->unbind_ag(pass,'P',1,act->spp,pass->spp);
				A->append_ag(&(A->nexthead),pass);
				//find_ag_gridpos(pass,run,&x,&y);
				//run->status[x][y]=G_NEXT;
				run->status[pass->x][pass->y]=G_NEXT;

				//find_ag_gridpos(act,run,&x,&y);
				run->grid[act->x][act->y]=NULL;
				run->status[act->x][act->y]=G_EMPTY;

				A->free_ag(act);

				break;
			case 2://Destroy passive - only append active
				A->unbind_ag(act,'A',1,act->spp,pass->spp);
				A->append_ag(&(A->nexthead),act);
				//find_ag_gridpos(act,run,&x,&y);
				run->status[act->x][act->y]=G_NEXT;


				//find_ag_gridpos(pass,run,&x,&y);
				run->grid[pass->x][pass->y]=NULL;
				run->status[pass->x][pass->y]=G_EMPTY;

				A->free_ag(pass);
				break;
			case 3://Destroy both
				printf("Destroying both parents after cleave!\nThis should never happen!\n");
				A->unbind_ag(act,'A',1,act->spp,pass->spp);
				A->unbind_ag(pass,'P',1,act->spp,pass->spp);
				A->free_ag(act);
				A->free_ag(pass);
				break;
			default://This can't be right can it?
				if(act->ft == act->it){
					act->i[act->it]--;
				}
				break;
			}
		}
	}
	if(!dac){
		act->i[act->it]++;
	}

	return dac;
}






int spatial_exec_step(stringPM *A, smsprun *run, s_ag *act, s_ag *pass){//, int x, int y){

	char *tmp;
	int dac=0;
	int safe_append=1;

	switch(*(act->i[act->it])){//*iptr[it]){

	case '$'://h-search
		//act->ft = act->it;
		char *cs;
		if(act->ft)
			cs = act->S;
		else
			cs = act->pass->S;
		tmp = HSearch(act->i[act->it],cs,A->blosum,&(act->it),&(act->ft),A->maxl);
		act->f[act->ft] = tmp;
		act->i[act->it]++;
		break;

	/*************
	 *   P_MOVE  *
	 *************/
	case '>':
			tmp=act->i[act->it];
			tmp++;
			switch(*tmp){
			case 'A':
				act->it = act->ft;
				act->i[act->it] = act->f[act->ft];
				act->i[act->it]++;
				break;
			case 'B':
				act->rt = act->ft;
				act->r[act->rt] = act->f[act->ft];
				act->i[act->it]++;
				break;
			case 'C':
				act->wt = act->ft;
				act->w[act->wt] = act->f[act->ft];
				act->i[act->it]++;
				break;
			default:
				act->it = act->ft;
				act->i[act->it] = act->f[act->ft];
				act->i[act->it]++;
				break;
			}
			break;


	/************
	 *   HCOPY  *
	 ************/
	case '='://h-copy
		if(A->hcopy(act)<0){
			A->unbind_ag(act,'A',1,act->spp,pass->spp);
			A->unbind_ag(pass,'P',1,act->spp,pass->spp);
		}
		break;


	/************
	 *   INC_R  *
	 ************/
	case '+'://h-copy
		if(A->granular_1==1){
			//printf("Incrementing read \n");
			/* Select the modifier */
			tmp=act->i[act->it];
			tmp++;
			switch(*tmp){
			case 'A':
				act->i[act->it]++;
				break;
			case 'B':
				act->r[act->rt]++;
				break;
			case 'C':
				act->w[act->wt]++;
				break;
			default:
				act->f[act->ft]++;
				break;
			}
		}
		act->i[act->it]++;
		break;



	/************
	 *  TOGGLE  *
	 ************/
	case '^'://p-toggle: toggle active pointer
			tmp=act->i[act->it];
			tmp++;
			switch(*tmp){
			case 'A':
				act->it = 1-act->it;
				break;
			case 'B':
				act->rt = 1-act->rt;
				break;
			case 'C':
				act->wt = 1-act->wt;
				break;
			default:
				act->ft = 1-act->ft;
				break;
			}
			act->i[act->it]++;
			break;

	/************
	 *  IFLABEL *
	 ************/
	case '?'://If-label
			act->i[act->it]=IfLabel(act->i[act->it],act->r[act->rt],act->S,A->blosum,A->maxl);
			break;


	/************
	 *  CLEAVE  *
	 ************/
	case '%':
			//Decide where to put the cleaved molecule
			if((dac = spatial_cleave(A,run,act))){//,x,y))){
				//TODO: Need to determine what safe_append is used for (after looking at cleave)
				safe_append=0;	//extract_ag(&nowhead,p);
			}
			break;

	/**************
	 *  TERMINATE *
	 **************/
	case 0:
	case '}'://ex-end - finish execution

#ifdef V_VERBOSE
			printf("Unbinding...\n");
#endif

			A->unbind_ag(act,'A',1,act->spp,pass->spp);
			run->status[act->x][act->y] = G_NEXT;

			A->unbind_ag(pass,'P',1,act->spp,pass->spp);
			//int xx,yy;
			//find_ag_gridpos(pass,run,&xx,&yy);
			run->status[pass->x][pass->y] = G_NEXT;

			break;

	default://Just increment the i-pointer
		act->i[act->it]++;
		break;
	}
#ifdef V_VERBOSE
	printf("Exec step - looks like:\n");
	print_exec(stdout,act,pass);
#endif


	//TODO: This action should be elsewhere - much harder to follow here
	if(safe_append){
		act->ect++;
		A->append_ag(&(A->nexthead),act);
		A->append_ag(&(A->nexthead),pass);
	}
	A->energy--;


	return 1;

}

int spatial_testdecay(stringPM *A, smsprun *run, s_ag *pag){

 	float prob = A->decayrate;//1./pow(65,2);//4./3.); //This is now done in load_decay...

	float rno = rand0to1();

	if(rno<prob){
		//unbind_ag(pag);


		s_ag *bag;
		switch(pag->status){
		case B_UNBOUND:
			bag = NULL;
			break;
		case B_ACTIVE:
			bag = pag->pass;
			break;
		case B_PASSIVE:
			bag = pag->exec;
			break;
		}
		int x,y;

		//find_ag_gridpos(pag,run,&x,&y);
		run->grid[pag->x][pag->y]=NULL;
		run->status[pag->x][pag->y]=G_EMPTY;

		A->free_ag(pag);

		if(bag!=NULL){

			//find_ag_gridpos(bag,run,&x,&y);
			run->grid[bag->x][bag->y]=NULL;
			run->status[bag->x][bag->y]=G_EMPTY;

			A->free_ag(bag);
		}

		return 1;
	}
	else
		return 0;
}


int smspatial_step(stringPM *A, smsprun *run){

	//again, we follow make_next, but are a little more careful with the binding and uncoupling
	s_ag *pag,*bag;
	int changed;

	A->energy += A->estep;

	while(A->nowhead!=NULL){

		pag = A->rand_ag(A->nowhead,-1);
		A->extract_ag(&A->nowhead,pag);
		changed = 0;

		//extract any partner:
		bag = NULL;

		switch(pag->status){
		case B_UNBOUND:
			break;
		case B_ACTIVE:
			bag = pag->pass;
			A->extract_ag(&(A->nowhead),bag);
			break;
		case B_PASSIVE:
			bag = pag->exec;
			A->extract_ag(&(A->nowhead),bag);
			break;
		}

		if(!spatial_testdecay(A,run,pag)){
			if(A->energy>0){
				switch(pag->status){
				case B_UNBOUND:
					//seek binding partner, set binding states.
					//changed = A->testbind(pag);

					//int x,y;
					align sw;

					//We can only bind neighbours in the spatial model
					//find_ag_gridpos(pag,run,&x,&y);
					run->status[pag->x][pag->y]=G_NEXT;

					if((bag = pick_partner(A,run,pag->x,pag->y))!=NULL){

						//TODO: We need to make sure that bag is in nowhead first!
						A->extract_ag(&(A->nowhead),bag);

						//Now we've found a potential partner, we can see if it binds:
						float bprob;
						bprob = A->get_sw(pag,bag,&sw);

						float rno;
						rno = rand0to1();
						if(rno<bprob){//Binding success!
							//figure out which is the executing string:
							A->set_exec(pag,bag,&sw);
							pag->nbind++;
							bag->nbind++;

							A->energy--;

							A->append_ag(&(A->nexthead),pag);
							A->append_ag(&(A->nexthead),bag);
							changed=1;
						}
					}

					break;
				case B_PASSIVE:


					//find_ag_gridpos(pag->exec,run,&x,&y);

					changed = spatial_exec_step(A,run,pag->exec,pag);//,x,y);

					break;
				case B_ACTIVE:

					//find_ag_gridpos(pag,run,&x,&y);
					changed = spatial_exec_step(A,run,pag,pag->pass);//,x,y);
					break;
				default:
					printf("ERROR: agent with unknown state encountered!\n");
				}
			}
			if(!changed){
				A->append_ag(&(A->nexthead),pag);
				if(bag!=NULL)
					A->append_ag(&(A->nexthead),bag);

			}
		}
	}

	A->update();
	update_grid(run);

	return 0;
}




int smspatial_init(char *fn, stringPM *A, smsprun **run){


	A->load(fn,NULL,0,1);

	*run = init_smprun(300,300);

	//Now we have to place each agent on the grid - use the makenext() model -
	while(A->nowhead!=NULL){
		s_ag *pag;
		pag = A->rand_ag(A->nowhead,-1);
		A->extract_ag(&(A->nowhead),pag);
    	int found = 0;
    	while(!found){
			int pos = (*run)->gridx * (*run)->gridy * rand0to1();

			int x = pos/(*run)->gridx;
			int y = pos%(*run)->gridx;

			if((*run)->grid[x][y]==0){

				//TODO: this command should be moved to the smspatial
		        //((uint8_t *)screen->pixels)[x + (y * sdlPitch)] = 0;

		        //Add the partner in the Moore neighborhood
		        int ffound=0;

		        while(!ffound){
					int xx,yy;
					//randy_Moore(const int X, const int Y, const int Xlim, const int Ylim, int *xout, int *yout){
					randy_Moore(x,y,(*run)->gridx,(*run)->gridy,&xx,&yy);
					if((*run)->grid[xx][yy]==0){

						ffound=found=1;

						//Place each agent on the list
						place_mol(pag,*run,x,y);

						//Move to the 'used' bucket
						A->append_ag(&(A->nexthead),pag);

						s_ag *bag;
						bag = A->rand_ag(A->nowhead,-1);
						if(bag != NULL){
							A->extract_ag(&(A->nowhead),bag);
							place_mol(bag,*run,xx,yy);
							A->append_ag(&(A->nexthead),bag);
						}
					}
		        }
			}
    	}
    }

	A->update();

	return 0;
}

int smspatial(int argc, char *argv[]) {

	printf("Hello spatial stringmol world\n");

	SMspp		SP;
	stringPM	A(&SP);

	smsprun *run;
	run = NULL;

	smspatial_init(argv[2],&A,&run);

	int bt=0,ct=0;
	ct = A.nagents(A.nowhead,-1);
	printf("Initialisation done, number of molecules is %d\n",ct);

	int iteration = 0;
//	while(A.nagents(A.nowhead,-1)){
	while(iteration < 100000){
		smspatial_step(&A,run);
		ct = A.nagents(A.nowhead,-1);
		bt = ct - A.nagents(A.nowhead,B_UNBOUND);
#ifdef DODEBUG
		printf("Nowhead is %d, Nexthead is %d\n",A.nowhead,A.nexthead);
		s_ag *p;
		p=A.nowhead;
		int mno=0;
		while(p!=NULL){
			int x,y;
			find_ag_gridpos(p,run,&x,&y);
			printf("%d, %d, [%d,%d]  status: %d, bound to %d / %d, prev = %d, next = %d\n",++mno,p,x,y,p->status,p->exec,p->pass,p->prev,p->next);
			p = p->next;
		}
#endif
		iteration++;
		if(!(iteration%100))
				printf("Step %d done, number of molecules is %d, nbound = %d\n",iteration,ct,bt);
//		if(iteration == 66396)
//			printf("Pauuuuse\n!");
		if((!(iteration%10000)) ||     iteration == 66396    ){
			FILE *fp;char fn[128];
			sprintf(fn,"splist%d.dat",iteration);
			fp = fopen(fn,"w");
			SP.print_spp_list(fp);
			fclose(fp);
		}
	}

	printf("FINISHED smspatial\n");
	fflush(stdout);
	return 0;
}
