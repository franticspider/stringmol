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


/* TODO: Description */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//utils
//extern "C" {
#include "randutil.h"
#include "params.h"
#include "hsort.h"
//}

//stringmol
#include "alignment.h"

//metabolism
#include "rules.h"
#include "agents_base.h"
#include "SMspp.h"
#include "stringPM.h"

//signal
#include "signalSM.h"

//TODO: Sort out this dependency nightmare!
// Writing PNGs
#include "lodepng.h"
#include <iostream>
//setup
#include "setupSM.h"

//microbial GA
#include "microbial_ga.h"

//tests
#include "tests.h"

//DEFINES
//Calculate the ancestry & epochs inline
//#define DO_ANCESTRY



//FORWARD FUNCTION DECLARATIONS
//extern const int maxl;//150;// 512;
//extern const int maxl0; //allow room for a terminating 0

void printsppct(stringPM *A, int t);


int joinsplists(int argc, char *argv[]){


	FILE *fin;
	int nlists = atoi(argv[3]);
	char fn[256];
	char seq[2001];
	char line[3000];
	SMspp		SP;
	stringPM	oA(&SP),*A;
	s_ag *pag;
	int i,spno,spold;
	int tt,tottime=0;

	A = &oA;


	pag = A->make_ag('A');//,1);
	pag->S =(char *) malloc(A->maxl0*sizeof(char));
//First lets do the loading:
	for(i=0;i<nlists;i++){
		spold=-1;
		sprintf(fn,"%s/splist%03d.dat",argv[2],i);
		if((fin=fopen(fn,"r"))!=NULL){

			while((fgets(line,A->maxl,fin))!=NULL){
				memset(pag->S,0,A->maxl0*sizeof(char));



				//p=strtok(line,",");
				//for(int t=0;t<6;t++){
				//	p = strtok(NULL,",");
				//}
				sscanf(line,"%d,%*d,%*d,%*f,%*d,%*d,%s",&spno,seq);
				if(spno!=spold){
					spold = spno;
					strncpy(pag->S,seq,strlen(seq));
					pag->len = strlen(pag->S);
					printf("adding string %s\n",pag->S);

					//No parents for these initial agents!
					pag->pp = A->spl->make_parents(NULL,NULL);

					//int stringPM::update_lineage(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp)
					A->update_lineage(pag,'I',1,NULL,NULL,0);
					//s = A->spl->getspp(pag,0);
					//	s->tspp = 0;
					//}
				}
			}
			fclose(fin);

			//now get the trial time


		}
		else{
			printf("Unable to open file %s for reading species data\n",fn);
			fflush(stdout);
			getchar();
			return 1;
		}


		//now get the trial time
		sprintf(fn,"%s/popdy%03d.dat",argv[2],i);
		if((fin=fopen(fn,"r"))!=NULL){
			while((fgets(line,A->maxl,fin))!=NULL){
				sscanf(line,"%d",&tt);
			}
			fclose(fin);
		}
		else{
			printf("Unable to open file %s for reading population data\n",fn);
			fflush(stdout);
			getchar();
			return 1;
		}
		tottime += tt;


	}

//second, let's summarize the diversity measures to a file

	printf("nspp = %d, tottime = %d\n",(int) SP.spp_count, tottime);
	return 0;
}


//Called by SmPm_AlifeXII
void check_pop_dy(char *fn,stringPM *A,int it){

	l_spp *s;
	FILE *fout;
	fout=fopen(fn,"a");


	A->get_ecosystem();

	fprintf(fout,"%d :",it);
	s=A->spl->species;
	while(s!=NULL){
		fprintf(fout,"\t%d",s->pf);
		s=s->next;
	}
	fprintf(fout,"\n");
	fclose(fout);
}


void add_spp(const int nag, stringPM *A, char *label, char symbol){
	int i;
	s_ag *pag;

	for(i=0;i<nag;i++){
	l_spp *species;

		pag = A->make_ag(symbol);//,1);

		pag->S =(char *) malloc(A->maxl0*sizeof(char));
		memset(pag->S,0,A->maxl0*sizeof(char));
		strncpy(pag->S,label,strlen(label));
		pag->len = strlen(pag->S);

		//No parents for these initial agents!
		pag->pp = A->spl->make_parents(NULL,NULL);

		if(!i){
			//int stringPM::update_lineage(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp)
			A->update_lineage(pag,'I',1,NULL,NULL,0);
			species = A->spl->getspp(pag,0,A->maxl0);
			species->tspp = 0;
		}
		pag->spp=species;
		if(A->nowhead == NULL){
			A->nowhead = pag;
		}
		else{
			A->append_ag(&(A->nowhead),pag);
		}
	}
}



int biomass_config(stringPM *A, char *fout, const int nbio, const int nag, SMspp *SP, char * frep){

	FILE *out;
	int j;

	//find out ntypes
	//find out len
	//find out pop

	//We need to make an array of species, and an array of strings.
	int *biosum;
	char **str;
	int *index;
	l_spp *ls;
	s_parent *pp;
	char symbol = 'A';
	int count = SP->spp_count;

	biosum = (int *) malloc(count * sizeof(int));
	index = (int *) malloc(count * sizeof(int));
	str = (char **) malloc(count *sizeof(char *));

	memset(biosum,0,count * sizeof(int));
	memset(index,0,count * sizeof(int));

	for(j=0;j<count;j++){
		str[j] = NULL;
	}

	ls = SP->species;
	while(ls != NULL){
		pp = ls->pp;
		while(pp!=NULL){
			if(ls->biomass){
				if(pp->pa != NULL)
					biosum[(pp->pa->spp)] += ls->biomass;
				if(pp->pp != NULL)
					biosum[(pp->pp->spp)] += ls->biomass;
			}
			pp=pp->next;
		}
		str[ls->spp] = (char *) malloc(A->maxl0*sizeof(char));
		memset(str[ls->spp],0,A->maxl0*sizeof(char));
		strncpy(str[ls->spp],ls->S,strlen(ls->S));

		ls = ls->next;
	}

	//indexx_int(count,biosum,index);
	idx_hsort_int(biosum,index,count);

	int sumsum=0,sumlim=0;
	for(j=0;j<count;j++){
		printf("%d, %s\n",biosum[index[j]],str[index[j]]);
		sumsum+=biosum[index[j]];
		if(j<nbio)
			sumlim+=biosum[index[count-j-1]];
	}

	out = fopen(frep,"a");
	fprintf(out,"%d\t%d\n",sumsum,sumlim);
	fclose(out);

	SP->clear_list();

	int fc=0;
	for(j=0;j<nbio;j++){
		//Need to check if *any* biomass has been made by this spp - if not, don't//We actually need to copy the species, cos we are going to clear the list...
		if(biosum[index[count-1-j]]){
			add_spp(nag,A,str[index[count-1-j]],symbol++);
			fc++;
		}
	}

	free(biosum);
	free(index);
	for(j=0;j<count;j++)
		if(str[j]!=NULL)
			free(str[j]);
	free(str);

	out = fopen(fout,"w");
	A->print_conf(out);
	fclose(out);
	return fc;
}



int random_config(stringPM *A, char *fout, const int nnew,const int nag){

	FILE *out;
	int i,j,idx,len=200;//,ntypes=nnew
	char label[A->maxl0];
	float rno;
	s_ag *pag;

	//find out ntypes
	//find out len
	//find out pop

	//write to outfile
	for(j=0;j<nnew;j++){

		memset(label,0,sizeof(char)*A->maxl0);
		for(i=0;i<len;i++){
			rno=rand0to1();
			idx=floor(rno*A->blosum->N);
			label[i]=A->blosum->key[idx];
		}

		for(i=0;i<nag;i++){
			l_spp *s;

			pag = A->make_ag('A'+j);//,1);

			pag->S =(char *) malloc(A->maxl0*sizeof(char));
			memset(pag->S,0,A->maxl0*sizeof(char));
			strncpy(pag->S,label,strlen(label));
			pag->len = strlen(pag->S);

			//No parents for these initial agents!
			pag->pp = A->spl->make_parents(NULL,NULL);

			if(!i){
				//int stringPM::update_lineage(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp)
				A->update_lineage(pag,'I',1,NULL,NULL,0);
				s = A->spl->getspp(pag,0,A->maxl0);
				s->tspp = 0;
			}
			pag->spp=s;
			if(A->nowhead == NULL){
				A->nowhead = pag;
			}
			else{
				A->append_ag(&(A->nowhead),pag);
			}
		}
	}


	//now let's sort out the biomass-based ones we are keeping:



	out = fopen(fout,"w");
	A->print_conf(out);
	fclose(out);
	return 0;
}


/* This function is intended for experiments with random molecules.
 * The idea is that novel mechanisms will arise if we don't code in the answer
 * STATUS: Experimental / Untested / Unpublished
 */
int origlife(int argc, char *argv[]){
	int i,div;

	SMspp		SP;
	stringPM	A(&SP);
	FILE *fp;
	int indefinite=1;

	long rseed = initmyrand(437);//-1);//437);
	//R.printasint();

	unsigned int nsteps = (int) A.nsteps;
	FILE *fsumm,*ftmp;
	char fitfile[128];

	//division values for summary
	//int divit,diven;//divct,

	char	fn[128],pfn[128],randfn[128];

	sprintf(fn,"%s.spatial.summary.dat",argv[1]);
	if((fsumm=fopen(fn,"w"))==NULL){
		printf("Coundlnt open %s\n",fn);
		getchar();
	}
	fprintf(fsumm,"Random seed is %ld\n",rseed);
	fflush(fsumm);


	ftmp = fopen("epochs.dat","w");
	fclose(ftmp);

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
	printf("NTRIALS = %d\n",rlim);



	//Read nsteps:
	//unsigned int maxnsteps=0;
	unsigned int maxnsteps {0};
	if((fp=fopen(argv[2],"r"))!=NULL){
		rerr = read_param_int(fp,"NSTEPS",&maxnsteps,1);
		switch(rerr){
		case 2:
			printf("Multiple NSTEPS specified. Check config file\n");
			getchar();
			exit(0);
		case 0:
			printf("Setting NSTEPS to %d\n",rlim);
			indefinite=0;
			break;
		default:
			printf("NSTEPS not sepcified;\nEach trial will run to extinction\n");
			break;
		}
		fclose(fp);
	}
	printf("NSTEPS = %d\n",maxnsteps);





	for(unsigned int rep=0;rep<rlim;rep++){

		sprintf(fitfile,"rep%03dfitness.dat",rep);
		ftmp=fopen(fitfile,"w");
		fclose(ftmp);
		for(unsigned int rr=0;rr<rlim;rr++){

			div=0;
			//divct = divit = diven = 0;

			//SP.clear_list();

			A.load(argv[2],NULL,0,1);

			test_adj(A.blosum);

			sprintf(randfn,"rep%03dAgents%03d.conf",rep,rr);
			if(rr){
				ftmp=fopen(fitfile,"a");
				fprintf(ftmp,"%d\t",rr-1);
				fclose(ftmp);
				biomass_config(&A, (char *) "bio.conf", 20, 20, &SP, fitfile);
				//possibly put fewer random mols in, or bind sites can "mop up" all successful spp
				random_config(&A, randfn, 30, 10);
			}
			else{
				SP.clear_list();
				random_config(&A, randfn, 30, 20);
			}

			//A.print_agents(stdout,"NOW",0);
			SP.print_spp_list(stdout);

			A.run_number=rr;
			sprintf(pfn,"rep%03dpopdy%03d.dat",rep,A.run_number);
			ftmp = fopen(pfn,"w");
			fclose(ftmp);


			int lastepoch=A.get_ecosystem(),thisepoch,nepochs=1;

			A.domut=1;
			nsteps=0;
			for(i=0;indefinite || nsteps <= maxnsteps;i++){

				A.extit = i;

				A.make_next();
				A.update();

				if(!(i%1000)){
					A.print_spp_count(stdout,0,-1);
					A.extit=i;//dummy line for breakpoint

					//printf("Printing species list\n");
					sprintf(fn,"rep%03dsplist%03d.dat",rep,A.run_number);
					if((fp = fopen(fn,"w"))!=NULL){
						SP.print_spp_list(fp);
						fclose(fp);
					}
					else{
						printf("Unable to write to %s\n",fn);
					}
				}

				if(!(i%1000)){
					printf("%03d At  time %d e=%d,div=%d\n",rr,i,(int)A.energy,div);
					printsppct(&A,i);
				}


				thisepoch = A.get_ecosystem();

				if(!(i%1000)){
					if(thisepoch != lastepoch){
						lastepoch = thisepoch;
						nepochs++;
					}
				}


				if(!A.nagents(A.nowhead,-1) || (!rr && nsteps >1500000) || nsteps >15000000){
					printf("DEATH\n");
					printf("At  time %d e=%d,div=%d\t",i,(int)A.energy,div);
					//A.print_agents_count(stdout);
					A.print_spp_count(stdout,0,-1);
					nsteps=i;

					break;
				}
				nsteps++;
				A.energy += A.estep;
			}
			printf("Finished - alls well!\nclear out memory now:\n");
			fflush(stdout);

			A.clearout();
		}
	}
	//A.print_propensity(fsumm);
	fclose(fsumm);
	return 0;
}


/*
 * NB: To get identical trials to those run for ALifeXII, do the following:
 * 1: Use the file "replicase.conf" as the input
 * 2: Fix the random number seed to 437
 * 3: #define DO_ANCESTRY to get ancestry files out...
 * 4: Run on a 32-bit linux slackware system, circa 2008 vintage...:)
 */
int SmPm_AlifeXII(int argc, char *argv[]){

	int i,div;

	SMspp		SP;
	stringPM	A(&SP);
	FILE *fp;
	int indefinite=1;

	long rseed = initmyrand(-1);//437);//-1);//437);
	//R.printasint();

	unsigned int nsteps;// = (int) A.nsteps;

	//Prime the printout:
	//FILE *fpo,*fpdiv;
	FILE *fsumm,*ftmp;

	//division values for summary
	int divct,divit,diven;

	char	fn[128],pfn[128];

	sprintf(fn,"%s.spatial.summary.dat",argv[1]);
	if((fsumm=fopen(fn,"w"))==NULL){
		printf("Coundlnt open %s\n",fn);
		getchar();
	}
	fprintf(fsumm,"Random seed is %ld\n",rseed);
	fflush(fsumm);


	ftmp = fopen("epochs.dat","w");
	fclose(ftmp);

	unsigned int rlim = 1;
	readordef_param_int(argv[2], "NTRIALS", &rlim, rlim, 1);



	//for parallel runs, set a number as a prefix to each trial - so you can have 000,001,...,099 in one process and 100,101,...,199 in another
	int proc = 0;
	//if(argc > 2){
	//	proc = atoi(argv[3]);
	//}
	printf("process number is: %d\n",proc);

	//Read nsteps:
	unsigned int maxnsteps=0;

	int nns = readordef_param_int(argv[2], "NSTEPS", &maxnsteps, -1, 1);
	if(nns==1){
		printf("NSTEPS was not specified. Simulations will run indefinitely");
		indefinite = 1;
	}
	else{
		indefinite = 0;
	}



	for(unsigned int rr=0;rr<rlim;rr++){

		div=0;

		divct = divit = diven = 0;

		SP.clear_list();

		//A.load(argv[2],NULL,0,1);
		if(!arg_load(&A, argc, argv))
			return 0;

		unsigned int dummy;
		int gg = readordef_param_int(argv[2], "GRANULAR", &dummy, -1, 1);
		if(gg==1){
			printf("GRANULAR was not specified. Simulations will use standard operators");
			A.granular_1 = 0;
		}
		else{

			printf("GRANULAR was specified. Simulations will use+standard operators");
			A.granular_1 = 1;
		}




		A.print_agents(stdout,"NOW",0);

		A.run_number=rr;
		sprintf(pfn,"popdy%d%02d.dat",proc,A.run_number);
		ftmp = fopen(pfn,"w");
		fclose(ftmp);

#ifdef DO_ANCESTRY
		int lastepoch=A.get_ecosystem(),thisepoch,nepochs=1;
#endif

		A.domut=1;
		nsteps=0;
		for(i=0;indefinite || nsteps <= maxnsteps;i++){

			A.extit = i;

			A.make_next();
			A.update();

			//A.sanity_check();
			//if(!(i%10000)){
			//	check_pop_dy(pfn,&A,i);
			//}


			if(!(i%1000)){
				A.print_spp_count(stdout,0,-1);

			}

			if(!(i%1000)){

				printf("%d%02d At  time %d e=%d,div=%d, mutrate = %0.9f & %0.9f\n",proc,rr,i,(int)A.energy,div,A.subrate,A.indelrate);
				//A.print_agents_count(stdout);

				printsppct(&A,i);

			}


#ifdef DO_ANCESTRY
			thisepoch = A.get_ecosystem();

			if(!(i%1000)){
				if(thisepoch != lastepoch){
					lastepoch = thisepoch;
					nepochs++;

					FILE *afp;

					sprintf(fn,"ancestry%d%02d_%07d.txt",proc,rr,i);
					afp = fopen(fn,"w");
					//A.print_spp_strings(NULL);
					A.print_spp_strings(afp);
					printf("Done ancestry printing\n");//afp is closed in this function

					sprintf(fn,"ancestry%d%02d_%07d.dot",proc,rr,i);
					afp = fopen(fn,"w");
					A.print_ancestry_dot(afp,i,10000);

				}


				//printf("Printing species list\n");
				sprintf(fn,"splist%d%02d.dat",proc,A.run_number);
				if((fp = fopen(fn,"w"))!=NULL){
					SP.print_spp_list(fp);
					fclose(fp);
				}
				else{
					printf("Unable to write to %s\n",fn);
				}
			}

			//PRINT OUT THE EPOCH DATA:
			ftmp = fopen("epochs.dat","a");
			fprintf(ftmp,"%d\t%d\t%d\t%d\n",rr,i,A.spp_count-1,nepochs);
			fclose(ftmp);

#endif

			if(!A.nagents(A.nowhead,-1) || (!rr && nsteps >15000000) || nsteps >15000000){
				printf("DEATH\n");
				//fprintf(fpdiv,"%d\t%d\t%d",div,i,(int)A.energy);
				//A.print_agents_count(fpdiv);
				//A.print_spp_count(fpdiv,0,-1);
				//fflush(fpdiv);


				printf("At  time %d e=%d,div=%d, mutrate = %0.9f & %0.9f\t",i,(int)A.energy,div,A.indelrate,A.subrate);
				//A.print_agents_count(stdout);
				A.print_spp_count(stdout,0,-1);
				nsteps=i;

				break;
			}
			nsteps++;

			A.energy += A.estep;
		}



		if(divct)
			fprintf(fsumm,"%d\t%d\t%d\t%f\n",divct,divit,diven,(float) divct/divit);
		else
			fprintf(fsumm,"%d\t%d\t%d\t%f\n",div,i,(int)A.energy,-1.);
		fflush(fsumm);

		//Terminate the printout
		//fflush(fpo);
		//fclose(fpo);

		//fflush(fpdiv);
		//fclose(fpdiv);

		printf("Finished - alls well!\nclear out memory now:\n");


		//printf("Printing species list\n");
		sprintf(fn,"splist%d%02d.dat",proc,A.run_number);
		if((fp = fopen(fn,"w"))!=NULL){
			SP.print_spp_list(fp);
			fclose(fp);
		}
		else{
			printf("Unable to write to %s\n",fn);
		}


		fflush(stdout);

		A.clearout();
	}

	A.print_propensity(fsumm);

	////print the rule firings:
	//A.printfr(fsumm,&R);

	fclose(fsumm);
	return 0;

}



void printmaxcode(FILE *fp, int *ct, stringPM *A){
	for(int i=0;i<A->blosum->N;i++){
		fprintf(fp,"%02d, %c: %d\n",i,A->blosum->key[i],ct[i]);
	}
}


int comass_AlifeXII(int argc, char *argv[]){

	int i,div;

	SMspp		SP;
	stringPM	A(&SP);
	FILE *fp;
	int indefinite {1};

	long rseed = initmyrand(437);//-1);//437);
	//R.printasint();

	unsigned int nsteps = A.nsteps;

	//Prime the printout:
	//FILE *fpo,*fpdiv;
	FILE *fsumm,*ftmp;

	//division values for summary
	int divct,divit,diven;

	char	fn[128],pfn[128];

	sprintf(fn,"%s.spatial.summary.dat",argv[1]);
	if((fsumm=fopen(fn,"w"))==NULL){
		printf("Coundlnt open %s\n",fn);
		getchar();
	}
	fprintf(fsumm,"Random seed is %ld\n",rseed);
	fflush(fsumm);


	ftmp = fopen("epochs.dat","w");
	fclose(ftmp);

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
			printf("NTRIALS not sepcifid;\nSetting NTRIALS to %d\n",rlim);
			break;
		}
		fclose(fp);
	}
	printf("NTRIALS = %d\n",rlim);




	//Read nsteps:
	unsigned int maxnsteps=0;
	if((fp=fopen(argv[2],"r"))!=NULL){
		rerr = read_param_int(fp,"NSTEPS",&maxnsteps,1);
		switch(rerr){
		case 2:
			printf("Multiple NSTEPS specified. Check config file\n");
			getchar();
			exit(0);
		case 0:
			printf("Setting NSTEPS to %d\n",rlim);
			indefinite=0;
			break;
		default:
			printf("NSTEPS not sepcified;\nEach trial will run to extinction\n");
			break;
		}
		fclose(fp);
	}
	printf("NSTEPS = %d\n",maxnsteps);
	fflush(stdout);

	//counter for max no. symbols at any time step.
	int *maxcode;


	for(unsigned int rr=0;rr<rlim;rr++){

		FILE *mc;
		sprintf(fn,"maxcodes%03d.txt",rr);
		mc = fopen(fn,"w");
		fclose(mc);

		div= divct = divit = diven = 0;

		SP.clear_list();

		A.load(argv[2],NULL,0,1);
		A.load_comass(argv[2],1);
		A.print_agents(stdout,"NOW",0);

		A.run_number=rr;
		sprintf(pfn,"popdy%03d.dat",A.run_number);
		ftmp = fopen(pfn,"w");
		fclose(ftmp);

		if(!rr)
			maxcode = (int *) malloc(A.blosum->N * sizeof(int));
		memset(maxcode,0,A.blosum->N*sizeof(int));

		int lastepoch=A.get_ecosystem(),thisepoch,nepochs=1;

		A.domut=0;
		nsteps=0;
		for(i=0;indefinite || nsteps <= maxnsteps;i++){

			A.extit = i;

			A.comass_make_next();
			A.update();


			if(!(i%1000)){
				A.print_spp_count(stdout,0,-1);
				A.extit=i;//dummy line for breakpoint

				//printf("Printing species list\n");
				sprintf(fn,"splist%03d.dat",A.run_number);
				if((fp = fopen(fn,"w"))!=NULL){
					SP.print_spp_list(fp);
					fclose(fp);
				}
				else{
					printf("Unable to write to %s\n",fn);
				}
			}


			if(!(i%1000)){
				setmaxcode(&A,maxcode);
				printf("%03d At  time %d e=%d,div=%d\n",rr,i,(int)A.energy,div);
				printsppct(&A,i);
				for(int k=0;k<A.blosum->N;k++){
					printf("%c:\t%d\t%d",A.blosum->key[k],A.mass[k],maxcode[k]);
					if(A.mass[k]<0)
						printf(" ERROR\n");
					else
						printf("\n");

				}
			}


			thisepoch = A.get_ecosystem();


			if(!(i%1000)){
				if(thisepoch != lastepoch){
					lastepoch = thisepoch;
					nepochs++;
					/*
					FILE *afp;

					sprintf(fn,"ancestry%03d_%07d.txt",rr,i);
					afp = fopen(fn,"w");
					//A.print_spp_strings(NULL);
					A.print_spp_strings(afp);
					printf("Done ancestry printing\n");//afp is closed in this function

					sprintf(fn,"ancestry%03d_%07d.dot",rr,i);
					afp = fopen(fn,"w");
					A.print_ancestry_dot(afp,i,10000);
					*/
				}
			}


			if(!A.nagents(A.nowhead,-1)){// || (!rr && nsteps >1500000) || nsteps >15000000){
				printf("DEATH\n");
				printf("At  time %d e=%d,div=%d\t",i,(int)A.energy,div);
				A.print_spp_count(stdout,0,-1);
				nsteps=i;
				break;
			}
			nsteps++;

			A.energy += A.estep;
		}

		//PRINT OUT THE EPOCH DATA:
		/*
		ftmp = fopen("epochs.dat","a");
		fprintf(ftmp,"%d\t%d\t%d\t%d\n",rr,i,A.spp_count-1,nepochs);
		fclose(ftmp);
		if(divct)
			fprintf(fsumm,"%d\t%d\t%d\t%f\n",divct,divit,diven,(float) divct/divit);
		else
			fprintf(fsumm,"%d\t%d\t%d\t%f\n",div,i,(int)A.energy,-1.);
		fflush(fsumm);
		printf("Finished - alls well!\nclear out memory now:\n");
		fflush(stdout);
		*/

		printmaxcode(stdout,maxcode,&A);
		sprintf(fn,"maxcodes%03d.txt",rr);
		mc = fopen(fn,"a");
		printmaxcode(mc,maxcode,&A);
		fclose(mc);


		A.clearout();
	}

	A.print_propensity(fsumm);

	fclose(fsumm);
	return 0;

}



float EAqnn_summary(int rr){
	float score = 0;
	FILE *fp;
	char cp[256];

	sprintf(cp,"cp popdy%03d.dat popdy.dat",rr);
	system(cp);


	//printf("Calling R:\n");
	system("R -q -f file.R");

	//printf("Calling ls:\n");
	//system("ls -ltrh");

	//printf("Calling R:\n");
	//system("R --slave -e \"library(Rstringmol);x <- read.popdy(infn);xqnn<-qnn.activity(x);qnn<-qnn.summarize(xqnn);\"");

	fp = fopen("qnn.txt","r");

	fscanf(fp,"%e",&score);
	printf("Read qnn, value is %f\n",score);
	fclose(fp);

	return score;
}



void print_ga(FILE *fp, int n, int i, double *eval, int **params, int np, int lifetime){
	int j;
	fprintf(fp,"%d\t%d\t%e\t%d",i,n,eval[n],lifetime);
	for(j=0;j<np;j++){
		fprintf(fp,"\t%d",params[n][j]);
	}
	fprintf(fp,"\n");

	fflush(fp);
}


int * paramsFromFile(char *fn, const int rr, const int N){
	FILE *fp;
	fp = fopen(fn,"r");
	char line[2000];
	char *ptr;
	int *vals;
	vals = (int *) malloc(N*sizeof(int));
	//grab the line specified by rr
	for(int i=0;i<=rr;i++){
		fgets(line,2000,fp);
	}
	ptr=strtok(line,"\t");
	ptr=strtok(NULL,"\t");
	for(int i=0;i<N;i++){
		sscanf(ptr,"%d",&(vals[i]));
		ptr=strtok(NULL,"\t");
	}
	return vals;
}


double evalFromFile(char *fn, int rr){
	FILE *fp;
	fp = fopen(fn,"r");
	char line[2000];
	int ival;
	double val;
	//grab the line specified by rr
	for(int i=0;i<=rr;i++){
		fgets(line,2000,fp);
	}

	sscanf(line,"%d",&ival);
	return val=ival;

}


int comass_GA(int argc, char *argv[]){

	/* The idea here is to evolve appropriate masses to generate the most evolutionary activity. What we'll do is take
	 * a config file, but have a GA component which contains the total masses of each symbol. So this will be very like
	 * comass_ALifeXII above, but will be more iterative.
	 *
	 * To make this work, we'll need two things: A GA (Inman Harvey's MGA), and a fitness measure (Hickinbotham and Droops EA measure)
	 *
	 * The procedure is as follows:
	 * 1: Init genomes:
	 *   Create random genomes
	 *   Run comass
	 *   evaluate using EA
	 * 2: Then iterate through the GA to get an optimum mass value.
	 */

	const int POPN = 20;		//The population size
	const int NEVALS = 5000;  	//The total number of evaluations
	const int MMP = 12;
	const int MAXMASS = (int)pow(2,MMP); //=4096//TODO: read this from the config file perhaps..?

	double * eval;
	int **params;

	/*TODO: use
	void init_randseed_config(int argc, char *argv[]); to set up randseed
	 */

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
		//TODO: load the full RNG state using load_mt (RNGFILE in config)


		rerr = read_param_int(fpr,"GAQNN",&qnnscoring,1);
		if(rerr)
			qnnscoring=1;
		seedin = stmp;
	}
	fclose(fpr);

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

	//fprintf(stdout,"(unsigned) random seed is %lu \n",rseed);




	eval  = (double *)malloc(POPN * sizeof(double));
	params = (int **) malloc(POPN * sizeof(int *));

	SMspp		SP;
	stringPM	A(&SP);
	struct runparams R;

	clearfiles(argv);
	setupSMol(R, argc, argv);

	FILE *gafp;
	gafp = fopen("garun.txt","w");
	fclose(gafp);
	//1: create a bunch of random jobbies and evaluate them!

	int lifetime;
	int **mutnet;
	mutnet = (int **) malloc(POPN*sizeof(int *));

	/*INITIALISE: iterate over each member of the population */
	for(int rr=0;rr<POPN;rr++){

		FILE *mc;char fn[256];
		sprintf(fn,"maxcodes%03d.txt",rr);
		mc = fopen(fn,"w");
		fclose(mc);

		A.spl->clear_list();

		//todo: check we aren't violating comass when we init - what do we do if we do???
		if(!arg_load(&A, argc, argv))
			return 0;
		//A.load(argv[2],NULL,0,1);

		if(argc<4){

			//TODO: Much of this is repeated in the GA below - need to put it in a function and call it once!
			params[rr] = randintarray(A.blosum->N,0,MAXMASS);

			mutnet[rr] = randboolarray(A.blosum->N*A.blosum->N);
			setmutnet(mutnet[rr],A.blosum);

			A.set_mass(params[rr]);
			//A.print_agents(stdout,"NOW",0);

			lifetime = run_one_comass_trial(rr, &A, params[rr], &R);

			if(qnnscoring){
				printf("Using QNN to score fitness!\n");
				if(lifetime>45000)
					eval[rr] = -1 * EAqnn_summary(rr);
				else
					eval[rr]=0;
			}
			else{
				printf("Using Lifetime to score fitness!\n");
				eval[rr]=lifetime;
			}



		}
		else{
			params[rr] = paramsFromFile(argv[3],rr,A.blosum->N);
			eval[rr] = evalFromFile(argv[3],rr);
			lifetime = eval[rr];//todo this is a hack so we can extend run7, based on lifetime- we can put this in later if we want
		}


		//void print_ga(FILE *fp, int n, int i, double *eval, int **params, int np)
		gafp = fopen("garun.txt","a");
		print_ga(gafp,rr,0,eval,params,A.blosum->N,lifetime);
		fclose(gafp);

		A.clearout();
	}

	//printf("Seed Evaluations:\n");
	//for(int e=0;e<POPN;e++){
	//	printf("%02d\t%f\n",e,eval[e]);
	//}


	for(int gg=POPN; gg<NEVALS+POPN; gg++){
		//FILE *mc;char fn[256],pfn[256];FILE *ftmp;

		//sprintf(fn,"maxcodes%03d.txt",gg);
		//mc = fopen(fn,"w");
		//fclose(mc);

		A.spl->clear_list();

		//todo: check we aren't violating comass when we init - what do we do if we do???
		A.load(argv[2],NULL,0,1);

		//int ga_step_int(int **pop, double *eval, const int POPSIZE, const int PARAMETERS);
		//int ll;
		//int rr = ga_step_int(params, eval, POPN, A.blosum->N, 0, MAXMASS, &ll);
		int rr = ga_step_int(params, eval, POPN, A.blosum->N, 0, MAXMASS,NULL);

		A.set_mass(params[rr]);
		//A.print_agents(stdout,"NOW",0);

		lifetime = run_one_comass_trial(gg, &A, params[rr], &R);

		//eval[rr] = -1*  EAqnn_summary(gg);
		if(qnnscoring){
			if(lifetime>45000)
				//eval[rr] = -1 * EAqnn_summary(rr);
				//MERGE NOTE, 23/1/15: possibly the following is the correct configuration:
				eval[rr] = -1 * EAqnn_summary(gg);

			else
				eval[rr]=0;
		}
		else{
			eval[rr]=lifetime;
		}

		A.clearout();

		gafp = fopen("garun.txt","a");
		print_ga(gafp,rr,gg,eval,params,A.blosum->N,lifetime);
		fclose(gafp);

		if(!(gg%10)){
			printf("Fitness at %d Evaluations:\n", gg);
			for(int e=0;e<POPN;e++){
				printf("%02d\t%f\n",e,eval[e]);
			}
		}



	}
	fclose(gafp);
	return 0;
}

float *generate_avg_conc(char *fn, int *en){

	FILE *fp;


	float *avg;
	avg = NULL;

	printf("Opening %s..\n",fn);
	if((fp=fopen(fn,"r"))!=NULL){

		const int navg = 20; //TODO: read this in from somewhere;
		int line_number[navg];
		float fitness[navg];
		char line[2000];
		float minfit = 0;


		/*Initialise register */
		for(int i=0;i<navg;i++){
			line_number[i]=-1;
			fitness[i]=0;
		}

		/*Read in the data*/
		int lno = 0;
		while((fgets(line,2000,fp))!=NULL){
			/* Get the fitness value */
			float fit;
			lno++;

			sscanf(line,"%*d %*d %e",&fit);
			fit = -fit;

			//NB: fitness is negative...
			if(fit>minfit){
				for(int j=0;j<navg;j++){
					if(fit > fitness[j]){
						if(line_number[j]>-1){
							//We need to shunt everything down..
							for(int jj=navg-1;jj>j;jj--){
								fitness[jj]=fitness[jj-1];
								line_number[jj]=line_number[jj-1];
							}
						}
						line_number[j]=lno;
						fitness[j]=fit;
						break;
					}
				}
				minfit = fitness[navg-1];
			}
		}


		/*Calculate the averages*/
		const char *opcodes="ABC$DEF%GH^IJK?LMN}OPQ>RST=UVWXYZ";
		const int nops = strlen(opcodes);
		*en = nops;

		avg = (float *) malloc(nops*sizeof(float));
		memset(avg,0,nops*sizeof(float));


		for(int j=0;j<navg;j++){
			rewind(fp);
			int l=1;
			while(l!=line_number[j]){//TODO: Can't quite get this right without the extra fgets after the while - this'll fail - line_number[j]==1 (unlikely)
				fgets(line,2000,fp);
				l++;
			}
			fgets(line,2000,fp);
			char *lptr,*p;
			lptr = line;
			//Tokenize the line
			p = strtok(lptr, "\t");
			//skip through the data to get to the concs
			printf("Parsing run %s, on line %d, fitness %f\n",p,line_number[j],fitness[j]);
			p = strtok(NULL, "\t");
			p = strtok(NULL, "\t");
			p = strtok(NULL, "\t");
			//p = strtok(NULL, "\t");
			for(int c = 0;c<nops;c++){
				int tconc;
				p = strtok(NULL, "\t");
				sscanf(p,"%d",&tconc);
				avg[c] += (float)tconc;

			}
		}

		for(int c = 0;c<nops;c++){
			avg[c] /= navg;
			printf("Opcode %c has avg conc %f\n",opcodes[c],avg[c]);
		}

		fclose(fp);
	}
	else{
		printf("Unable to open garun file %s\n",fn);
		fflush(stdout);
	}

	return avg;
}



int * load_concs(int row, char *fn){

	/* We're going to go through a directory of files, using the path as a guide */
	//char fn[300];
	float *aconc;
	int *conc;
	int nops;
	//sprintf(fn,"%s\garun%03d.txt",path,row,&nops);

	//TODO: There are other ways to achieve this!
	aconc = generate_avg_conc(fn,&nops);

	conc = (int *) malloc(nops*sizeof(int));

	for(int n=0;n<nops;n++){
		conc[n] = floor(aconc[n]);
	}




	return conc;
}


int comass_GA_boostwinners(int argc, char *argv[]){

	/* Experiment 3 in comass paper, submitted to ecal 2015
	 * 1: take a high-fitness run with low conc of '='
	 * 2: boost conc of '=' by 1 order of magnitude
	 * 3: run the simulation
	 * 4: calc. qnn value - does it go up or down significantly...?
	 */

	FILE *gafp;

	init_randseed_config(argc,argv);

	SMspp		SP;
	stringPM	A(&SP);
	struct runparams R;

	clearfiles(argv);
	setupSMol(R, argc, argv);

	int crow=1;
	int *concs;

	char sfn[300];
	sprintf(sfn,"%s.summary.txt",argv[3]);

	FILE *sfp;
	sfp = fopen(sfn,"w");


	gafp = fopen("boostgarun.txt","w");
	fclose(gafp);

	//Go through each row in the concentrations file.
	//while( (concs = load_concs(crow,argv[3])) != NULL){
	concs = load_concs(crow,argv[3]);


	//TODO: increase the conc of '=':(check the right index is 26!)
	//concs[26] *= 10;
	for(int cc =0;cc<33;cc++){
		if(concs[cc]<400){
			concs[cc] *= 10;
			printf("Changed conc of opcode %d, to %d\n",cc,concs[cc]);
		}
	}


	for(int run=0;run<20;run++){

		for(int n=0;n<33;n++){
			if(n)fprintf(sfp,"\t");
			fprintf(sfp,"%d",concs[n]);
		}
		fprintf(sfp,"\n");


		float eval =0.;
		//FILE *mc;char fn[256],pfn[256];FILE *ftmp;

		//Reload the config:
		A.spl->clear_list();
		A.load(argv[2],NULL,0,1);


		A.set_mass(concs);

		int lifetime = run_one_comass_trial(crow, &A, concs, &R);

		if(lifetime>45000)
			eval = -1 * EAqnn_summary(crow);
			/*MERGE NOTE, 23/1/15: possibly the following is the correct configuration:
			eval[rr] = -1 * EAqnn_summary(gg);
			*/
		else
			eval=0;

		A.clearout();

		gafp = fopen("boostgarun.txt","a");
		//print_ga(gafp,crow,crow,eval,concs,A.blosum->N,lifetime);
		fprintf(gafp,"%d\t%d\t%e\t%d",crow,crow,eval,lifetime);
		for(int j=0;j<A.blosum->N;j++){
			fprintf(gafp,"\t%d",concs[j]);
		}
		fprintf(gafp,"\n");
		fflush(gafp);
		fclose(gafp);

		printf("Fitness for entry %d: %f, lifetime = %d\n", crow, eval, lifetime);

		crow++;
	}
	fclose(sfp);
	return 0;
}



int energetic_AlifeXII(int argc, char *argv[]){

	int i,div;

	SMspp		SP;
	stringPM	A(&SP);
	FILE *fp;
	int indefinite=1;

	long rseed = initmyrand(437);//-1);//437);
	//R.printasint();

	unsigned int nsteps = A.nsteps;

	//Prime the printout:
	//FILE *fpo,*fpdiv;
	FILE *fsumm,*ftmp;

	//division values for summary
	int divct,divit,diven;

	char	fn[128],pfn[128];

	sprintf(fn,"%s.spatial.summary.dat",argv[1]);
	if((fsumm=fopen(fn,"w"))==NULL){
		printf("Coundlnt open %s\n",fn);
		getchar();
	}
	fprintf(fsumm,"Random seed is %ld\n",rseed);
	fflush(fsumm);


	ftmp = fopen("epochs.dat","w");
	fclose(ftmp);

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
			printf("NTRIALS not sepcifid;\nSetting NTRIALS to %d\n",rlim);
			break;
		}
		fclose(fp);
	}
	printf("NTRIALS = %d\n",rlim);




	//Read nsteps:
	unsigned int maxnsteps=0;
	if((fp=fopen(argv[2],"r"))!=NULL){
		rerr = read_param_int(fp,"NSTEPS",&maxnsteps,1);
		switch(rerr){
		case 2:
			printf("Multiple NSTEPS specified. Check config file\n");
			getchar();
			exit(0);
		case 0:
			printf("Setting NSTEPS to %d\n",rlim);
			indefinite=0;
			break;
		default:
			printf("NSTEPS not sepcified;\nEach trial will run to extinction\n");
			break;
		}
		fclose(fp);
	}
	printf("NSTEPS = %d\n",maxnsteps);
	fflush(stdout);




	for(unsigned int rr=0;rr<rlim;rr++){

		div=0;

		divct = divit = diven = 0;

		SP.clear_list();

		A.load(argv[2],NULL,0,1);

		//test_adj(A.blosum);

		A.print_agents(stdout,"NOW",0);

		A.run_number=rr;
		sprintf(pfn,"popdy%03d.dat",A.run_number);
		ftmp = fopen(pfn,"w");
		fclose(ftmp);


		int lastepoch=A.get_ecosystem(),thisepoch,nepochs=1;

		A.domut=1;
		nsteps=0;
		for(i=0;indefinite || nsteps <= maxnsteps;i++){

			A.extit = i;

			A.energetic_make_next();
			A.update();

			//A.sanity_check();
			//if(!(i%10000)){
			//	check_pop_dy(pfn,&A,i);
			//}


			if(!(i%1000)){
				A.print_spp_count(stdout,0,-1);
				A.extit=i;//dummy line for breakpoint

				//printf("Printing species list\n");
				sprintf(fn,"splist%03d.dat",A.run_number);
				if((fp = fopen(fn,"w"))!=NULL){
					SP.print_spp_list(fp);
					fclose(fp);
				}
				else{
					printf("Unable to write to %s\n",fn);
				}
			}

			if(!(i%1000)){

				printf("%03d At  time %d e=%d,div=%d\n",rr,i,(int)A.energy,div);
				//A.print_agents_count(stdout);

				printsppct(&A,i);

			}


			thisepoch = A.get_ecosystem();

			if(!(i%1000)){
				if(thisepoch != lastepoch){
					lastepoch = thisepoch;
					nepochs++;

					/*
					FILE *afp;

					sprintf(fn,"ancestry%03d_%07d.txt",rr,i);
					afp = fopen(fn,"w");
					//A.print_spp_strings(NULL);
					A.print_spp_strings(afp);
					printf("Done ancestry printing\n");//afp is closed in this function

					sprintf(fn,"ancestry%03d_%07d.dot",rr,i);
					afp = fopen(fn,"w");
					A.print_ancestry_dot(afp,i,10000);
					*/
				}
			}


			if(!A.nagents(A.nowhead,-1) || (!rr && nsteps >1500000) || nsteps >15000000){
				printf("DEATH\n");
				printf("At  time %d e=%d,div=%d\t",i,(int)A.energy,div);
				//A.print_agents_count(stdout);
				A.print_spp_count(stdout,0,-1);
				nsteps=i;

				break;
			}
			nsteps++;

			//NB: used 50 energy here previously - make sure it is set in config!
			A.energy += A.estep;
		}

		//PRINT OUT THE EPOCH DATA:

		ftmp = fopen("epochs.dat","a");
		fprintf(ftmp,"%d\t%d\t%d\t%d\n",rr,i,A.spp_count-1,nepochs);
		fclose(ftmp);


		if(divct)
			fprintf(fsumm,"%d\t%d\t%d\t%f\n",divct,divit,diven,(float) divct/divit);
		else
			fprintf(fsumm,"%d\t%d\t%d\t%f\n",div,i,(int)A.energy,-1.);
		fflush(fsumm);

		printf("Finished - alls well!\nclear out memory now:\n");
		fflush(stdout);

		A.clearout();
	}

	A.print_propensity(fsumm);

	fclose(fsumm);
	return 0;

}




//Let's check what the pointers are doing:
void check_sagll(stringPM *A,int gclock, int c){
	FILE *fpx;
	char xfn[128];
	int y;
	sprintf(xfn,"list%03d_%02d.dat",gclock,c);
	fpx = fopen(xfn,"w");
	s_ag *pa;


	fprintf(fpx,"\nT %d CELL %d:\n",gclock,c);
	printf("\nCELL %d:\n",c);
	pa = A->nowhead;
	y=0;
	while(pa!=NULL){
		fprintf(fpx,"%d,\t%p,\t%d,\t%d,\t%p,\t%p,\t%p,\t%p\n",++y,pa,pa->idx,pa->status,pa->prev,pa->next,pa->exec,pa->pass);
		printf("%d,\t%p,\t%d,\t%d,\t%p,\t%p,\t%p,\t%p\n",y,pa,pa->idx,pa->status,pa->prev,pa->next,pa->exec,pa->pass);
		pa=pa->next;
	}

	fflush(fpx);
	fclose(fpx);
}




int SmPm_conpop(int argc, char *argv[]){


	//SMspp		SP;			//Global Species list:

	int 	NCON = 4;	//Number of containers

	const int 	MAXCONSTEPS = 10000000;	//Number of runsTODO: This should be many more!
	int c,c2,r=0;			//counters
	FILE *fp;
	float *score;

	score = (float *) malloc(NCON*sizeof(float));

	//437 divides at 167000
	FILE *fpr;
	unsigned long seedin = 0;
	if((fpr=fopen(argv[2],"r"))!=NULL){
		unsigned int stmp;
		int rerr = read_param_int(fpr,"RANDSEED",&stmp,1);
		//TODO: load the full RNG state using load_mt (RNGFILE in config)


		seedin = stmp;

		unsigned int nctmp;
		rerr = read_param_int(fpr,"NCONTAINERS",&nctmp,1);
		if(!rerr){
			NCON = nctmp;
		}

	}
	fclose(fpr);

	unsigned long rseed;


	if(seedin){
		rseed = longinitmyrand(&seedin);//437);//-1);//437);
	}
	else{
		rseed = longinitmyrand(NULL);
	}

	printf("rseed = %lu, seedin = %lu, NCON = %d\n",rseed,seedin,NCON);


	SMspp		*SP[NCON];
	stringPM	*A[NCON];	//Array of containers;

	//Prime the printout:
	FILE *fpdiv,*fsumm,*ftmp;

	char	fn[128],pfn[128];


	sprintf(fn,"conpopdy2.dat");
	fpdiv = fopen(fn,"w");
	fclose(fpdiv);



	sprintf(fn,"fitness.dat");
	fpdiv = fopen(fn,"w");
	fclose(fpdiv);

	sprintf(fn,"%s.spatial.summary.dat",argv[2]);
	if((fsumm=fopen(fn,"w"))==NULL){
		printf("Coundln't open %s\n",fn);
		getchar();
	}
	fprintf(fsumm,"Random seed is %ld\n",rseed);
	fflush(fsumm);


	ftmp = fopen("epochs.dat","w");
	fclose(ftmp);


	/** Set up the system */
	for(c=0;c<NCON;c++){
		SP[c] = new SMspp();
		A[c] = new stringPM(SP[c]);

		//A[c]->signal = signal;
		//A[c]->load(argv[2],NULL,0,1);
		if(!arg_load(A[c], argc, argv))
			return 0;


		A[c]->biomass=A[c]->bstart = 0;
		A[c]->domut=1;

		//sprintf(fn,"numbers%03d.dat",c);
		//fpo = fopen(fn,"w");

		//sprintf(fn,"divs%03d.dat",c);
		//fpdiv = fopen(fn,"w");
		//A[c]->print_agents_header(fpo);
		//A[c]->print_agents(stdout,"NOW",0);
		//fprintf(fpdiv,"div\te\t");
		//A[c]->print_agents_header(fpdiv);


		sprintf(pfn,"popdy%03d.dat",c);
		ftmp = fopen(pfn,"w");
		fclose(ftmp);

		//A[c]->set_epochs();
		A[c]->extit=0;
		A[c]->run_number=c;
		A[c]->energy=0;

	}

	r = NCON;

	//We're going to run NRUNSindividual containers, but in parallel
	int gclock=0;
	while(r<=MAXCONSTEPS){ //TODO: amend this so that immortal containers are possible

		for(c=0;c<NCON;c++){
			A[c]->make_next();
			A[c]->update();
			A[c]->extit++;
			//Restore the 20-each model - we'll do selection on replenishment
			A[c]->energy += A[c]->estep;
		}

		if(!(gclock%1000)){
			sprintf(fn,"conpopdy2.dat");
			fpdiv = fopen(fn,"a");

			printf("%d",gclock);
			fprintf(fpdiv,"%d",gclock);
			for(c=0;c<NCON;c++){
				printf("\t%d",c);
			}
			printf("\nCount:");
			for(c=0;c<NCON;c++){
				printf("\t%d",A[c]->nagents(A[c]->nowhead,-1));
				fprintf(fpdiv,"\t%d",A[c]->nagents(A[c]->nowhead,-1));
				printsppct(A[c],gclock);
				score[c] = ctspp(A[c],3);
			}
			//printf("\nScore:");
			//for(c=0;c<NCON;c++){
			//	printf("\t%g",score[c]);
			//}
			//printf("\n=======\n");

			fprintf(fpdiv,"\n");
			printf("\n");
			fclose(fpdiv);

			//Print the fitness
			sprintf(fn,"fitness.dat");
			fpdiv=fopen(fn,"a");
			fprintf(fpdiv,"%d",gclock);
			for(c=0;c<NCON;c++){
				fprintf(fpdiv,"\t%f",score[c]);
			}
			fprintf(fpdiv,"\n");
			fclose(fpdiv);

			for(c=0;c<NCON;c++){
				/*
				sprintf(fn,"conpop%02d_%07d.conf",c,gclock);
				fpdiv = fopen(fn,"w");
				A[c]->print_conf(fpdiv);
				fclose(fpdiv);
				*/

				/* We may need to reinstate this - but for now....
				printf("Generating Evostats for container %d\n",c);

				//Generate the stats we need for evolutionary distance from the origin...
				evostats(argv[2],A[c],&spp_matches,self,gvm);//,gvm,self);
				printf("Finished Generating Evostats for container %d\n",c);
				*/
			}

		}
		//Longer-term diagnostics:
		if(!(gclock%MAXCONSTEPS)){
			printf("Printing species lists\n");
			for(c=0;c<NCON;c++){
				sprintf(fn,"splistc%02dt%d.dat",c,gclock);
				fp = fopen(fn,"w");
				SP[c]->print_spp_list(fp);
				fclose(fp);//TODO: Possibly empty the list here??
			}
		}



		////Clear one out to test replenishment
		////Note: something in here stops the popn self-maintaining...
		//if(gclock==1000){
		//	A[8]->clearout();
		//	A[8]->load_table(argv[1]);
		//}

		//Let's check what the pointers are doing:
		//if(!(gclock%500)||gclock>=999){
			//check_sagll(A[1],gclock,1);
			//check_sagll(A[8],gclock,8);
		//}

		////DEBUG AGENT SHARING:
		//if(gclock==436){
		//	A[8]->share_agents(&(A[12]->nowhead));
		//}

		for(c=0;c<NCON;c++){
			if(!(A[c]->nagents(A[c]->nowhead,-1))){

				//Print the ancestry of this cell:
				/*
				FILE *afp;
				sprintf(fn,"ancestry%03d_%07d.dot",A[c]->r,gclock);
				afp = fopen(fn,"w");
				A[c]->print_ancestry_dot(afp,gclock,10000);
				*/



				////randomly pick another cell
				//c2=c;
				//while(c2==c)
				//	c2 = floor((float) NCON *rand0to1());

				//Pick the cell with the most agents in it.
				/*
				printf("Selecting fullest container...\n");
				int cmax = A[0]->nagents(A[0]->nowhead,-1);
				c2=0;
				for(int cc=1;cc<NCON;cc++){
					if(cc!=c){
						int tp = A[0]->nagents(A[cc]->nowhead,-1);
						if(tp>cmax){
							cmax=tp;
							c2=cc;
						}
					}
				}*/

				//Pick the cell with the highest score:
				printf("At time %d, Cell %d has died. Selecting fittest container for division...\n",gclock,c);
				/*
				set_energy(NCON, 20, A, signal,ea,score);
				for(int cc=0;cc<NCON;cc++){
					A[cc]->energy = 20;
				}
				*/


				/** TODO: this is the way we select the fittest cell - but we want random for now.
				c2=0;
				float cmax = score[c2];
				for(int cc=1;cc<NCON;cc++){
					if(cc!=c){
						//int tp = A[0]->nagents(A[cc]->nowhead,-1);
						if(score[cc]>cmax){
							cmax=score[cc];
							c2=cc;
						}
					}
				}*/

				c2 = c;
				while(c2==c){
					c2 = floor(NCON*rand0to1());
					if(c2==NCON){
						printf("Alert - rand0to1 returned 1!\n");
						c2 = NCON-1;
					}
					else{
						printf("Cell %d has been randomly selected to replace cell %d\n",c2,c);
					}
				}


				printf("Cell %d is fittest\n",c2);


				A[c]->share_agents(&(A[c2]->nowhead));
				A[c2]->energy /=2;
				A[c]->energy=A[c2]->energy;
				A[c]->run_number = r++;

				//clear the species list
				A[c]->spl->clear_list();
				//clear the list of seen reactions (because we can't track the spp nos
				free_swlist(&(A[c]->swlist));

				//TODO: refill the species list (set parents to NULL)
				s_ag *pag;
				l_spp *s;
				for(pag = A[c]->nowhead;pag!=NULL;pag=pag->next){
					pag->pp = A[c]->spl->make_parents(NULL,NULL);
					A[c]->update_lineage(pag,'I',1,NULL,NULL,0);
					s = A[c]->spl->getspp(pag,A[c]->extit,A[c]->maxl0);
					s->tspp = 0;
					pag->spp=s;
				}




				//////////////////////
				printf("Moved half of container %d (run %d) to container %d (run %d)\n",c2,A[c2]->run_number,c,A[c]->run_number);


				//sprintf(fn,"popdy%03d.dat",A[c]->r);
				fpdiv = fopen("splitstory.dat","a");
				fprintf(fpdiv,"%d: Run %d started in container %d from run %d in container %d\n",gclock,A[c]->run_number,c,A[c2]->run_number,c2);
				fclose(fpdiv);


				////print out the new pointer set
				//check_sagll(A[c],gclock+10000,c);
				//check_sagll(A[c2],gclock+10000,c2);
			}
		}

		gclock++;
	}

	//Print out the end of the ancestries
	for(c=0;c<NCON;c++){
		if(!(A[c]->nagents(A[c]->nowhead,-1))){

			//Print the ancestry of this cell:
			FILE *afp;
			sprintf(fn,"ancestry%03d_%07d.dot",A[c]->run_number,gclock);
			afp = fopen(fn,"w");
			A[c]->print_ancestry_dot(afp,gclock,10000);
		}
	}

	free(score);
	fclose(fsumm);
	return 0;

}



void check_setup(int argc, char *argv[]){
	//for now, all we'll do is check the parameters:

	test_config_settings(argc, argv, 0);
	return;

}


void SmPm_1on1(int argc,char *argv[]){

	SMspp		SP;
	stringPM A(&SP);
	s_ag *p,*a,*b;
	int i = 0;
	char c;
	char fn[256];
	FILE *fp;

	if(!arg_load(&A, argc, argv))
		return;
	/*switch(argc){
	case 3:
		printf("Traditional config\n");
		A.load(argv[2],NULL,1,1);
		break;
	case 4:
		printf("youShare-compatible config\n");
		A.load(argv[2],argv[3],1,1);
		break;
	default:
		printf("Error: wrong number of arguments - try 2 or 3\n");
		return;
	}*/

	A.verbose_bind=1;

	int maxn = 10000;
	int N = maxn<A.nsteps?maxn:A.nsteps;

	for(int n=0;n<N;n++){
		A.energy = 50;
		A.extit = i++;

		A.make_next();
		A.update();

		p = A.nowhead;
		while(p!=NULL){
		switch(p->status){
			case B_UNBOUND:
				break;
			case B_ACTIVE:
				a = p;
				b = p->pass;
				A.print_exec(stdout,a,b);
				break;
			case B_PASSIVE:
				//a = p->exec;
				//b = p;
				//A.print_exec(stdout,a,b);
				break;
			}
			p = p->next;
		}


		A.print_spp_count(stdout,0,-1);

		printf("Printing species list\n");
		sprintf(fn,"splist.dat");
		fp = fopen(fn,"w");
		SP.print_spp_list(fp);
		fclose(fp);

		if(argc==3){
			c=getchar();

			if(c=='q')
				return;
		}
	}

}


void swdist(int argc, char *argv[]){

	SMspp		SP;
	stringPM	A(&SP);
	FILE *fp;

	long rseed = initmyrand(437);

	FILE *fsumm,*ftmp;

	char	fn[128];

	sprintf(fn,"%s.spatial.summary.dat",argv[1]);
	if((fsumm=fopen(fn,"w"))==NULL){
		printf("Coundln't open %s\n",fn);
		getchar();
	}
	fprintf(fsumm,"Random seed is %ld\n",rseed);
	fflush(fsumm);


	ftmp = fopen("epochs.dat","w");
	fclose(ftmp);

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
			printf("NTRIALS not sepcifid;\nSetting NTRIALS to %d\n",rlim);
			break;
		}
		fclose(fp);
	}
	printf("NTRIALS = %d\n",rlim);



	//Read nsteps:
	unsigned int maxnsteps=0;
	if((fp=fopen(argv[2],"r"))!=NULL){
		rerr = read_param_int(fp,"NSTEPS",&maxnsteps,1);
		switch(rerr){
		case 2:
			printf("Multiple NSTEPS specified. Check config file\n");
			getchar();
			exit(0);
		case 0:
			printf("Setting NSTEPS to %d\n",rlim);
			break;
		default:
			printf("NSTEPS not sepcified;\nEach trial will run to extinction\n");
			break;
		}
		fclose(fp);
	}
	printf("NSTEPS = %d\n",maxnsteps);

	char ofn[512];
	char line[3000];
	FILE *fin, *outfile;

	A.load(argv[2],NULL,0,1);
	align sw;

	char s1[A.maxl0], s2[A.maxl0];

	if((fin= fopen(argv[3], "r"))!=NULL){
		sprintf(ofn,"%s.out",argv[3]);

		if((outfile= fopen(ofn, "w"))==NULL){
			printf("error opening %s",ofn);
		}

		fgets(line,A.maxl,fin);
		printf("%s\n",line);
		fprintf(outfile,"%s\n", line);

		while((fgets(line,A.maxl,fin))!=NULL){

			sscanf(line,"%*s\t%*s\t%*s\t%s\t%s",(char *) &s1, (char *) &s2);

			printf("read %s and %s\n",s1,s2);

			SmithWatermanV2(s1,s2,&sw,A.blosum,0);


			int l = strlen(line);
			line[l-1]='\0';


			printf("sw score is %f\n",sw.score);
			fprintf(outfile,"%s\t%f\t%d\t%d\n", line, sw.score,(int) strlen(s1), (int) strlen(s2));
			fflush(outfile);

			fprintf(stdout,"%s\t%f\t%d\t%d\n", line, sw.score,(int) strlen(s1),(int) strlen(s2));
			fflush(stdout);

		}


		fclose(fin);
		fclose(outfile);
	}
	else{
		printf("Couldn't open file %s",argv[2]);
	}



}



int speigpipette(stringPM *A, const int nmols, const int nrep, char *repstring, int massval){

	s_ag *pag;
	int count=0;
	int replen = strlen(repstring);

	//now replenish the replicases.
	for(int i=0;i<nrep;i++){

		l_spp *s;

		pag = A->make_ag('R');//,1);

		pag->S =(char *) malloc(A->maxl0*sizeof(char));
		pag->label = 'R';
		memset(pag->S,0,A->maxl0*sizeof(char));
		strncpy(pag->S,repstring,replen);
		pag->len = strlen(pag->S);

		//No parents for these initial agents!
		pag->pp = A->spl->make_parents(NULL,NULL);

		if(!i){
			A->update_lineage(pag,'R',1,NULL,NULL,0);
			s = A->spl->getspp(pag,0,A->maxl0);
			//TODO: tidy up handling of seed species, but for now:
			s->tspp = 0;
		}

		//Record that the replicase copies itself
		pag->spp=s;
		A->append_ag(&(A->nexthead),pag);
	}



	while(A->nowhead!=NULL && count < nmols){
		pag = A->rand_ag(A->nowhead,-1);
		A->extract_ag(&(A->nowhead),pag);
		A->update_lineage(pag,'R',1,NULL,NULL,0);

		//safest & quickest to destroy the replicases and replenish.
		if(!(strncmp(pag->spp->S,repstring,replen))){
			A->free_ag(pag);
			continue;
		}
		else{
			if(pag->status!=B_UNBOUND){
				A->free_ag(pag);
				continue;
			}
		}
		A->update_lineage(pag,'M',1,NULL,NULL,0);
		A->append_ag(&(A->nexthead),pag);
		count++;
	}



	A->update();
	for(int i=0;i<A->blosum->N;i++){
		A->mass[i]=massval;
	}
	for(pag=A->nowhead;pag!=NULL;pag=pag->next){
		A->update_mass(pag->S,strlen(pag->S),-1,1);
	}
	return 1;
}



int speigmonst(int argc, char *argv[]){

	int i,div;

	SMspp		SP;
	stringPM	A(&SP);
	FILE *fp;
	int indefinite=1;

	long rseed = initmyrand(437);//-1);//437);
	//R.printasint();

	unsigned int nsteps = A.nsteps;

	//Prime the printout:
	//FILE *fpo,*fpdiv;
	FILE *fsumm,*ftmp;

	//division values for summary
	int divct,divit,diven;

	char	fn[128],pfn[128];

	sprintf(fn,"%s.spatial.summary.dat",argv[1]);
	if((fsumm=fopen(fn,"w"))==NULL){
		printf("Coundlnt open %s\n",fn);
		getchar();
	}
	fprintf(fsumm,"Random seed is %ld\n",rseed);
	fflush(fsumm);


	ftmp = fopen("epochs.dat","w");
	fclose(ftmp);

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
			printf("NTRIALS not sepcifid;\nSetting NTRIALS to %d\n",rlim);
			break;
		}
		fclose(fp);
	}
	printf("NTRIALS = %d\n",rlim);




	//Read nsteps:
	unsigned int maxnsteps=0;
	if((fp=fopen(argv[2],"r"))!=NULL){
		rerr = read_param_int(fp,"NSTEPS",&maxnsteps,1);
		switch(rerr){
		case 2:
			printf("Multiple NSTEPS specified. Check config file\n");
			getchar();
			exit(0);
		case 0:
			printf("Setting NSTEPS to %d\n",rlim);
			indefinite=0;
			break;
		default:
			printf("NSTEPS not sepcified;\nEach trial will run to extinction\n");
			break;
		}
		fclose(fp);
	}
	printf("NSTEPS = %d\n",maxnsteps);
	fflush(stdout);

	//counter for max no. symbols at any time step.
	int *maxcode;


	char *repstring;
	repstring = (char*)malloc(A.maxl0 * sizeof(char));
	for(unsigned int rr=0;rr<rlim;rr++){

		FILE *mc;
		sprintf(fn,"maxcodes%03d.txt",rr);
		mc = fopen(fn,"w");
		fclose(mc);

		div= divct = divit = diven = 0;

		if(!rr){
			A.load(argv[2],NULL,0,1);
			A.load_comass(argv[2],1);
			strncpy(repstring,SP.species->next->S,A.maxl0);
		}
		else{
			SP.clear_list();
			speigpipette(&A, 50, 50, repstring, 3000);//TODO: fix the massval thing (once comass_ga experiments are finished)
		}

		A.energy=20;
		A.print_agents(stdout,"NOW",0);

		A.run_number=rr;
		sprintf(pfn,"popdy%03d.dat",A.run_number);
		ftmp = fopen(pfn,"w");
		fclose(ftmp);

		if(!rr)
			maxcode = (int *) malloc(A.blosum->N * sizeof(int));
		memset(maxcode,0,A.blosum->N*sizeof(int));

		int lastepoch=A.get_ecosystem(),thisepoch,nepochs=1;

		A.domut=1;

		//turn decay off!
		A.dodecay=0;


		nsteps=0;
		for(i=0;indefinite || nsteps <= maxnsteps;i++){

			A.extit = i;

			A.comass_make_next();
			A.update();


			if(!(i%1000)){
				A.print_spp_count(stdout,0,-1);
				A.extit=i;//dummy line for breakpoint

				//printf("Printing species list\n");
				sprintf(fn,"splist%03d.dat",A.run_number);
				if((fp = fopen(fn,"w"))!=NULL){
					SP.print_spp_list(fp);
					fclose(fp);
				}
				else{
					printf("Unable to write to %s\n",fn);
				}
			}


			if(!(i%1000)){
				setmaxcode(&A,maxcode);
				printf("%03d At  time %d e=%d,div=%d\n",rr,i,(int)A.energy,div);
				printsppct(&A,i);
				for(int k=0;k<A.blosum->N;k++){
					printf("%c:\t%d\t%d",A.blosum->key[k],A.mass[k],maxcode[k]);
					if(A.mass[k]<0)
						printf(" ERROR\n");
					else
						printf("\n");

				}
			}


			thisepoch = A.get_ecosystem();


			if(!(i%1000)){
				if(thisepoch != lastepoch){
					lastepoch = thisepoch;
					nepochs++;
					/*
					FILE *afp;

					sprintf(fn,"ancestry%03d_%07d.txt",rr,i);
					afp = fopen(fn,"w");
					//A.print_spp_strings(NULL);
					A.print_spp_strings(afp);
					printf("Done ancestry printing\n");//afp is closed in this function

					sprintf(fn,"ancestry%03d_%07d.dot",rr,i);
					afp = fopen(fn,"w");
					A.print_ancestry_dot(afp,i,10000);
					*/
				}
			}


			if(!A.nagents(A.nowhead,-1)){// || (!rr && nsteps >1500000) || nsteps >15000000){
				printf("DEATH\n");
				printf("At  time %d e=%d,div=%d\t",i,(int)A.energy,div);
				A.print_spp_count(stdout,0,B_UNBOUND);
				nsteps=i;
				break;
			}
			nsteps++;

			A.energy += A.estep;
		}

		printf("\n+++++++++++++++++\nUNBOUND SPECIES:");
		A.print_spp_count(stdout,0,B_UNBOUND);
		printf("\n+++++++++++++++++\nUNBOUND SPECIES:");

		//PRINT OUT THE EPOCH DATA:
		/*
		ftmp = fopen("epochs.dat","a");
		fprintf(ftmp,"%d\t%d\t%d\t%d\n",rr,i,A.spp_count-1,nepochs);
		fclose(ftmp);
		if(divct)
			fprintf(fsumm,"%d\t%d\t%d\t%f\n",divct,divit,diven,(float) divct/divit);
		else
			fprintf(fsumm,"%d\t%d\t%d\t%f\n",div,i,(int)A.energy,-1.);
		fflush(fsumm);
		printf("Finished - alls well!\nclear out memory now:\n");
		fflush(stdout);
		*/

		printmaxcode(stdout,maxcode,&A);
		sprintf(fn,"maxcodes%03d.txt",rr);
		mc = fopen(fn,"a");
		printmaxcode(mc,maxcode,&A);
		fclose(mc);


		//A.clearout();
	}

	A.print_propensity(fsumm);

	fclose(fsumm);
	return 0;
}


int test_config(){



	return 0;

}






int main(int argc, char *argv[]) {

	if(argc>1){


		printf("Hello There! Let's run stringmol!\n");

		printf("Argc = %d\n",argc);

		int trial = atoi(argv[1]);
		switch(trial){
		case 0: //Check's a molecules interaction
			SmPm_1on1(argc,argv);
			break;

		case 1:
			//printf("First thing to do is run two pops on the alifeXII model\n(Running two to check mem alloc)\n");
			SmPm_AlifeXII(argc,argv);
			break;

		case 2:
			//Now let's try a container population
			SmPm_conpop(argc,argv);
			break;

		case 3:
			//Now let's try a container population
			origlife(argc,argv);
			break;

		case 4:
			comass_GA(argc,argv);
			break;

		case 5:
			joinsplists(argc,argv);
			break;

		case 6:
			//printf("First thing to do is run two pops on the alifeXII model\n(Running two to check mem alloc)\n");
			energetic_AlifeXII(argc,argv);
			break;

		case 7:
			swdist(argc,argv);
			break;

		/*************************************************/
		case 8:
			comass_AlifeXII(argc,argv);
			//comass_GA(argc,argv);
			break;

		/*************************************************/
		case 9://UNFINISHED: stringmol emulation of Speigelmen's Monster
			speigmonst(argc,argv);
			break;

		/*************************************************/
		case 10://TODO: This is probably obsolete now we have test_all() (case 99)
			check_setup(argc,argv);
			break;

		/*************************************************/
		case 33://Spatial stringmol experiments, summer 2016
			smspatial(argc,argv);
			break;

		/*************************************************/
		case 34://Analyse spatial stringmol experiments - generate ancestry
			smspatial_ancestry(argc,argv);
			break;

		/*************************************************/
		case 35://Analyse spatial stringmol experiments - generate length images
			smspatial_lengthpicsfromlogs(argc,argv);
			break;

		/*************************************************/
		case 36://Analyse spatial stringmol experiments - generate community from config file
			smspatial_community(argc,argv);
			break;


		/*************************************************/
		case 44://For final experiment in the comass experiment, spring 2015
			comass_GA_boostwinners(argc,argv);
			break;

		/*************************************************/
		case 99:
			test_all(argc,argv);
			break;


		/*************************************************/
		case 901:

			test_rand_config(argc,argv);
			break;

		}
		printf("Finished!\n");
		fflush(stdout);
	}
	else{
		printf("\nUSAGE for TestSM\n\nThe first argument after TestSM should be the trial type, followed by the arguments needed to run it\n\n");
		printf("TRIAL TYPES LIST (NUMBERS IN BRACKETS ARE EXPERIMENTAL!)\n");
		printf("NAME             TTYPE        ARGUMENTS\n\n");
		printf("1 on 1            0           1: .conf;  (2: .mtx)\n");
		printf("ALife XII         1           1: .conf;  (2: .mtx)\n");
		printf("Con Pop           2           1: .conf;  (2: .mtx)\n");
		printf("Origlife         (3)\n");
		printf("Comass GA         4           1: .conf;\n");
		printf("Joinsplists      (5)\n");
		printf("Energetic ALXII  (6)\n");
		printf("swdist           (7)\n");
		printf("Comass ALXII     (8)\n");
		printf("speigmonst       (9)\n");
		printf("Check setup       10          1: .conf;  (2: .mtx)\n");
		printf("Spatial Stringmol(33)         1: .conf;  (2: .mtx)\n");
		printf("Comass GA boost   44          1: .conf;   2: boost.dat\n\n\n");
		printf("Test All          99          1: .conf;  (2: .mtx)\n\n\n");
		printf("Test RNG         (901)       (1: .conf;) (2: .mtx)\n\n\n");
	}

	return 0;
}
