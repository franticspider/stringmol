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

//#include "trigutil.h"
extern "C" {

	#include "memoryutil.h" //for memerror
	#include "randutil.h"

	#include "params.h"
}

#include "rules.h"
#include "agents_base.h"



int agents_base::eqn_prop(const int n){

	float agarea = (float) M_PI*pow((float) agrad,2);
	float cellarea = M_PI*pow((float) vcellrad-(2.*(float) move),2);
	float arearatio = agarea/cellarea;


	if(n){
		//reactant coverage = 1 - ( 1- (area of reactant/area of cell) )^(number of reactants)
		float cov = 1.-pow(1.-arearatio,n) ;

		//printf("%d\t%f\t%f\t%f\t%f\n",n,cov,agarea,cellarea,arearatio);

		float rno = rand0to1();

		if(rno<cov)
			return 1;
		else
			return 0;
	}
	else{
		//printf("Zero mean coverage\n");
		return 0;
	}
}




//creators and destructors
agents_base::agents_base(){

	//Need to nullify these because stringPM doesn't use them
	//aat=NULL;
	//aac=NULL;
	//adc=NULL;
	//aro=NULL;

	//for preset
	bmax = 100;
	bct = (int *) malloc(bmax * sizeof(int));
	bpp = (int *) malloc(bmax * sizeof(int));
	memset(bct,0,bmax*sizeof(int));
	memset(bpp,0,bmax*sizeof(int));
	fr= NULL;
	tr = NULL;

	/*
	int i,j,k,found;
	float x1,y1,x2,y2;
	float dist,rad = 400-(9.9*2);
	for(i=0;i<bmax;i++){
		printf("Doing %d\n",i);
		for(j=0;j<100000;j++){
			if(!j%10000){
				printf("*");fflush(stdout);
			}
			found=0;
			rand_in_rad(rad,&x1,&y1);
			for(k=0;k<i;k++){
				rand_in_rad(rad,&x2,&y2);
				dist = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
				if(dist<9.9){
					found = 1;
				}
			}
			note_propensity(i,found);
		}
		printf("\n");
	}

	print_propensity(stdout);
	exit(7);
	*/

	preset();

}

agents_base::~agents_base(){

	clearout(0);

	free(bct);
	free(bpp);
}


void agents_base::note_propensity(int N,int X){
	if(N<bmax){
		bct[N]++;
		if(X)
			bpp[N]++;
	}

}


void agents_base::print_propensity(FILE *fp){
	int i;
	fprintf(fp,"\nPropensity table\n");
	for(i=0;i<bmax;i++)
		fprintf(fp,"%d\t%d\t%d\n",i,bct[i],bpp[i]);
	fflush(fp);
}


void agents_base::preset(){

	ifxhead = NULL;
	//btab = NULL;
	//com = NULL;
	//dcom = NULL;


	cellrad = 2500;
	agrad = 10;
	move = 0;
	energy = 0;
	//divtime = 0;
	vcellrad = 0;
}

//Load from a file - this should be common - but load agents can't be overwritten!!
//int agents_base::load(char *fn){

	//NB: Error checking required!!

	//load_params(fn);
	//load_influx(fn);

	//Load Agents now!

	//load_division(fn);
	//load_replenish(fn);
	//return 1;
//}



// VJH - added this function for youShare - some alterations on merge
void agents_base::load(const char *fn, char *fninput, int test=0, int verbose=0){

	load_params(fn,test,verbose);
	load_influx(fn);

	// VJH load_agents(fn,test);
	load_agents(fn,fninput,test);
}





int agents_base::load_params(const char *fn, int test, int verbose){

	FILE *fp;
	int e,err = 0;

	if((fp=fopen(fn,"r"))!=NULL){
		float tmpen;
		e=  read_param_float(fp,"CELLRAD",&cellrad, verbose);
		if(e>1)err++;

		vcellrad = cellrad;
		e=  read_param_float(fp,"AGRAD",&agrad, verbose);
		if(e>1)err++;

		//err +=  read_param_float(fp,"MOVE",&move);
		e=  read_param_float(fp,"ENERGY",&tmpen, verbose);
		if(e>1)err++;
		energy = (int) tmpen;

		e=  read_param_float(fp,"NSTEPS",&nsteps, verbose);
		if(e>1)err++;
		//err += 	read_param_float(fp,"DIVTIME",&divtime);

		if(err){
			printf("Some error reading config file\n");
			fclose(fp);
			return 0;
		}
		else{
			if(test){//Make sure things will collide
				vcellrad=cellrad=agrad;
			}
			fclose(fp);
			return 1;
		}
	}
	else{
		printf("Unable to open file %s\n",fn);
		fflush(stdout);
		return 0;
	}

}


int agents_base::load_influx(const char *fn){
	const int maxl = 128;
	FILE *fp;
	char line[maxl];
	char label[maxl];
	s_ix *pix;
	int n,start,stop;
	float prob;
	char lab;

	//WHY???
	////reset the tracking count
	//ntt=0;

	if((fp=fopen(fn,"r"))!=NULL){
		while((fgets(line,maxl,fp))!=NULL){
			memset(label,0,maxl);
			sscanf(line,"%s",label);
			//printf("line = %s",line);
			if(!strncmp(line,"INFLUX",6)){
				//sscanf(line,"%*s %c %d %d %d",&lab,&n,&start,&stop);
				//printf("INFLUX: %c %d %d %d\n",lab,n,start,stop);
				//pix = make_influx(lab,n,start,stop);

				sscanf(line,"%*s %c %d %f %d %d",&lab,&n,&prob,&start,&stop);
				printf("INFLUX: %c %d %f %d %d\n",lab,n,prob,start,stop);
				pix = make_influx(lab,n,prob,start,stop);

				if(ifxhead == NULL){
					ifxhead = pix;
				}
				else{
					append_ix(&ifxhead,pix);
					}

			}
		}
		fclose(fp);
		return 1;
	}
	else{
		printf("Unable to open file %s\n",fn);
		fflush(stdout);
		return 0;
	}

}



/*
int agents_base::load_division(char *fn){
	const int maxl = 128;
	FILE *fp;
	char line[maxl];
	char label[maxl];
	int i,n;
	char lab;

	if((fp=fopen(fn,"r"))!=NULL){
		while((fgets(line,maxl,fp))!=NULL){
			memset(label,0,maxl);
			sscanf(line,"%s",label);
			//printf("line = %s",line);
			if(!strncmp(line,"DIVIDE",6)){

				sscanf(line,"%*s %c %d",&lab,&n);
				printf("DIVIDE if %c >= %d\n",lab,n);

				if(aat!=NULL){
					for(i=0;i<ntt;i++){
						if(lab==aat[i])
							adc[i]=n;
					}
				}
				else{
					printf("WARNING: aat not initialised for setting division conditions\n");
					fflush(stdout);
				}
			}


			//if(!strncmp(line,"DCOMPLEX",8)){
			//	sscanf(line,"%*s %c %c",&lab,&lab2);
			//	dcom[nttindex(lab)][nttindex(lab2)]=1;
			//}

		}
		fclose(fp);
		return 1;
	}
	else{
		printf("Unable to open file %s\n",fn);
		fflush(stdout);
		return 0;
	}

}
*/



/*
int agents_base::load_replenish(char *fn){
	const int maxl = 128;
	FILE *fp;
	char line[maxl];
	char label[maxl];
	int i,n;
	char lab,lab2;

	if((fp=fopen(fn,"r"))!=NULL){
		while((fgets(line,maxl,fp))!=NULL){
			memset(label,0,maxl);
			sscanf(line,"%s",label);
			//printf("line = %s",line);
			if(!strncmp(line,"REPLENISH",9)){

				sscanf(line,"%*s %c %d",&lab,&n);
				printf("REPLENISH %c %d times\n",lab,n);

				for(i=0;i<ntt;i++){
					if(lab==aat[i])
						aro[i]=n;
				}
			}
			if(!strncmp(line,"RCOMPLEX",8)){
				sscanf(line,"%*s %c %c",&lab,&lab2);
				com[nttindex(lab)][nttindex(lab2)]=1;
			}

		}
		fclose(fp);

		return 1;
	}
	else{
		printf("Unable to open file %s\n",fn);
		fflush(stdout);
		return 0;
	}
}
*/


/*
int agents_base::nttindex(const int label){

	// '*' means 'empty' or 'no agent' in this code - see rules.cpp!
	if(label=='*')
		return -1;
	for(int i=0;i<ntt;i++){
		if(label==aat[i])
			return i;
	}
	return -1;

}
*/

/*
void agents_base::make_com(){
	int i;
	if(com==NULL){
		com = (int **)  mymalloc(ntt,sizeof(int *));
		for(i=0;i<ntt;i++){
			com[i] = (int *) mymalloc(ntt,sizeof(int));
			memset(com[i],0,ntt*sizeof(int));
			com[i][i]=1;
		}
	}
}
*/

/*
void agents_base::make_dcom(){
	int i;
	if(dcom==NULL){
		dcom = (int **)  mymalloc(ntt,sizeof(int *));
		for(i=0;i<ntt;i++){
			dcom[i] = (int *) mymalloc(ntt,sizeof(int));
			memset(dcom[i],0,ntt*sizeof(int));
			dcom[i][i]=1;
		}
	}
}
*/


/*
void agents_base::make_btab(rules *rset){

	int i,r,A,B;
	i=0;
	if(btab==NULL){
		btab = (int **)  mymalloc(ntt,sizeof(int *));
		for(i=0;i<ntt;i++){
			btab[i] = (int *) mymalloc(ntt,sizeof(int));
			memset(btab[i],0,ntt*sizeof(int));
		}
	}

	//now go through the rule table...
	for(r=0;r<rset->nr;r++){
		if((A = nttindex(rset->rset[r][0]))>-1){
			if((B = nttindex(rset->rset[r][1]))>-1){
				btab[A][B]=1;
				btab[B][A]=1;
			}
		}
	}
}
*/



int agents_base::append_ix(s_ix **list, s_ix *ax){
	s_ix *pix;

	//printf("appending: list = %p, ag = %p\n",*list,ag);
	if(*list==NULL){
		*list=ax;
		//printf("appended: list = %p, ag = %p\n",*list,ag);
	}
	else{
		pix = *list;
		while(pix->next != NULL){
			pix = pix->next;
		}
		pix->next = ax;
		//ag->prev = pag;
	}
	return 0;
}







s_ix * agents_base::make_influx(int lab, int n, float prob, int start, int stop){
	s_ix *ix;

	if((ix = (s_ix *) mymalloc(1,sizeof(s_ix)))!=NULL){
		ix->n = n;
		ix->prob = prob;
		ix->start = start;
		ix->stop = stop;
		ix->label = lab;
		ix->next = NULL;
		return ix;
	}
	else{
		return NULL;
	}

}




void agents_base::clearout(int verbose){

	s_ix	*ixp,*ixp2;

	if(verbose){printf("Starting agents_base clearout..");fflush(stdout);}

	ixp=ifxhead;
	while(ixp!=NULL){
		ixp2=ixp->next;
		free(ixp);
		ixp=ixp2;
	}

	preset();

	if(verbose){printf("....finished\n");fflush(stdout);}
}



