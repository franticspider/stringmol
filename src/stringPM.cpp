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
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "memoryutil.h"
#include "randutil.h"
#include "params.h"
#include "hsort.h"

//string stuff
#include "stringmanip.h"
#include "alignment.h"
#include "instructions.h"

//metabolism stuff
#include "rules.h"
#include "agents_base.h"
#include "SMspp.h"
#include "stringPM.h"

//For debugging:
//#define VERBOSE
//#define V_VERBOSE

//Sticky stringmol:
//#define BIND_ALL

//For decay experiments in stringmol ALife paper:
//#define LONG_DECAY
//#define UNB_DECAY_ONLY

//Brutal hack for line length
//extern const int maxl = 2000;
//extern const int maxl0 = maxl+1; //allow room for a terminating 0


void print_status(FILE *fp,s_bind st){
	switch(st){
	case B_UNBOUND:
		fprintf(fp,"UNBOUND");
		break;
	case B_ACTIVE:
		fprintf(fp,"ACTIVE");
		break;
	case B_PASSIVE:
		fprintf(fp,"PASSIVE");
		break;
	default:
		fprintf(fp,"UNDEFINED STATUS");
	}
	fflush(fp);
}


stringPM::stringPM(SMspp * pSP){

	if(pSP!=NULL)
		spl = pSP;
	else{
		spl = new SMspp();
	}

	dodecay=1;

	swlist=NULL;

	blosum = NULL;
	blosum = (swt *) malloc(sizeof(swt));
	blosum->N=0;
	blosum->T=NULL;
	blosum->key=NULL;
	preset();
	agct = 0;
	//species=NULL;
	spp_count=1;
	verbose_bind=0;

	//Set defaults:

	maxl = 2000;
	maxl0 = maxl+1; //allow room for a terminating 0
	estep = 20;

	signal = NULL;

	/** This is a toggle to turn the '+' operator on and off
	 *
	 * It can be set using the 'GRANULAR' flag in the config file
	 *
	 */
	granular_1=0;

	splprint = 10000;

}


stringPM::~stringPM() {
	clearout(0);
}


void stringPM::preset(){

	nowhead = NULL;
	nexthead = NULL;

	//note - it's possible that this'll be called twice - but will do no harm!
	agents_base::preset();
}

/*
int stringPM::load(char *fn){

	load_params(fn);
	load_influx(fn);

	load_agents(fn);

	make_com();
	make_dcom();

	load_division(fn);
	load_replenish(fn);
	return 1;
}
*/


void stringPM::testprop(){

	int i,j,k,found;
	float x1,y1,x2,y2;
	float dist,rad = 400-(9.9*2);


	memset(bct,0,bmax*sizeof(int));
	memset(bpp,0,bmax*sizeof(int));

	for(i=0;i<bmax;i++){
		printf("Doing %d\n",i);
		for(j=0;j<100000;j++){
			if(!(j%10000)){
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


}

/* Make sure the masseage memory is freed - so you can't do printf("%s",parse_error(1)); - the memory will not be freed.
 *
 */
char * stringPM::parse_error(int errno){

	char * message;
	message = (char *) malloc(100*sizeof(char));

	memset(message,0,100*sizeof(char));

	sprintf(message,"Unspecified error");

	return message;

}


float stringPM::load_mut(char *fn, int verbose){

	FILE *fp;
	float mut;
	char *emsg;
	int finderr=1;

	if((fp=fopen(fn,"r"))!=NULL){

		finderr=read_param_float(fp,"MUTATE",&mut,verbose);
		fclose(fp);

		switch(finderr){
		case 0:
			subrate = mut;
			indelrate = mut;
			if(mut<FLT_MIN)
				domut=0;
			else
				domut=1;
			if(verbose)printf("MUTATE setting: subrate = %f, indelrate = %f\n",subrate,indelrate);
			return 0;
		case 1:
			printf("No mutation rate found in config file\nUsing ALifeXII values instead\n");
			indelrate = 0.0000000306125;
			subrate=0.00001;
			domut=1;
			return 0;
			break;
		default: //Some other error
			printf("Error %d(%s) in loading mutation scheme\n",finderr,emsg=parse_error(finderr));
			free(emsg);
			return 1;
			break;

		}
	}
	return 2;
}




float stringPM::load_decay(char *fn, int verbose){

	FILE *fp;
	float dec;
	int finderr=1;

	if((fp=fopen(fn,"r"))!=NULL){

		finderr=read_param_float(fp,"DECAY",&dec,verbose);
		fclose(fp);

		switch(finderr){
		case 0:
			decayrate = dec;
			if(verbose)printf("DECAY rate setting = %f\n",decayrate);
			return 0;
		case 1:
			printf("No mutation rate found in config file\nUsing ALifeXII values instead\n");
			decayrate = 1.0/pow(65,2);
			return 0;
			break;
		default: //Some other error
			printf("Error %d in loading mutation scheme\n",finderr);
			return 1;
			break;

		}
	}
	return 2;
}




//TODO: This should be in alignment.cpp
int stringPM::load_table_matrix(char *fn){
	FILE *fp;
	char line[maxl];
	char label[maxl];
	char *p;
	int i,j;
	printf("File name is %s",fn);
	if((fp=fopen(fn,"r"))!=NULL){
		fgets(line,maxl,fp);
		memset(label,0,maxl);
		sscanf(line,"%s",label);

		//Set N and key:
		blosum->N = strlen(label);
		blosum->key = (char *)malloc((blosum->N+1) * sizeof(char));
		memset(blosum->key,0,(blosum->N+1)*sizeof(char));
		strncpy(blosum->key,label,blosum->N);

		//Set up the data structure
		blosum->T = (float **) malloc((blosum->N +1) * sizeof(float *));
		for(i=0;i<blosum->N+1;i++){
			blosum->T[i] = (float *) malloc(blosum->N * sizeof(float));
		}

		for(i=0;i<blosum->N+1;i++){
			fgets(line,maxl,fp);
			p=strtok(line,", \t");
			for(j=0;j<blosum->N;j++){
				sscanf(p,"%f",&(blosum->T[i][j]));
				p = strtok(NULL,", \t");
			}
		}

		blosum->adj = (int **) malloc((blosum->N) * sizeof(int *));
		for(i=0;i<blosum->N;i++){
			blosum->adj[i] = (int *) malloc(blosum->N * sizeof(int));
		}
		for(i=0;i<blosum->N;i++){
			fgets(line,maxl,fp);
			p=strtok(line,", \t");
			for(j=0;j<blosum->N;j++){
				sscanf(p,"%d",&(blosum->adj[i][j]));
				p = strtok(NULL,", \t");
			}
		}

		fclose(fp);
		return 0;
	}
	else
		return 60;
}


int stringPM::load_table(char *fn){


	//const int maxl=256;
	FILE *fp,*fp2;
	char fn2[maxl];
	char line[maxl];
	char label[maxl];
	int i,found,found2;


	if((fp=fopen(fn,"r"))!=NULL){
		found = 0;
		while((fgets(line,maxl,fp))!=NULL){
			memset(label,0,maxl);
			sscanf(line,"%s",label);
			//printf("line = %s",line);
			if(!strncmp(line,"USING",5)){
				sscanf(line,"%*s %s",fn2);
				strcpy(swt_fn,fn2);
				found = 1;
				break;
			}
			if(!strncmp(line,"SUBMAT",6)){
				sscanf(line,"%*s %s",fn2);
				strcpy(swt_fn,fn2);
				found = 2;
				break;
			}
		}
		fclose(fp);
		if(!found){
			printf("No MIS or MTX file name found - using default values\n");
			sprintf(swt_fn,"UsingDefaults!");
			found = 3;
		}


		switch(found){

		case 1:
			if((fp2=fopen(fn2,"r"))!=NULL){
				found2 = 0;
				//BUG: For some reason, when we open this again, we aren't at the beginning of the file..
				//can't understand why this is suddenly happening, but let's put a rewind in to fix it for now!
				rewind(fp2);
				while((fgets(line,maxl,fp2))!=NULL){
					memset(label,0,maxl);
					sscanf(line,"%s",label);
					//printf("line = %s",line);
					if(!strncmp(line,"SET",3)){
						sscanf(line,"%*s %s",fn2); //using fn2 to temporarily hold the alphabet...
						blosum->N = strlen(fn2);
						blosum->key = (char *)malloc((blosum->N+1) * sizeof(char));
						memset(blosum->key,0,(blosum->N+1)*sizeof(char));
						strncpy(blosum->key,fn2,blosum->N);

						blosum->T = (float **) malloc((blosum->N +1) * sizeof(float *));
						for(i=0;i<blosum->N+1;i++){
							blosum->T[i] = (float *) malloc(blosum->N * sizeof(float));
						}

						table_from_string(blosum->T,blosum->key,blosum->N);

						found2 = 1;
						fclose(fp2);

						//SUCCESS
						return 0;
					}
				}
				if(!found2){
					fclose(fp2);
					printf("No MIS string found in MIS file\n");
					return 2;
				}

			}
			else{
				printf("Unable to open MIS file %s",fn2);
				return 3;
			}
			break;
		case 2:
			return load_table_matrix(swt_fn);
			return 0;
			break;
		case 3:
			blosum =  default_table();
			return 0;
			break;
		}
		print_swt(stdout,blosum);
	}
	else
		return 4;

	return 5;


}




/**Suggest: load_agents(char *fn, char *fntab, int test, int verbose)
   if fntab == null: load_table(fn)
   else load_table_matrix(fntab)

*/
int stringPM::load_agents(char *fn, char *fntab, int test, int verbose){


	//testprop();

	int dmxl = maxl;
	int err = readordef_param_int(fn,"MAXLEN", &maxl, dmxl, 0);
	if(err>1){
		printf("ERROR %d on loading max line length (MAXLEN)\n",err);
		exit(0);
	}
	maxl0 = maxl+1;


	err = readordef_param_int(fn,"SPLPRINT", &splprint, splprint, 0);
	if(err>1){
		printf("ERROR %d on loading specieslist print frequency\n",err);
		exit(0);
	}
	maxl0 = maxl+1;


	const int llen = maxl0 + 256;
	FILE *fp;
	char line[llen];
	char label[llen];
	char code;
	char lab;
	int i,nag;
	s_ag *pag;
	int ntt = 0;



	int destep = estep;
	int estep_err = readordef_param_int(fn,"ESTEP", &estep, destep, 0);
	if(estep_err>1){
		printf("ERROR %d on loading energy per time step (ESTEP)\n",estep_err);
		exit(0);
	}


	//ALSO important to load blosum table here too!
	int table_err;
	if(fntab)
		table_err = load_table_matrix(fntab);
	else
		table_err = load_table(fn);
	if(table_err){
		printf("ERROR %d on loading table\n",table_err);
	}

	int mut_err = load_mut(fn,verbose);
	if(mut_err){
		printf("ERROR %d on loading mutation rate\n",mut_err);
	}

	int decay_err = load_decay(fn,verbose);
	if(decay_err){
		printf("ERROR %d on loading decay rate\n",decay_err);
	}

	ntt = 0;

	if((fp=fopen(fn,"r"))!=NULL){
		while((fgets(line,llen,fp))!=NULL){
			memset(label,0,llen);
			sscanf(line,"%s",label);
			//printf("line = %s",line);
			if(!strncmp(line,"AGENT",5)){
				ntt++;
			}
		}

		//Now set up the tracking arrays:
		//aat = (int *) mymalloc(ntt,sizeof(int));
		//aac = (int *) mymalloc(ntt,sizeof(int));
		//adc = (int *) mymalloc(ntt,sizeof(int));
		//aro = (int *) mymalloc(ntt,sizeof(int));
		//memset(adc,0,ntt*sizeof(int));
		//memset(aro,0,ntt*sizeof(int));
		rewind(fp);
		//int a=0;

		while((fgets(line,llen,fp))!=NULL){


			memset(label,0,llen);
			sscanf(line,"%s",label);
			//printf("line = %s",line);
			if(!strncmp(line,"AGENT",5)){



				memset(label,0,llen);
				sscanf(line,"%*s %s %d %c",label,&nag,&code);

				if(test){
					nag=1;
				}

				//make the agents
				for(i=0;i<nag;i++){
					l_spp *s;

					pag = make_ag(lab,1);

					pag->S =(char *) malloc(maxl0*sizeof(char));
					pag->label = code;
					memset(pag->S,0,maxl0*sizeof(char));
					strncpy(pag->S,label,strlen(label));
					pag->len = strlen(pag->S);


					//No parents for these initial agents!
					pag->pp = spl->make_parents(NULL,NULL);

					if(!i){
						//Record the species
						//s = spl->make_spp(pag);
						//s->tspp=0;
						//spl->prepend_spp(s); //append_lspp(s);

						//int stringPM::update_lineage(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp)
						update_lineage(pag,'I',1,NULL,NULL,0);
						s = spl->getspp(pag,extit,maxl0);
						//TODO: tidy up handling of seed species, but for now:
						s->tspp = 0;
					}

					//Record that the replicase copies itself
					pag->spp=s;
					//pag->paspp=pag->ppspp=1;


					//printf("agent %d, %c, %0.3f %0.3f\n",i,pag->label,pag->x,pag->y);
					if(nowhead == NULL){
						nowhead = pag;
					}
					else{
						append_ag(&nowhead,pag);
					}
				}


			}
		}
		//close
		fclose(fp);
		return 1;
	}
	else{
		printf("Unable to open file %s\n",fn);
		return 0;
	}
}



int stringPM::count_spp(){
	int i,found,finished;
	int nag, *done;
	s_ag *pa;
	int spc;
	nag=nagents(nowhead,-1);
	int sppcount=0;

	done = (int *) malloc(nag*sizeof(int));
	memset(done,0,nag*sizeof(int));

	do{
		i = 0;
		found=0;
		finished = 1;
		for(i=0,pa=nowhead;i<nag;i++,pa=pa->next){
			if(!done[i]){
				if(!found){
					//Set the species
					spc = pa->spp->spp;
					//Check it off the list
					done[i]=1;
					sppcount++;
					finished=0;
					found=1;
				}
				else{
					if(pa->spp->spp==spc){
						done[i]=1;
					}
				}
			}
		}
	}while(!finished);
	free(done);
	return sppcount;
}




void stringPM::print_spp_count(FILE *fp,int style, int state){

	//char fn[128];
	s_ag *pa;
	int finished = 0;
	int nag,*done;
	int i,found;
	int *spno, *spct, nspp = count_spp();

	nag = nagents(nowhead,-1);

	done = (int *) malloc(nag*sizeof(int));
	memset(done,0,nag*sizeof(int));

	spno = (int *) malloc(nspp*sizeof(int));
	memset(spno,0,nspp*sizeof(int));
	spct = (int *) malloc(nspp*sizeof(int));
	memset(spct,0,nspp*sizeof(int));
	int spidx=0;

	do{
		i = 0;
		found=0;
		finished = 1;
		for(i=0,pa=nowhead;i<nag;i++,pa=pa->next){
			if(state==-1 || pa->status == state){
				if(!done[i]){
					if(!found){
						done[i]=1;
						spno[spidx]=pa->spp->spp;
						spct[spidx]=1;
						finished=0;
						found=1;
					}
					else{
						if(pa->spp->spp==spno[spidx]){
							done[i]=1;
							spct[spidx]++;
						}
					}
				}
			}
		}
		spidx++;
	}while(!finished);

	int *spindx;
	spindx = (int *) malloc(nspp*sizeof(int));
	for(int ii=0;ii<nspp;ii++){
		spindx[ii]=ii;
	}

	if(nspp>1){

		//print_hsort_data(spno,spindx,nspp,stdout);
		idx_hsort_int(spno,spindx,nspp);
	}
	else{
		spindx[0]=0;
	}

	switch(style){
	case 0: //Horizontal
		fprintf(fp,"-------");
		for(i=0;i<nspp;i++){
			fprintf(fp,"-------");
		}
		fprintf(fp,"\n%d\t",(int) extit);
		for(i=0;i<nspp;i++){
			fprintf(fp,"%03d\t",spno[spindx[i]]);
		}
		fprintf(fp,"\n%d\t",(int) extit);
		for(i=0;i<nspp;i++){
			fprintf(fp,"%03d\t",spct[spindx[i]]);
		}
		fprintf(fp,"\n");
		break;
	case 1: //Vertical
		for(i=0;i<nspp;i++){
			fprintf(fp,"%d,%d,%d\n",(int) extit,spno[spindx[i]],spct[spindx[i]]);
		}
		break;
	}


	free(done);
	free(spno);
	free(spct);
	free(spindx);
}









void stringPM::get_spp_count(int state){

	//char fn[128];
	s_ag *pa;
	int nag,*done;
	int *spno, *spct, nspp = count_spp();

	nag = nagents(nowhead,-1);

	l_spp *pA;

	done = (int *) malloc(nag*sizeof(int));
	memset(done,0,nag*sizeof(int));

	spno = (int *) malloc(nspp*sizeof(int));
	memset(spno,0,nspp*sizeof(int));
	spct = (int *) malloc(nspp*sizeof(int));
	memset(spct,0,nspp*sizeof(int));

	//Set the spl->count param to zero
	for(pA=spl->species;pA!=NULL;pA=pA->next){
		pA->count = 0;
	}


	s_ag *p;
	for(p = nowhead;p!=NULL;p=p->next){
		if(state==-1 || pa->status == state){
			p->spp->count++;
		}
	}


}














int stringPM::countcomp(char E, char P){

	int count = 0;
	s_ag *pag;
	pag = nowhead;
	while(pag!=NULL){
		switch(pag->status){
		case B_ACTIVE:
			if(pag->label == E)
				if(pag->pass->label==P)
					count++;
			break;
		default:
			break;
		}
		pag = pag->next;
	}

	return count;


}


int stringPM::append_ag(s_ag **list, s_ag *ag){
	s_ag *pag;

	//printf("appending: list = %p, ag = %p\n",*list,ag);
	if(*list==NULL){
		*list=ag;
		//printf("appended: list = %p, ag = %p\n",*list,ag);
	}
	else{
		pag = *list;
		while(pag->next != NULL){
			pag = pag->next;
		}
		pag->next = ag;
		ag->prev = pag;
	}

	return 0;
}



int stringPM::extract_ag(s_ag **list, s_ag *ag){

	//printf("extracting: list = %p, ag = %p\n",*list,ag);
	if(ag == *list){
		*list = ag->next;
		if(*list !=NULL)
			(*list)->prev = NULL;
	}
	else{
		if(ag->prev==NULL)
			printf("No previous member of the list!\n");
		ag->prev->next = ag->next;
		if(ag->next != NULL)
			ag->next->prev = ag->prev;
	}
	ag->prev=NULL;
	ag->next=NULL;
	return 0;
}


int stringPM::hasdied(){

	if(nowhead == NULL)
		return 1;
	else
		return 0;

	/*
	int i;

	int sum=0;
	//NOTE: This is hard coded!
	for(i=0;i<ntt;i++)
		if(aat[i]=='P')
			sum += aac[i];

	if(!sum)
		return 1;

	return 0;
	*/
}



/*
void stringPM::move_ag(s_ag *a1){

	float rad = sqrt(pow(a1->x,2) + pow(a1->y,2));
	float phi;
	float edge = cellrad-(move*2);
	if (rad<edge){
		phi = rand0to1()*pi()*2;
	}
	else{
	    phi = atan2(a1->x,a1->y) + (pi()/2.);
	}


	a1->x = a1->x + (move * cos(phi));
	a1->y = a1->y - (move * sin(phi));

	//Now check if we've moved outside
	rad = sqrt(pow(a1->x,2) + pow(a1->y,2));

	if(rad>edge){
		float x,y;
		rand_in_rad(cellrad-(move*2),&x,&y);
		a1->x = x;
		a1->y = y;
	}
}
*/


int stringPM::free_ag(s_ag *pag){

	if(pag->S != NULL){
		//printf("destroying agent %d, code = %s\n",pag->idx,pag->S);
		free(pag->S);
	}

	free(pag);
	pag = NULL;

	return 0;
}

int stringPM::nagents(s_ag *head, int state){
	s_ag *pag;
	int count=0;
	pag = head;
	while(pag!=NULL){
		switch(state){
		case -1:
			count++;
			break;
		default:
			if(pag->status == state)
				count++;
				/* no break */
		}
		pag=pag->next;
	}
	return count;
}



s_ag * stringPM::rand_ag(s_ag *head, int state){
	int count = nagents(head,state);
	int i,pos;
	s_ag *pag,**arr;

	if(!count)
		return NULL;

	pag=NULL;

	switch(state){
	case -1:
		while(pag==NULL){
			pos = (int) (count * rand0to1());
			//printf("count = %d, pos = %d\n",count,pos);
			pag = head;
			for(i=0;i<pos;i++){
				pag = pag->next;
			}
		}
		break;
	case B_UNBOUND:
	case B_ACTIVE:
	case B_PASSIVE:
		arr = (s_ag **) malloc(count*sizeof(s_ag *));
		i=0;
		pag=head;
		while(pag!=NULL){
			if(pag->status==state)
				arr[i++]=pag;
			pag=pag->next;
		}
		pos=count;
		while(pos==count){
			pos = (int) (count * rand0to1());
		}
		pag=arr[pos];
		free(arr);
		break;
	default:
		pag = NULL;

	}
	return pag;
}




s_ag * stringPM::make_ag(int alab, int randpos){

	s_ag *ag;

	//printf("Spatial make_ag called\n");fflush(stdout);

	if((ag = (s_ag *) mymalloc(1,sizeof(s_ag)))!=NULL){
		ag->label=alab;
		ag->next = NULL;
		ag->prev = NULL;
		ag->exec = NULL;
		ag->pass = NULL;
		ag->S = NULL;
		ag->spp = NULL;
		ag->status = B_UNBOUND;
		ag->idx = agct++;ag->nbind=0;ag->ect=0;
		ag->biomass=0;
		return ag;
	}
	else{
		printf("mymalloc error\n");fflush(stdout);
		getchar();
		return NULL;
	}
}



//Diagnoistics / outputs
void stringPM::print_agents(FILE *fp, const char *spec, int verbose){

	s_ag *pag;

	if(!strncmp("NOW",spec,strlen("NOW"))){
		pag = nowhead;
	}

	if(!strncmp("NEXT",spec,strlen("NEXT"))){
		pag = nexthead;
	}

	while(pag!=NULL){
		if(verbose){
			switch(pag->status){
				case B_UNBOUND:
					fprintf(fp,"Agent %6d,\texec=%4d\tnbind=%3d\tUNBOUND, %s\n",pag->idx,pag->ect,pag->nbind,pag->S);
					break;
				case B_ACTIVE:
					print_exec(fp,pag,pag->pass);
					break;
				default:
					break;
				}
		}
		else{
			switch(pag->status){
				case B_UNBOUND:
					fprintf(fp,"Agent %6d,\texec=%4d\tnbind=%3d\tUNBOUND, %s\n",pag->idx,pag->ect,pag->nbind,pag->S);
					break;
				case B_ACTIVE:
					fprintf(fp,"Agent %6d, \texec=%4d\tnbind=%3d\tACTIVE,  %s\n",pag->idx,pag->ect,pag->nbind,pag->S);
					fprintf(fp,"Agent %6d, \texec=%4d\tnbind=%3d\tPASSIVE, %s\n",pag->pass->idx,pag->pass->ect,pag->pass->nbind,pag->pass->S);
					//print_exec(stdout,pag,pag->pass);
					break;
				default:
					break;
				}
		}
		pag=pag->next;
	}
}

int stringPM::print_agent_idx(FILE *fp, int det, int idx){
	s_ag *pag;
	pag = nowhead;
	while(pag!=NULL){
		if(pag->idx == idx)
			switch(pag->status){
			case B_UNBOUND:
				printf("Agent %6d,\texec=%4d\tnbind=%3d\tUNBOUND, %s\n",pag->idx,pag->ect,pag->nbind,pag->S);
				return 1;
			case B_ACTIVE:
				if(det)
					print_exec(stdout,pag,pag->pass);
				else{
					printf("Agent %6d, \texec=%4d\tnbind=%3d\tACTIVE,  %s\n",pag->idx,pag->ect,pag->nbind,pag->S);
					printf("Agent %6d, \texec=%4d\tnbind=%3d\tPASSIVE, %s\n",pag->pass->idx,pag->pass->ect,pag->pass->nbind,pag->pass->S);
				}
				//print_exec(stdout,pag,pag->pass);
				return 1;
			case B_PASSIVE:
				if(det){
					print_exec(stdout,pag->exec,pag);
				}
				else{
					printf("Agent %6d, \texec=%4d\tnbind=%3d\tACTIVE,  %s\n",pag->exec->idx,pag->exec->ect,pag->exec->nbind,pag->exec->S);
					printf("Agent %6d, \texec=%4d\tnbind=%3d\tPASSIVE, %s\n",pag->idx,pag->ect,pag->nbind,pag->S);
				}
				//print_exec(stdout,pag,pag->pass);
				return 1;
			}

		pag=pag->next;
	}
	return 0;
}

/*
float stringPM::close_dist(s_ag *a1, s_ag *a2){

	return sqrt(pow(a1->x-a2->x,2) + pow(a1->y-a2->y,2));
}
*/


/*
void stringPM::update_aac(){
	int a;
	s_ag *ag;
	ag = nowhead;
	memset(aac,0,ntt*sizeof(int));
	while(ag!=NULL){
		a=0;
		while(aat[a]!=ag->label && a <= ntt){
			a++;
			}
		if(a<ntt){
			aac[a]++;
		}
		ag=ag->next;
	}
}*/


int stringPM::aspatial_coverage(int n){

	return eqn_prop(n);
}


int stringPM::aspbind(rules *rset, s_ag *pag){

	int found=0;
	int count=0;
	int
		frule,
		rule = -1;
	float
		rno;
	s_ag *bag,*wbag,*pag2;

	bag=nowhead;

	while(bag!=NULL){
		//if Partner
		if((rule = rset->getrule("lhs",bag->label,pag->label))>-1){
			wbag = bag;
			count++;
		}
		bag=bag->next;
	}

	if(aspatial_coverage(count)){
		found = 1;
	}
	note_propensity(count,found);

	float prno = rand0to1();
	int pr2=0,
		r2 = floor(count * 0.9999 * prno);
	if(found){
		frule = -1;
		pag2 = nowhead;
		while(pag2!=NULL){

			if((rule = rset->getrule("lhs",pag2->label,pag->label))>-1){
				if(pr2==r2){
					frule = rule;
					break;
				}
				pr2++;
			}

			pag2=pag2->next;
		}

		if(frule ==-1){
			printf("Fucked up!\n");
		}
		rno = rand0to1();
		tr[frule]++;

		//printf("rno = %0.6f\tfrp=%0.6f\t",rno,rset->rval[frule] );
		if(rno<rset->rval[frule] && energy >= 3*(-1*rset->re[frule])){
			//change label and move pag:
			//printf("Bonding %c to %c now, making %c\n",pag->label,pag2->label,rset->rset[frule][2]);
			//fflush(stdout);

			//printf("Fired\n");

			//record the firing
			fr[frule]++;

			energy += rset->re[frule];

			//pag->label = rset->rset[frule][2];
			append_ag(&nexthead,pag);

			//delete the second agent;
			extract_ag(&nowhead,pag2);
			free_ag(pag2);
		}
		else{
			found=0;
			//printf("Didn't fire\n");
		}
	}
	return found;
}


int stringPM::testdecay(s_ag *pag){

#ifdef LONG_DECAY
	//VARIABLE DECAY RATE BASED ON LENGTH OF STRING (ECAL 2009)
	float len = strlen(pag->S);
	float prob = 1./pow(len,2);//4./3.);
#else
	//CONSTANT DECAY RATE to match ECAL (ALife 2010 and on)
 	float prob = decayrate;//1./pow(65,2);//4./3.); //This is now done in load_decay...
#endif


	float rno = rand0to1();

#ifdef UNB_DECAY_ONLY
	if(rno<prob && pag->status == B_UNBOUND){
#else
	if(rno<prob && dodecay){
#endif
		//unbind_ag(pag);
		free_ag(pag);
		return 1;
	}
	else
		return 0;
}

void stringPM::get_string_comp(s_ag *pag){
/*
	if(pag->comp != NULL){
		free(pag->comp);
		pag->comp = NULL;
	}
	pag->comp = string_comp(pag->S);
	*/
}



float stringPM::sigalign(char *str){

	int L_sig = strlen(signal);
	int L_str = strlen(str);
	align sw;
	char *comp;
	float bprob;

	comp = string_comp(str);
	bprob = SmithWatermanV2(comp,signal,&sw,blosum,0);
	free(comp);
	int l = L_str>L_sig?L_str:L_sig;
	int la = sw.e1-sw.s1 < sw.e2-sw.s2 ? sw.e1-sw.s1 : sw.e2-sw.s2;
	if(la<=2)
		bprob=0;
	else{
		float s = sw.score<l-1.124? sw.score : l-1.124;
		bprob = s/(l-1.124);
	}

	return bprob;

}



float stringPM::get_bprob(align *sw){
	float bprob = 0.;
	//This is the old bind prob, with a modifier for short strings:

	int l = sw->e1-sw->s1 < sw->e2-sw->s2 ? sw->e1-sw->s1 : sw->e2-sw->s2;
	if(l<=2)
		bprob=0;
	else{
		//bprob = pow(sw->score,l)/pow(l,l);
		//BRUTAL HACK:
		float s = sw->score<l-1.124? sw->score : l-1.124;
		bprob = s/(l-1.124);
	}

	//Here is the new bind prob:

	//float sless1 = sw->score -1.;
	//if(sless1<0){
	//	printf("Negative score encountered!\n");
	//	sless1=0;
	//}
	//bprob = pow(sless1/l,1/l);

	return bprob;
}


float stringPM::get_sw(s_ag *a1, s_ag *a2, align *sw){

	float bprob;
	char *comp;
	s_sw *swa;

	//SUGGEST: pass in pointer to the species - not its index
	swa = read_sw(swlist,a1->spp->spp,a2->spp->spp);

	if(swa==NULL){

		//get_string_comp(a1);
		comp = string_comp(a1->S);

		//

		bprob = SmithWatermanV2(comp,a2->S,sw,blosum,0);
		//bprob = SmithWaterman(comp,a2->S,sw,blosum,0);

		free(comp);

		align sw2;

		bprob = SmithWatermanV2(a1->S,a2->S,&sw2,blosum,0);


		int lA = strlen(a1->S);
		int lB = strlen(a2->S);
		int L = lA>lB?lA:lB;

		float spp_sim=sw2.score/L;


		//SUGGEST: pass in pointer to the species - not its index
		store_sw(&swlist,sw,a1->spp->spp,a2->spp->spp);
	}
	else{
		load_sw(swa,sw);
	}

	bprob = get_bprob(sw);


	if(verbose_bind){
		printf("Alignment:\nm1: %d to %d\nm2: %d to %d\nscore = %f\nProb = %f = %E\n",sw->s1,sw->e1,sw->s2,sw->e2,sw->score,bprob,bprob);
	}


	return bprob;
}


void stringPM::print_ptr_offset(FILE *fp, char *S, char *p,int F, char c){
	int i,n=p-S;
	if(n<0){
		printf("Problem calculating pointer location\n");
		fflush(stdout);
	}
	for(i=0;i<n;i++)
		fprintf(fp," ");
	fprintf(fp,"%c\n",F?c-32:c);
}


void stringPM::print_exec(FILE *fp, s_ag *act, s_ag *pas){

	//Diagnostics to screen for passive:
	if(!strlen(pas->S))
		printf("Zero length passive string\n");
	fprintf(fp,"%6d:\n%s\n",pas->idx,pas->S);
	print_ptr_offset(fp,pas->S,act->i[0],1-act->it,'i');
	print_ptr_offset(fp,pas->S,act->f[0],1-act->ft,'f');
	print_ptr_offset(fp,pas->S,act->r[0],1-act->rt,'r');
	print_ptr_offset(fp,pas->S,act->w[0],1-act->wt,'w');

	//Diagnostics to screen for active:
	if(!strlen(act->S))
		printf("Zero length active string\n");
	fprintf(fp,"%6d:\n%s\n",act->idx,act->S);
	print_ptr_offset(fp,act->S,act->i[1],act->it,'i');
	print_ptr_offset(fp,act->S,act->f[1],act->ft,'f');
	print_ptr_offset(fp,act->S,act->r[1],act->rt,'r');
	print_ptr_offset(fp,act->S,act->w[1],act->wt,'w');

	act->len = strlen(act->S);
	pas->len = strlen(pas->S);


	if(pas->len<=maxl){
		//printf("Passive string length = %d\n",pas->len);
	}
	else
		printf("Passive string length = %d - TOO LONG\n",pas->len);

	if(act->len<=maxl){
		//printf("Active  string length = %d\n",act->len);
	}
	else
		printf("Active  string length = %d - TOO LONG\n",act->len);


}



void stringPM::set_exec(s_ag *A, s_ag *B, align *sw){

	s_ag *active,*passive;
	int active_idx,passive_idx;

	if(sw->s1>=sw->s2){
		active = A;
		passive = B;
		active_idx = sw->s1;
		passive_idx = sw->s2;
	}
	else{
		active = B;
		passive = A;
		active_idx = sw->s2;
		passive_idx = sw->s1;
	}

	active->status = B_ACTIVE;
	active->exec = NULL;
	active->pass = passive;

	active->biomass = 0;

	passive->status = B_PASSIVE;
	passive->exec = active;
	passive->pass = NULL;



	active->f[0] = active->i[0] = active->r[0] = active->w[0] = &(passive->S[passive_idx]);//&(passive->S[sw->s1]);
	active->f[1] = active->i[1] = active->r[1] = active->w[1] = &(active->S[active_idx]);//&(active->S[sw->s2]);
	active->ft   = active->it   = active->rt   = active->wt = 1;


	passive->f[0] = passive->i[0] = passive->r[0] = passive->w[0] = 0;
	passive->f[1] = passive->i[1] = passive->r[1] = passive->w[1] = 0;
	passive->ft   = passive->it   = passive->rt   = passive->wt = 0;




#ifdef V_VERBOSE
	printf("Bind finished - looks like:\n");
	print_exec(stdout,active,passive);
#endif
}


int stringPM::unbind_ag(s_ag * pag, char sptype, int update, l_spp *pa, l_spp *pp){

	int found;
	int mass=0;

	if(pag->status==B_ACTIVE){
		mass = pag->biomass;
		pag->biomass = 0;
	}


	pag->status = B_UNBOUND;
	pag->pass = NULL;
	pag->exec = NULL;

	pag->ect=0;

	pag->f[0] = pag->i[0] = pag->r[0] = pag->w[0] = 0;
	pag->f[1] = pag->i[1] = pag->r[1] = pag->w[1] = 0;
	pag->ft   = pag->it   = pag->rt   = pag->wt = 0;

	found = update_lineage(pag,sptype,update,pa,pp,mass);
	return found;
}


int stringPM::testbind(s_ag *pag){

	int found=0;
	int count=0;
	float
		rno;
	s_ag *bag;
	align sw;
	float bprob =0.;

	bag=nowhead;

	count = nagents(nowhead,B_UNBOUND);

	//Calc propensity
	if(eqn_prop(count)){
		found = 1;
	}

	if(found){
		bag = rand_ag(nowhead,B_UNBOUND);
#ifndef BIND_ALL
		bprob = get_sw(pag,bag,&sw);
#else
		bprob =1.0;
		sw.match = 1;		// the number of matching characters.
		sw.score = 1; 	// the score of the match
		sw.prob = 1.0;		// the probability of the match - used for determining events based on the score/match
		sw.s1=0;			// start of the match in string 1
		sw.e1=1;			// end of the match in string 1
		sw.s2=0;			// start of the match in string 2
		sw.e2=1;			// end of the match in string 2
#endif
		rno = rand0to1();
		if(rno<bprob){//Binding success!
			//figure out which is the executing string:
			set_exec(pag,bag,&sw);
			pag->nbind++;
			bag->nbind++;

			append_ag(&nexthead,pag);

			//uptime the second agent;
			extract_ag(&nowhead,bag);
			append_ag(&nexthead,bag);
			energy--;
		}
		else{
			found = 0;
		}
	}
	note_propensity(count,found);
	return found;
}





int stringPM::h_pos(s_ag *pag, char head){

	char *ph;
	char *ps;

	if(pag->status != B_ACTIVE)
		printf("ERROR: attempting head position for inactive string");

	switch(head){
	case 'w':
		ph = pag->w[pag->wt];
		if(pag->wt)
			ps = pag->S;
		else
			ps = pag->pass->S;
		break;
	case 'f':
		ph = pag->f[pag->ft];
		if(pag->ft)
			ps = pag->S;
		else
			ps = pag->pass->S;
		break;
	case 'i':
		ph = pag->i[pag->it];
		if(pag->it)
			ps = pag->S;
		else
			ps = pag->pass->S;
		break;
	case 'r':
		ph = pag->r[pag->rt];
		if(pag->rt)
			ps = pag->S;
		else
			ps = pag->pass->S;
		break;
	}

	return ph-ps;

}


int stringPM::rewind_bad_ptrs(s_ag* act){

	int plen,alen,pdist;
	char *ps;

	//PUT DANGLING POINTERS AT THE *END* OF THE STRINGS:
	//DO THE PASSIVE POINTERS FIRST:
	plen = strlen(act->pass->S);
	if(plen){
		ps = act->pass->S;

		pdist = act->i[0]-ps;
		if(pdist>plen || pdist<0)
			act->i[0]=ps+plen;

		pdist = act->r[0]-ps;
		if(pdist>plen || pdist<0)
			act->r[0]=ps+plen;

		pdist = act->w[0]-ps;
		if(pdist>plen || pdist<0)
			act->w[0]=ps+plen;

		pdist = act->f[0]-ps;
		if(pdist>plen || pdist<0)
			act->f[0]=ps+plen;
	}
	else{//Toggle everything off this string...
		act->i[0]=act->pass->S;
		act->r[0]=act->pass->S;
		act->w[0]=act->pass->S;
		act->f[0]=act->pass->S;

		act->it=1;
		act->rt=1;
		act->wt=1;
		act->ft=1;
	}


	//DO THE ACTIVE POINTERS NOW
	alen = strlen(act->S);
	ps = act->S;

	if(alen){
		pdist = act->i[1]-ps;
		if(pdist>alen || pdist<0)
			act->i[1]=ps+alen;

		pdist = act->r[1]-ps;
		if(pdist>alen || pdist<0)
			act->r[1]=ps+alen;

		pdist = act->w[1]-ps;
		if(pdist>alen || pdist<0)
			act->w[1]=ps+alen;

		pdist = act->f[1]-ps;
		if(pdist>alen || pdist<0)
			act->f[1]=ps+alen;
	}
	else{//Toggle everything off this string...
		act->i[1]=act->S;
		act->r[1]=act->S;
		act->w[1]=act->S;
		act->f[1]=act->S;

		if(plen){
			act->it=0;
			act->rt=0;
			act->wt=0;
			act->ft=0;
		}
	}

	//TODO: error checking on this!
	return 0;
}


int stringPM::check_ptrs(s_ag* act){

#ifdef VERBOSE
	print_exec(stdout,act,act->pass);
#endif
	//int len,pdist;
	//char *ps;

	//Sort pointers out first - even if there's going to be an error!
	rewind_bad_ptrs(act);

	//Step 1: make sure act and pass *have* strings...
	if(!strlen(act->S)){
#ifdef VERBOSE
		printf("Zero length active string - dissoc\n");
#endif
		return 1;
		//if(!strlen(act->pass->S)){
		//	printf("Zero length active and passive strings - destroy\n");
		//	return 3;
		//}
	}
	if(!strlen(act->pass->S)){
#ifdef VERBOSE
		printf("Zero length passive string - dissoc\n");
#endif
		return 2;
	}



#ifdef VERBOSE
	print_exec(stdout,act,act->pass);
#endif

	return 0;

}


int stringPM::hcopy(s_ag *act){

	//s_ag *pass;
	//pass = act->pass;
	int cidx;
	float rno;
	int safe = 1;// this gets set to zero if any of the tests fail..
	int doprint = 0;
	e_mut mut = M_NONE;


	//MUTATION RATES:
	//THESE ARE HARD-CODED FOR NOW - THEY SHOULD BE DERIVED FROM THE BLOSUM SOMEHOW...
	//const float indelrate = 0.0005,subrate=0.375;//0.0749/2
	//const float indelrate = 0.0000306125,subrate=0.01;//0.02
	//const float indelrate = 0.00006125,subrate=0.05;//0.02
	//const float indelrate = 0.000125,subrate=0.1;//0.02
	//const float indelrate = 0.000005,subrate=0.0375;//0.0749/2

	//const float indelrate = 0., 		subrate = 0.;

	if(!domut){
		indelrate = subrate =0;
	}

	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);

	int p;
	if( (p = h_pos(act,'w'))>=maxl){
		printf("Write head out of bounds: %d\n",p);
		//just to make sure no damage is done:
		if(act->wt)
			act->S[maxl]='\0';
		else
			act->pass->S[maxl]='\0';

		act->i[act->it]++;

		safe = 0;
		return -1;
	}
	if(h_pos(act,'r')>=maxl){
		printf("Read head out of bounds\n");

		act->i[act->it]++;

		safe = 0;
		return -2;
	}




	if(*(act->r[act->rt]) == 0){
		//possibly return a negative value and initiate a b
		safe = 0;
		//return -3;
	}
	//if(h_pos(act,'w')>=maxl){
	//	//We are at the end of a copy...
	//	//...so just increment *R
	//	act->r[act->rt]++;
	//	safe = 0;
	//}

	if(safe){
		rno=rand0to1();
		if(rno<indelrate){//INDEL
			doprint=1;
			//should follow the blosum table for this....
			rno=rand0to1();
			if(rno<0.5){//insert
				mut = M_INSERT;
				//first do a straight copy..
				*(act->w[act->wt])=*(act->r[act->rt]);

				//no need to test for granular here since we are inserting...
				act->w[act->wt]++;

				//Then pick a random instruction:
				cidx = (float) rand0to1() * blosum->N;

				//insert the random instruction
				*(act->w[act->wt])=blosum->key[cidx];
				if(granular_1==0){
					act->w[act->wt]++;
				}
			}
			else{//delete



				act->i[act->it]++;

				mut = M_DELETE;
				//simply increment the read head without doing anything else
			}

			if(granular_1==0){
				act->r[act->rt]++;
			}
			//act->i[act->it]++;
		}
		else{
			if(rno<subrate+indelrate){//INCREMENTAL MUTATION
				doprint=1;
				cidx = sym_from_adj(*(act->r[act->rt]),blosum);
				*(act->w[act->wt])=cidx;
			}
			else{//NO MUTATION
				*(act->w[act->wt])=*(act->r[act->rt]);
			}
			if(granular_1==0){
				act->w[act->wt]++;
				act->r[act->rt]++;
			}
		}
	}
	//update lengths
	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);


	act->i[act->it]++;


#ifdef VERBOSE
	if(mut)
	printf("Mutant event %d. new string is:\n%s\n\n",mut,act->wt?act->S:act->pass->S);
#endif
	act->biomass++;
	biomass++;
	return 0;
}




int stringPM::cleave(s_ag *act){

	int dac = 0,cpy;
	s_ag *c,*pass,*csite;

	pass = act->pass;

	//pick the mol containing the cleave site:
	csite = act->ft?act:pass;

	if(act->f[act->ft]-csite->S < csite->len){

		//1: MAKE THE NEW MOLECULE FROM THE CLEAVE POINT

		//Can't really say what the label is easily - for ECAL, it's always pass
		c = make_ag(pass->label,1);

		//Copy the cleaved string to the agent
		char *cs;

		c->S =(char *) malloc(maxl0*sizeof(char));
		memset(c->S,0,maxl0*sizeof(char));

		cs = csite->S;
		cpy = strlen(cs);
		//Check that we aren't creating a zero-length molecule:
		if(!cpy){
			printf("WARNING: Zero length molecule being created!\n");
		}

		//Make the parent structure: ALL DONE NOW IN update_lineage
		//c->pp = splist->make_parents(act->spp,pass->spp);

		cpy -= act->f[act->ft]-cs;

		strncpy(c->S,act->f[act->ft],cpy);
		c->len = strlen(c->S);
#ifdef VERBOSE
		printf("String %d created:\n%s\n",c->idx,c->S);
#endif

		//ALL DONE NOW IN update_lineage
		//Fill in the birth certificate:
		//Parents now set in update_lineage
		//c->paspp = act->spp;
		//c->ppspp = pass->spp;

		//Check the lineage
		update_lineage(c,'C',1,act->spp,pass->spp,act->biomass);
		act->biomass=0; //reset this; we might continue to make stuff!

		//append the agent to nexthead
		append_ag(&nexthead,c);

		//2: HEAL THE PARENT

		memset(act->f[act->ft],0,cpy*sizeof(char));

		csite->len = strlen(csite->S);

		if((dac = check_ptrs(act))){
			switch(dac){
			case 1://Destroy active - only append passive
				unbind_ag(pass,'P',1,act->spp,pass->spp);
				append_ag(&nexthead,pass);
				free_ag(act);
				break;
			case 2://Destroy passive - only append active
				unbind_ag(act,'A',1,act->spp,pass->spp);
				append_ag(&nexthead,act);
				free_ag(pass);
				break;
			case 3://Destroy both
				printf("This should never happen\n");
				unbind_ag(act,'A',1,act->spp,pass->spp);
				unbind_ag(pass,'P',1,act->spp,pass->spp);
				free_ag(act);
				free_ag(pass);
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
		//dac=-1;
	}

	return dac;
}




int stringPM::exec_step(s_ag *act, s_ag *pass){//pset *p,char *s1, swt *T){

	int finished=0;
	char *tmp;
	int dac=0;
	int safe_append=1;



	//if(act->idx==1222 || pass->idx==1222){
	//	print_exec(stdout,act,act->pass);
	//	fflush(stdout);
	//}

	//The following is unstable - explosions/catastrophe - possibly because of the nature of the decay??
	//Now moved energy to calling function
	//if(eqn_prop(energy)){
	//if(energy){

	switch(*(act->i[act->it])){//*iptr[it]){

	case '$'://h-search
		//act->ft = act->it;
		char *cs;
		if(act->ft)
			cs = act->S;
		else
			cs = act->pass->S;
		tmp = HSearch(act->i[act->it],cs,blosum,&(act->it),&(act->ft),maxl);
		act->f[act->ft] = tmp;
		act->i[act->it]++;
		break;

		/*
		//put the active flow head on the string where the active ip is:
		act->ft = act->it;
		char *cs;
		if(act->it)
			cs = act->S;
		else
			cs = act->pass->S;
		act->f[act->ft] = HSearch(act->i[act->it],cs,blosum);
		act->i[act->it]++;
		break;
		*/

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
		if(hcopy(act)<0){
			unbind_ag(act,'A',1,act->spp,pass->spp);
			unbind_ag(pass,'P',1,act->spp,pass->spp);
			finished = 1;
		}
		break;


	/************
	 *   INC_R  *
	 ************/
	case '+'://h-copy
		if(granular_1==1){
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
			act->i[act->it]=IfLabel(act->i[act->it],act->r[act->rt],act->S,blosum,maxl);
			break;


	/************
	 *  CLEAVE  *
	 ************/
	case '%':
			if((dac = cleave(act))){
				//unbind_ag(act);
				//unbind_ag(pass);
				finished = 1;
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
			unbind_ag(act,'A',1,act->spp,pass->spp);
			unbind_ag(pass,'P',1,act->spp,pass->spp);

			finished = 1;
			break;

	default://Just increment the i-pointer
		act->i[act->it]++;
		break;
	}
#ifdef V_VERBOSE
	printf("Exec step - looks like:\n");
	print_exec(stdout,act,pass);
#endif


	if(safe_append){
		act->ect++;
		append_ag(&nexthead,act);
		append_ag(&nexthead,pass);
	}
	energy--;



	return 1;
	//return finished;
}











/* HOW make_next WORKS:
 *
 * -While there's some agents left in "now"
 *      -pick an agent at random
 *      -if there are any rules applying to the agent
 *          -find the maximum number of rules on the lhs of all rules applying to the agent = maxl
 *          switch(maxl)
 *          -case 2: We have a bonding rule to test
 *          	-go through all other agents.
 *          		-If closer than the current closest agent && there's a bonding rule for the two
 *						-set 2nd agent to agent b.
 *						-update closest
 *				- if we have a reactive pair
 *					-if we've enough energy and a random number is less than reaction prob
 *						-do the reaction and put the product in the 'next' bucket
 * 			-case 1: Dissocn or decay rule to test (assuming only one of these rules)
 *				- get the appropriate rule and find if its dissoc or decay
 *				- if dissoc
 *					- do the dissoc
 *				- else
 *					- do the decay
 *		- if no reaction has been carried out
 *			- move the agent to the next bucket
 *
 */
//A make_next that doesn't use the "rules" class
void stringPM::make_next(){
	s_ag *pag,*bag;
	s_bind bb;
	int changed;

	//SUGGEST: write function to count what's around (saves rechecking every time)
	//countstates();

	//SUGGEST: IF WE ARE PUTTING ENERGY IN FROM OUTSIDE, AND NO ENERGY CAN BE PRODUCED
	//while(nowhead!=NULL && energy){
	//}
	//if(nowhead != NULL){
	//	pag = nowhead;
	//	append_ag(pag); //check that pag->next is also appended
	// NB: will have to have a different scheme for decay, since that doesn't require energy - sample from a binomial?

	while(nowhead!=NULL){

		pag = rand_ag(nowhead,-1);
		extract_ag(&nowhead,pag);
		changed = 0;

		//if(0){//extit>=10){
		//	if(pag->status==B_PASSIVE)
		//		//if(pag->exec->idx==1041){
		//			print_exec(stdout,pag->exec,pag);
		//		//}
		//	if(pag->status==B_ACTIVE){
		//		//if(pag->idx==1041)
		//			print_exec(stdout,pag,pag->pass);
		//	}
		//	fflush(stdout);
		//}

		//extract any partner:
		bag = NULL;
		bb = pag->status;
		switch(pag->status){
		case B_UNBOUND:
			break;
		case B_ACTIVE:
			//if(pag->idx == 214 && extit > 79)
			//	printf("Active 214\n");
			bag = pag->pass;
			extract_ag(&nowhead,bag);
			break;
		case B_PASSIVE:
			bag = pag->exec;
			extract_ag(&nowhead,bag);
			break;
		}

		if(pag->label=='P' && pag->status != B_UNBOUND  ){
			if(print_agent_idx(stdout,1,pag->idx))
				fflush(stdout);
		}

		/*
		if(pag->idx == 354 || pag->idx == 214){
			printf("%d ",(int) extit);
			if(print_agent_idx(stdout,1,354))
				fflush(stdout);
			if(print_agent_idx(stdout,1,214))
				fflush(stdout);
		}*/

		int dc = testdecay(pag);
		if(dc){//we must check what else needs to be destroyed...
			if(bag!=NULL){
				free_ag(bag);
			}
		}
		else{
			if(energy>0){
				switch(pag->status){
				case B_UNBOUND:
					//seek binding partner, set binding states.
					changed = testbind(pag);
					//if(!changed)
					//	changed = testdecay(pag);

					break;
				case B_PASSIVE:

					//extract_ag(&nowhead,pag->exec);
					changed = exec_step(pag->exec,pag);

					break;
				case B_ACTIVE:

					//extract_ag(&nowhead,pag->pass);
					changed = exec_step(pag,pag->pass);
					break;
				default:
					printf("ERROR: agent with unknown state encountered!\n");
				}
			}
			if(!changed){
				append_ag(&nexthead,pag);
				if(bag!=NULL)
					append_ag(&nexthead,bag);

			}
		}

		/*
		//Infinite loop at low nos - energy ALWAYS available - therefore always changed, therefore no decay!
		if(!changed){
			//if(changed)
			//	printf("Changed!");

			switch(pag->status){
			case B_UNBOUND:
				//if(pag->status == B_UNBOUND)
				changed = testdecay(pag);
				if(!changed)
					append_ag(&nexthead,pag);
				break;
			case B_PASSIVE:
				bag = pag->exec;
				changed = testdecay(pag);
				if(!changed){
					append_ag(&nexthead,pag->exec);
					append_ag(&nexthead,pag);
				}
				else{
					free_ag(bag);
					//unbind_ag(bag);
					//append_ag(&nexthead,bag);
				}
				break;
			case B_ACTIVE:
				bag = pag->pass;
				changed = testdecay(pag);
				if(!changed){
					append_ag(&nexthead,pag->pass);
					append_ag(&nexthead,pag);
				}
				else{
					free_ag(bag);
					//unbind_ag(bag);
					//append_ag(&nexthead,bag);
				}
				break;
			}
		}
		*/

		//if(pag->status == B_ACTIVE && pag->idx == 669){
		//	print_exec(stdout,pag,pag->pas);
		//}

		/*
		//For debugging:
		if(extit==18574){
			if(pag->status==B_PASSIVE)
				if(pag->exec->idx==668){
					print_exec(stdout,pag->exec,pag);
					print_status(stdout,pag->status);
				}
			print_status(stdout,pag->status);
			printf("ag = %d, ",pag->idx);
			print_status(stdout,pag->status);
			printf("\n");
			print_status(stdout,pag->status);
			//print_agent_idx(pag->idx);
			printf(" ");
			print_agent_idx(stdout,1,1222);
			printf(" \n");
			fflush(stdout);
		}*/

	}
}




int stringPM::speig_hcopy(s_ag *act){

	//s_ag *pass;
	//pass = act->pass;
	int cidx;
	float rno;
	int safe = 1;// this gets set to zero if any of the tests fail..
	int doprint = 0;
	e_mut mut = M_NONE;

	if(!domut){
		indelrate = subrate =0;
	}

	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);

	int p;
	if( (p = h_pos(act,'w'))>=maxl){
		printf("Write head out of bounds: %d\n",p);
		//just to make sure no damage is done:
		if(act->wt)
			act->S[maxl]='\0';
		else
			act->pass->S[maxl]='\0';

		act->i[act->it]++;
		safe = 0;
		return -1;
	}
	if(h_pos(act,'r')>=maxl){
		printf("Read head out of bounds\n");
		act->i[act->it]++;
		safe = 0;
		return -2;
	}




	if(*(act->r[act->rt]) == 0){
		//possibly return a negative value and initiate a b
		safe = 0;
		//return -3;
	}

	if(safe){

		//see if we are overwriting or not:
		int rm,wm=-1;
		if(*(act->w[act->wt])){
			wm=tab_idx(*(act->w[act->wt]),blosum);
		}

		const float speig_idrate = 0.001;
		float winc=rand0to1();
		float rinc=rand0to1();

		//todo: make sure no increments happen if the symbol (or mutant) is not available

		rno=rand0to1();
		if(rno<subrate){//INCREMENTAL MUTATION
			doprint=1;
			cidx = sym_from_adj(*(act->r[act->rt]),blosum);
			rm = tab_idx(cidx,blosum);
			if(mass[rm]){
				*(act->w[act->wt])=cidx;
				if(winc>speig_idrate)
					act->w[act->wt]++;

				if(rinc>speig_idrate)
					act->r[act->rt]++;//possible deletion here...

				if(!(wm<0)){
					mass[wm]++;
				}
				mass[rm]--;
			}
		}
		else{//NO MUTATION (but possible sub via comass effects)
			//cidx = sym_from_adj(*(act->r[act->rt]),blosum);
			rm = tab_idx(*(act->r[act->rt]),blosum);
			if(mass[rm]){
				*(act->w[act->wt])=*(act->r[act->rt]);

				if(winc>speig_idrate)
					act->w[act->wt]++;

				if(rinc>speig_idrate)
					act->r[act->rt]++;

				if(!(wm<0)){
					mass[wm]++;
				}
				mass[rm]--;
			}
			else{
				cidx = sym_from_adj(*(act->r[act->rt]),blosum);
				rm = tab_idx(cidx,blosum);
				if(mass[rm]){
					if(winc>speig_idrate)
						act->w[act->wt]++;

					if(rinc>speig_idrate)
						act->r[act->rt]++;

					act->w[act->wt]++;
					if(!(wm<0)){
						mass[wm]++;
					}
					mass[rm]--;
				}
			}
		}

	}
	//update lengths
	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);
	act->i[act->it]++;

#ifdef VERBOSE
	if(mut)
	printf("Mutant event %d. new string is:\n%s\n\n",mut,act->wt?act->S:act->pass->S);
#endif
	act->biomass++;
	biomass++;
	return 0;
}


/////////////////////////////////////////////////////////start of comass stuff



int stringPM::comass_hcopy(s_ag *act){

	//s_ag *pass;
	//pass = act->pass;
	int cidx;
	float rno;
	int safe = 1;// this gets set to zero if any of the tests fail..
	int doprint = 0;
	e_mut mut = M_NONE;


	//MUTATION RATES:
	//THESE ARE HARD-CODED FOR NOW - THEY SHOULD BE DERIVED FROM THE BLOSUM SOMEHOW...
	//const float indelrate = 0.0005,subrate=0.375;//0.0749/2
	//const float indelrate = 0.0000306125,subrate=0.01;//0.02
	//const float indelrate = 0.00006125,subrate=0.05;//0.02
	//const float indelrate = 0.000125,subrate=0.1;//0.02
	//const float indelrate = 0.000005,subrate=0.0375;//0.0749/2

	//const float indelrate = 0., 		subrate = 0.;

	if(!domut){
		indelrate = subrate =0;
	}
	//Make sure the recorded lengths are current
	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);

	int p;

	//Check positions of the pointers
	if( (p = h_pos(act,'w'))>=maxl){
		printf("Write head out of bounds: %d\n",p);
		//just to make sure no damage is done:
		if(act->wt)
			act->S[maxl]='\0';
		else
			act->pass->S[maxl]='\0';

		act->i[act->it]++;
		safe = 0;
		return -1;
	}
	if(h_pos(act,'r')>=maxl){
		printf("Read head out of bounds\n");
		act->i[act->it]++;
		safe = 0;
		return -2;
	}

	//Check that we aren't off the end of the string, but within the allocated memory:
	if(*(act->r[act->rt]) == 0){
		//possibly return a negative value and initiate a b
		safe = 0;
		//return -3;
	}
	//if(h_pos(act,'w')>=maxl){
	//	//We are at the end of a copy...
	//	//...so just increment *R
	//	act->r[act->rt]++;
	//	safe = 0;
	//}

	if(safe){

		//see if we are overwriting or not:
		int rm,wm=-1;
		if(*(act->w[act->wt])){
			wm=tab_idx(*(act->w[act->wt]),blosum);
		}


		rno=rand0to1();
		if(rno<indelrate){//INDEL - we should never be doing this in comass!
			doprint=1;
			//should follow the blosum table for this....
			rno=rand0to1();
			if(rno<0.5){//insert
				mut = M_INSERT;
				//first do a straight copy..
				*(act->w[act->wt])=*(act->r[act->rt]);
				act->w[act->wt]++;

				//Then pick a random instruction:
				rm = (float) rand0to1() * blosum->N;

				//Check there's mass for this symbol:
				if(mass[rm]){
					//insert the random instruction
					*(act->w[act->wt])=blosum->key[rm];
					if(!(wm<0)){
						mass[wm]++;
					}
					mass[rm]--;
				}
				act->w[act->wt]++;
			}
			else{//delete
				act->i[act->it]++;
				mut = M_DELETE;
				//simply increment the read head without doing anything else
			}
			act->r[act->rt]++;
			//act->i[act->it]++;
		}
		else{
			if(rno<subrate+indelrate){//INCREMENTAL MUTATION - we should never be doing this either
				doprint=1;
				cidx = sym_from_adj(*(act->r[act->rt]),blosum);
				rm = tab_idx(cidx,blosum);
				if(mass[rm]){
					*(act->w[act->wt])=cidx;
					act->w[act->wt]++;
					if(!(wm<0)){
						mass[wm]++;
					}
					mass[rm]--;
				}
				else{//
					//simply increment the read head without doing anything else
				}
				act->r[act->rt]++;//possible deletion here...
			}
			else{//NO MUTATION (but possible sub via comass effects)
				//cidx = sym_from_adj(*(act->r[act->rt]),blosum);
				rm = tab_idx(*(act->r[act->rt]),blosum);
				if(mass[rm]){
					*(act->w[act->wt])=*(act->r[act->rt]);
					act->w[act->wt]++;
					if(!(wm<0)){
						mass[wm]++;
					}
					mass[rm]--;
				}
				else{
					cidx = sym_from_adj(*(act->r[act->rt]),blosum);
					rm = tab_idx(cidx,blosum);
					if(mass[rm]){
						*(act->w[act->wt])=cidx;
						act->w[act->wt]++;
						if(!(wm<0)){
							mass[wm]++;
						}
						mass[rm]--;
					}
				}
				act->r[act->rt]++;
			}
		}
	}
	//update lengths
	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);
	act->i[act->it]++;

#ifdef VERBOSE
	if(mut)
	printf("Mutant event %d. new string is:\n%s\n\n",mut,act->wt?act->S:act->pass->S);
#endif
	act->biomass++;
	biomass++;
	return 0;
}


/*
int stringPM::comass_hcopy(s_ag *act){

	//s_ag *pass;
	//pass = act->pass;
	int cidx;
	float rno;
	int safe = 1;// this gets set to zero if any of the tests fail..
	int doprint = 0;
	e_mut mut = M_NONE;
	int doindel;


	//MUTATION RATES:
	//THESE ARE HARD-CODED FOR NOW - THEY SHOULD BE DERIVED FROM THE BLOSUM SOMEHOW...
	//const float indelrate = 0.0005,subrate=0.375;//0.0749/2
	//const float indelrate = 0.0000306125,subrate=0.01;//0.02
	//const float indelrate = 0.00006125,subrate=0.05;//0.02
	//const float indelrate = 0.000125,subrate=0.1;//0.02
	//const float indelrate = 0.000005,subrate=0.0375;//0.0749/2

	//const float indelrate = 0., 		subrate = 0.;

	if(!domut){
		indelrate = subrate =0;
	}

	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);

	int p;
	if( (p = h_pos(act,'w'))>=maxl){
		printf("Write head out of bounds: %d\n",p);
		//just to make sure no damage is done:
		if(act->wt)
			act->S[maxl]='\0';
		else
			act->pass->S[maxl]='\0';

		act->i[act->it]++;
		safe = 0;
		return -1;
	}
	if(h_pos(act,'r')>=maxl){
		printf("Read head out of bounds\n");
		act->i[act->it]++;
		safe = 0;
		return -2;
	}




	if(*(act->r[act->rt]) == 0){
		//possibly return a negative value and initiate a b
		safe = 0;
		//return -3;
	}

	if(safe){
		//we are going to use conservation of mass to drive mutation.
		int rm = tab_idx(*(act->r[act->rt]),blosum);
		int wm = -1;

		doindel=0;

		//see if we are overwriting or not:
		if(*(act->w[act->wt])){
			wm=tab_idx(*(act->w[act->wt]),blosum);
		}

		if(!mass[rm]){//if there's no mass left, pick a neighbour
			cidx = sym_from_adj(*(act->r[act->rt]),blosum);
			rm = tab_idx(cidx,blosum);
			if(mass[rm]==0)
				doindel = 1;
		}
		else{
			cidx = *(act->r[act->rt]);
		}

		if(doindel){
			rno=rand0to1();
			if(rno>0.5){//insert;

				//Then pick a random instruction:
				cidx = (float) rand0to1() * blosum->N;

				//insert the random instruction
				*(act->w[act->wt])=blosum->key[cidx];

				act->r[act->rt]++;
				act->w[act->wt]++;
			}
			else{//"delete" by skipping
				act->r[act->rt]++;
			}
		}
		else{
			*(act->w[act->wt]) = cidx;
			mass[rm]--;
			if(wm>-1)
				mass[wm]++;

			act->r[act->rt]++;
			act->w[act->wt]++;
		}
	}
	//update lengths
	act->len = strlen(act->S);
	act->pass->len = strlen(act->pass->S);
	act->i[act->it]++;

#ifdef VERBOSE
	if(mut)
	printf("Mutant event %d. new string is:\n%s\n\n",mut,act->wt?act->S:act->pass->S);
#endif
	act->biomass++;
	biomass++;
	return 0;
}
*/


int stringPM::comass_exec_step(s_ag *act, s_ag *pass){//pset *p,char *s1, swt *T){

	int finished=0;
	char *tmp;
	int dac=0;
	int safe_append=1;



	//if(act->idx==1222 || pass->idx==1222){
	//	print_exec(stdout,act,act->pass);
	//	fflush(stdout);
	//}

	//The following is unstable - explosions/catastrophe - possibly because of the nature of the decay??
	//Now moved energy to calling function
	//if(eqn_prop(energy)){
	//if(energy){

	switch(*(act->i[act->it])){//*iptr[it]){

	case '$'://h-search
		//act->ft = act->it;
		char *cs;
		if(act->ft)
			cs = act->S;
		else
			cs = act->pass->S;
		tmp = HSearch(act->i[act->it],cs,blosum,&(act->it),&(act->ft),maxl);
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
		if(comass_hcopy(act)<0){
			unbind_ag(act,'A',1,act->spp,pass->spp);
			unbind_ag(pass,'P',1,act->spp,pass->spp);
			finished = 1;
		}
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
			act->i[act->it]=IfLabel(act->i[act->it],act->r[act->rt],act->S,blosum,maxl);
			break;


	/************
	 *  CLEAVE  *
	 ************/
	case '%':
			if((dac = cleave(act))){
				//unbind_ag(act);
				//unbind_ag(pass);
				finished = 1;
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
			unbind_ag(act,'A',1,act->spp,pass->spp);
			unbind_ag(pass,'P',1,act->spp,pass->spp);

			finished = 1;
			break;

	default://Just increment the i-pointer
		act->i[act->it]++;
		break;
	}
#ifdef V_VERBOSE
	printf("Exec step - looks like:\n");
	print_exec(stdout,act,pass);
#endif


	if(safe_append){
		act->ect++;
		append_ag(&nexthead,act);
		append_ag(&nexthead,pass);
	}
	energy--;



	return 1;
	//return finished;
}



//This to be called AFTER agents and blosum have been loaded.
//Sets the mass for everything to a single value, read from config.
int stringPM::load_comass(char *fn, int verbose){
	//int massval = 2000;
	FILE *fp;
	int finderr;
	int i;
	int massval;

	mass = (int*) malloc(blosum->N * sizeof(int));

	if((fp=fopen(fn,"r"))!=NULL){
		finderr=read_param_int(fp,"MASS",&massval,verbose);
		fclose(fp);
		//Load the massvalue into the mass table;
		for(i=0;i<blosum->N;i++){
			mass[i]=massval;
		}
		s_ag *pag;
		for(pag=nowhead;pag!=NULL;pag=pag->next){
			update_mass(pag->S,strlen(pag->S),-1, 1);
		}
	}
	return 0;
}


int stringPM::set_mass(int *param){
	//int massval = 2000;
	FILE *fp;
	int finderr;
	int i;
	int massval;

	mass = (int*) malloc(blosum->N * sizeof(int));

	for(i=0;i<blosum->N;i++){
		mass[i]=param[i];
	}
	s_ag *pag;
	for(pag=nowhead;pag!=NULL;pag=pag->next){
		if(update_mass(pag->S,strlen(pag->S),-1, 1)){
			//update_lineage(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp, int mass){
			update_lineage(pag,'I',1,NULL,NULL,0);
		}
	}

	return 0;
}





int stringPM::update_mass(char *S, int len, int val, const int doconcat){
	int c;
	int updated=0;
	for(int i=0;i<len;i++){
		c = tab_idx(S[i],blosum);
		mass[c] += val;
		//Strategy: if there isn't the mass available for the string, delete the code
		if(doconcat && mass[c]<0){
			mass[c] -= val;
			for(int j = i+1; j< len; j++){
				S[j-1] = S[j];
			}
			S[--len]=0;
			i--; //decrement i so that we are checking the left-shifted character next
			updated=1;
		}
	}
	return updated;
}


int stringPM::comass_free_ag(s_ag *pag){

	if(pag->S != NULL){
		update_mass(pag->S,strlen(pag->S),1,0);//TODO: the last argument specifies whether to concatenate if no mass available, but since we are adding, it'll never be called
		//printf("destroying agent %d, code = %s\n",pag->idx,pag->S);
		free(pag->S);
	}

	free(pag);
	pag = NULL;

	return 0;
}


int stringPM::comass_testdecay(s_ag *pag){

#ifdef LONG_DECAY
	//VARIABLE DECAY RATE BASED ON LENGTH OF STRING (ECAL 2009)
	float len = strlen(pag->S);
	float prob = 1./pow(len,2);//4./3.);
#else
	//CONSTANT DECAY RATE to match ECAL (ALife 2010 and on)
	float prob = 1./pow(65,2);//4./3.);
#endif


	float rno = rand0to1();

#ifdef UNB_DECAY_ONLY
	if(rno<prob && pag->status == B_UNBOUND){
#else
	if(rno<prob && dodecay){
#endif
		//unbind_ag(pag);
		comass_free_ag(pag);
		return 1;
	}
	else
		return 0;
}

void stringPM::comass_make_next(){
	s_ag *pag,*bag;
	s_bind bb;
	int changed;

	//SUGGEST: write function to count what's around (saves rechecking every time)
	//countstates();

	//SUGGEST: IF WE ARE PUTTING ENERGY IN FROM OUTSIDE, AND NO ENERGY CAN BE PRODUCED
	//while(nowhead!=NULL && energy){
	//}
	//if(nowhead != NULL){
	//	pag = nowhead;
	//	append_ag(pag); //check that pag->next is also appended
	// NB: will have to have a different scheme for decay, since that doesn't require energy - sample from a binomial?

	while(nowhead!=NULL){

		pag = rand_ag(nowhead,-1);
		extract_ag(&nowhead,pag);
		changed = 0;

		//if(0){//extit>=10){
		//	if(pag->status==B_PASSIVE)
		//		//if(pag->exec->idx==1041){
		//			print_exec(stdout,pag->exec,pag);
		//		//}
		//	if(pag->status==B_ACTIVE){
		//		//if(pag->idx==1041)
		//			print_exec(stdout,pag,pag->pass);
		//	}
		//	fflush(stdout);
		//}

		//extract any partner:
		bag = NULL;
		bb = pag->status;
		switch(pag->status){
		case B_UNBOUND:
			break;
		case B_ACTIVE:
			//if(pag->idx == 214 && extit > 79)
			//	printf("Active 214\n");
			bag = pag->pass;
			extract_ag(&nowhead,bag);
			break;
		case B_PASSIVE:
			bag = pag->exec;
			extract_ag(&nowhead,bag);
			break;
		}

		if(pag->label=='P' && pag->status != B_UNBOUND  ){
			if(print_agent_idx(stdout,1,pag->idx))
				fflush(stdout);
		}

		int dc = comass_testdecay(pag);
		if(dc){//we must check what else needs to be destroyed...
			if(bag!=NULL){
				comass_free_ag(bag);
			}
		}
		else{
			if(energy>0){
				switch(pag->status){
				case B_UNBOUND:
					//seek binding partner, set binding states.
					changed = testbind(pag);
					//if(!changed)
					//	changed = testdecay(pag);

					break;
				case B_PASSIVE:

					//extract_ag(&nowhead,pag->exec);
					changed = comass_exec_step(pag->exec,pag);

					break;
				case B_ACTIVE:

					//extract_ag(&nowhead,pag->pass);
					changed = comass_exec_step(pag,pag->pass);
					break;
				default:
					printf("ERROR: agent with unknown state encountered!\n");
				}
			}
			if(!changed){
				append_ag(&nexthead,pag);
				if(bag!=NULL)
					append_ag(&nexthead,bag);

			}
		}
	}
}





/////////////////////////////////////////////////////////end of comass stuff




/////////////////////////////////////////////////////////start of energetic stuff


int stringPM::energetic_exec_step(s_ag *act, s_ag *pass){//pset *p,char *s1, swt *T){

	int finished=0;
	char *tmp;
	int dac=0;
	int safe_append=1;

	if(energy>0){
			switch(*(act->i[act->it])){//*iptr[it]){

				case '$'://h-search
					//act->ft = act->it;
					char *cs;
					if(act->ft)
						cs = act->S;
					else
						cs = act->pass->S;
					tmp = HSearch(act->i[act->it],cs,blosum,&(act->it),&(act->ft),maxl);
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
					if(hcopy(act)<0){
						unbind_ag(act,'A',1,act->spp,pass->spp);
						unbind_ag(pass,'P',1,act->spp,pass->spp);
						finished = 1;
					}
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
						act->i[act->it]=IfLabel(act->i[act->it],act->r[act->rt],act->S,blosum,maxl);
						break;


				/************
				 *  CLEAVE  *
				 ************/
				case '%':
						if((dac = cleave(act))){
							finished = 1;
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
						unbind_ag(act,'A',1,act->spp,pass->spp);
						unbind_ag(pass,'P',1,act->spp,pass->spp);

						finished = 1;
						break;

				default://Just increment the i-pointer
					act->i[act->it]++;
					break;
		}
	}
	else{//Just increment the i-pointer
		act->i[act->it]++;
	}

#ifdef V_VERBOSE
	printf("Exec step - looks like:\n");
	print_exec(stdout,act,pass);
#endif


	if(safe_append){
		act->ect++;
		append_ag(&nexthead,act);
		append_ag(&nexthead,pass);
	}
	if(energy>0)
		energy--;



	return 1;
	//return finished;
}



int stringPM::energetic_testbind(s_ag *pag){

	int found=0;
	int count=0;
	const int nrglim = 3000;
	float
		rno;
	s_ag *bag;
	align sw;
	float bprob =0.;

	bag=nowhead;

	count = nagents(nowhead,B_UNBOUND);

	//Calc propensity
	if(eqn_prop(count)){
		found = 1;
	}

	if(found){
		bag = rand_ag(nowhead,B_UNBOUND);
#ifndef BIND_ALL
		bprob = get_sw(pag,bag,&sw);
#else
		bprob =1.0;
		sw.match = 1;		// the number of matching characters.
		sw.score = 1; 	// the score of the match
		sw.prob = 1.0;		// the probability of the match - used for determining events based on the score/match
		sw.s1=0;			// start of the match in string 1
		sw.e1=1;			// end of the match in string 1
		sw.s2=0;			// start of the match in string 2
		sw.e2=1;			// end of the match in string 2
#endif
		rno = rand0to1();
		if(rno<bprob && ((1-bprob)*nrglim < energy   ) ){//Binding success!
			//figure out which is the executing string:
			set_exec(pag,bag,&sw);
			pag->nbind++;
			bag->nbind++;

			append_ag(&nexthead,pag);

			//uptime the second agent;
			extract_ag(&nowhead,bag);
			append_ag(&nexthead,bag);
			energy--;
		}
		else{
			found = 0;
		}
	}
	note_propensity(count,found);
	return found;
}


void stringPM::energetic_make_next(){
	s_ag *pag,*bag;
	s_bind bb;
	int changed;
	while(nowhead!=NULL){

		pag = rand_ag(nowhead,-1);
		extract_ag(&nowhead,pag);
		changed = 0;

		//extract any partner:
		bag = NULL;
		bb = pag->status;
		switch(pag->status){
		case B_UNBOUND:
			break;
		case B_ACTIVE:
			//if(pag->idx == 214 && extit > 79)
			//	printf("Active 214\n");
			bag = pag->pass;
			extract_ag(&nowhead,bag);
			break;
		case B_PASSIVE:
			bag = pag->exec;
			extract_ag(&nowhead,bag);
			break;
		}

		if(pag->label=='P' && pag->status != B_UNBOUND  ){
			if(print_agent_idx(stdout,1,pag->idx))
				fflush(stdout);
		}

		int dc = testdecay(pag);
		if(dc){//we must check what else needs to be destroyed...
			if(bag!=NULL){
				free_ag(bag);
			}
		}
		else{


			switch(pag->status){
			case B_UNBOUND:
				if(energy>0){
					//seek binding partner, set binding states.
					changed = energetic_testbind(pag);
					//if(!changed)
					//	changed = testdecay(pag);
				}

				break;
			case B_PASSIVE:

				//extract_ag(&nowhead,pag->exec);
				changed = energetic_exec_step(pag->exec,pag);

				break;
			case B_ACTIVE:

				//extract_ag(&nowhead,pag->pass);
				changed = energetic_exec_step(pag,pag->pass);
				break;
			default:
				printf("ERROR: agent with unknown state encountered!\n");
			}


			if(!changed){
				append_ag(&nexthead,pag);
				if(bag!=NULL)
					append_ag(&nexthead,bag);

			}
		}
	}
}
/////////////////////////////////////////////////////////end of energetic stuff













void stringPM::update(){
	//s_ag *pag;
	nowhead = nexthead;
	nexthead = NULL;
	//pag = nowhead;
	//while(pag!=NULL){
	//	//move_ag(pag);
	//	pag = pag->next;
	//}
	//update_aac();
}



/*
//INFLUX STUFF
void stringPM::influx(int t){

	s_ix *pix;
	s_ag *pag;
	float rno;

	int irx=0;

	pix = ifxhead;

	while(pix!=NULL){

		if(t>=pix->start && t<=pix->stop){
			for(int i=0;i<pix->n;i++){
				rno = rand0to1();
				irt[irx]++;
				if(rno<pix->prob){
					irf[irx]++;
					//pag = make_ag(pix->label,1);

					pag = make_ag(pix->label,0);
					//SUGGEST: load and write the string for the influx'd dudes

					append_ag(&nowhead,pag);
				}
			}
		}
		pix=pix->next;
	}
	update_aac();
}


void stringPM::influx_special(int t){

	//These are the variables we tinker with!
	int step = 3.* log(0.5)/log(1-0.00004);//100000;
	float gprob = 0.02;
	float lprob = 0.0225;
	s_ag *pag;
	float prob;
	int lab;


	if((t/step)%2){
		lab='I'; //may as well be something we can detect!

		//starvation:
		//prob=0;

		//lactose
		prob=lprob;
		//lab='I';
	}
	else{
		//glucose
		prob=gprob;
		lab='A';
	}

	float rno = rand0to1();
	//int ct = aac_count(lab);
	if(rno<prob){//&& ct<500){
		//pag = make_ag(lab,1);


		pag = make_ag(lab,0);

		append_ag(&nowhead,pag);
	}
}



void stringPM::replenish_operons(){
	int i;
	s_ag *pag;
	pag=nowhead;

	int Jc,Oc,qc,dc;
	Jc=Oc=qc=dc=0;
	int count,a,j;
	update_aac();
	for(i=0;i<ntt;i++){
		if(aro[i]>-1){
			//count operons
			count=0;
			for(j=0;j<ntt;j++)
				if(com[i][j])
					count+=aac[j];
			if(aro[i]>count){
				for(a=0;a<aro[i]-count;a++){
					pag = make_ag(aat[i],1);
					append_ag(&nowhead,pag);
				}
				//aac[i]+=aro[i]-count;
			}
		}
	}

	update_aac();
}

void stringPM::divide(){

	float rand;
	s_ag *pag,*pag2;
	pag=nowhead;
	while(pag!=NULL){
		rand = rand0to1();
		if(rand<0.5){//Delete it.
			pag2=pag;
			pag=pag->next;
			extract_ag(&nowhead,pag2);
			free_ag(pag2);
		}
		else{
			//rand_in_rad(cellrad,&(pag->x),&(pag->y));
			pag=pag->next;
		}
	}
	printf("After div: ");
	update_aac();
	print_agents_count(stdout);
	replenish_operons();
	energy = energy/2.;
}
*/

void stringPM::free_swt(swt *pSWT, int verbose){


	//print for debug.
	int i,j;
	if(verbose){
		printf("   ");
		for(i=0;i<pSWT->N;i++)
			printf("%c   ",pSWT->key[i]);
		printf("\n");
		for(i=0;i<pSWT->N;i++){
			printf("%c  ",pSWT->key[i]);
			for(j=0;j<pSWT->N;j++){
				printf("%0.2f ",pSWT->T[i][j]);
			}
			printf("\n");
		}
		//Print the indel value
		printf("-   ");
		for(i=0;i<pSWT->N;i++){
			printf("%02d ",(int) pSWT->T[pSWT->N][i]);
		}
		fflush(stdout);
	}

	float *pT;
	int *pI;

	//Free elements of pSWT
	if(pSWT->T != NULL){
		free(pSWT->key);
		pSWT->key= NULL;
		int i;
		for(i=0;i<pSWT->N+1;i++){
			pT = pSWT->T[i];
			free(pT);
		}
		free(pSWT->T);
		pSWT->T=NULL;

		for(i=0;i<pSWT->N;i++){
			pI = pSWT->adj[i];
			free(pI);
		}
		free(pSWT->adj);
		pSWT->adj = NULL;
	}


}

/** Have to specifically ask for verbose output now... */
void stringPM::clearout(int verbose){

	s_ag	*agp,*agp2;

	if(verbose){
		printf("Starting spatial clearout..");fflush(stdout);
	}


	free_swt(blosum,verbose);

	agp = nowhead;
	while(agp!=NULL){
		agp2=agp->next;
		free_ag(agp);
		agp=agp2;
	}

	agp=nexthead;
	while(agp!=NULL){
		agp2=agp->next;
		free_ag(agp);
		agp=agp2;
	}

	/*
	//TODO: put this in SMspp->clearout();
	l_spp *sp,*sp2;
	sp = spl->species;
	while(sp!=NULL){
		sp2=sp->next;
		free_lspp(sp);
		sp=sp2;
	}
	spl->species = NULL;
	spp_count=1;
	*/

	//TODO: we don't need to do this, it happens anyway
	//agents_base::clearout();

	free_swlist(&swlist);

	preset();

}


void stringPM::sanity_check(){
	s_ag *pag;
	pag = nowhead;
	while(pag!=NULL){

		if(pag->S[maxl]!=0){
			printf("Problem with string %d, status ",pag->idx);
			printf(" out of range: %s", pag->S);

		}

		pag=pag->next;

	}
}

//SPECIES ANALYSIS FUNCTIONS

/*
l_spp * stringPM::make_lspp(s_ag *a,fixthis){

	l_spp *sp;

	//printf("Spatial make_ag called\n");fflush(stdout);

	if((sp = (l_spp *) mymalloc(1,sizeof(l_spp)))!=NULL){

	    sp->sp = splist->getspp(a);//, a->pp);//->pa, a->pp->pp );

	    // THE FOLLOWING ARE NOW IN THE SMspp CLASS
		//sp->S=(char *) malloc(maxl0+1*sizeof(char));
		//memset(sp->S,0,maxl0*sizeof(char));
		//l = strlen(a->S);
		//strncpy(sp->S,a->S,l);
		//sp->paspp = a->paspp;
		//sp->ppspp = a->ppspp;
		//sp->spp=spp_count++;
		//sp->sptype=0;

		sp->count=0;
		sp->first=-1;
		sp->next=NULL;
		sp->tspp=extit;
		return sp;
	}
	else{
		printf("Error allocating in make_spp\n");
		return NULL;
	}
}
*/

void stringPM::free_lspp(l_spp *sp){
	//free(sp->S);
	free(sp);
}

/*
int stringPM::append_lspp(l_spp *sp){

	l_spp *p;
	//printf("appending: list = %p, ag = %p\n",*list,ag);
	if(species==NULL){
		species=sp;
		//printf("appended: list = %p, ag = %p\n",*list,ag);
	}
	else{
		p = species;
		while(p->next != NULL){
			p = p->next;
		}
		p->next = sp;
	}
	return 0;
}
*/


//This is called from CLEAVE, otherwise we can't tell if its in the middle of being constructed....
int stringPM::update_lineage(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp, int mass){

	l_spp *sp;
	sp = spl->species;
	int found=0;
	int novel=0;

	//Step 1: sort out if it is a new species:
	/* It'd be nice to do this, but we don't know what the species IS yet!
	while(sp!=NULL&&!found){

		if(p->spp->spp == sp->spp){//We've found a matching species index.
			if(add){

				if(!strcmp(sp->S,p->S)){//we've found a match on the string
					//sp->count++;
					break;
				}
				else{//Need a second check to see if the string matches anything else
					for(sp2=spl->species;sp2!=NULL;sp2=sp2->next){
						if(!strcmp(sp2->S,p->S)){//we've found a match on the string
							sp->count++;
							novel=1;
						}
						break;
					}
				}
			}
		}
	}
	*/

	while(sp!=NULL&&!found){
		if(!strcmp(sp->S,p->S)){//we've found a match on the string
			found=sp->spp;
			if(add){
				//Need to check whether this is a dissociating partner in a reaction that has changed:
				if(p->spp!=NULL){//We haven't decided what the spp is yet
					if(p->spp->spp != sp->spp){
						sp->count++;
						novel=1;
					}
					//novel=0 IF dissociating and no new spp are produced.
				}
				else{
					sp->count++;
					novel=1;
				}
			}
			break;
		}
		sp = sp->next;
	}

	if(add){//Only do this if we are adding to the list (not just checking reaction-space
		if(!found){
			sp = spl->make_spp(p,extit,maxl0);
			if(signal!=NULL){
				sp->sig_sc = sigalign(sp->S);
			}
			else
				sp->sig_sc=0;

			sp->sptype = sptype;
			sp->biomass += mass;
			spl->prepend_spp(sp); //append_lspp(sp);
			novel=1;
		}

		//Now sort out the parentage:
		p->spp = sp;//->spp;
		if(novel)
			p->pp = spl->get_parents(p->spp, paspp, ppspp);

	}

	return found;
}


l_spp * stringPM::get_spp(int n){
	l_spp * sp;
	sp = spl->species;
	while(sp!=NULL){
		if(sp->spp == n)
			return sp;
		sp = sp->next;
	}

	printf("Species %d \n",n);
	return NULL;

}


int stringPM::get_ecosystem(){
	l_spp *sp;
	s_ag *ag;
	int maxsp;
	int maxct;

	//reset the flags
	sp=spl->species;
	while(sp!=NULL){
		sp->pf=sp->anc=0;
		sp = sp->next;
	}

	//Find the current ecosystem
	ag=nowhead;
	while(ag!=NULL){
		sp=spl->species;
		while(sp!=NULL){
			if(sp->spp == ag->spp->spp){
				sp->pf++;
				break;
			}
			sp=sp->next;
		}
		ag=ag->next;
	}

	maxsp= spl->species->spp;
	maxct=spl->species->pf;
	sp=spl->species;
	while(sp!=NULL){
		if(sp->pf > maxct){
			maxsp = sp->spp;
			maxct = sp->pf;
		}
		sp=sp->next;
	}

	return maxsp;
}



void stringPM::print_ancestry_dot(FILE *fp, int time,int step){

	l_spp *sp;
	s_ag *ag;

	const int minno = 1; //There must be more than this many individuals present for the ancestry to be plotted.....

	if(fp==NULL)
		fp=fopen("ancestry.dot","w");

	//reset the flags
	sp=spl->species;
	while(sp!=NULL){
		sp->pf=sp->anc=0;
		sp = sp->next;
	}

	//Find the current ecosystem
	ag=nowhead;
	while(ag!=NULL){
		sp=spl->species;
		while(sp!=NULL){
			if(sp->spp == ag->spp->spp){
				sp->pf++;
				break;
			}
			sp=sp->next;
		}
		ag=ag->next;
	}

	//Now trace it's origins...
	int finished = 0;
	l_spp *sp2;
	while(!finished){
		finished=1;
		sp=spl->species;
		while(sp!=NULL){
			if(sp->pf > minno || sp->anc){

				if(sp->pp->pp != NULL){
					sp2 = get_spp(sp->pp->pp->spp);
					if(sp2==NULL){
						printf("NULL sp2\n");fflush(stdout);
					}
					//if(sp2->sptype){//This means its not a "start" type...
						//if(!sp2->pf){
					if(!sp2->anc){
						sp2->anc=1;
						finished=0;
					}
				}
					//}
				if(sp->pp->pp != NULL){
					sp2 = get_spp(sp->pp->pa->spp);
					//if(sp2->sptype){//This means its not a "start" type...
					if(!sp2->anc){
						sp2->anc=1;
						finished=0;
					}
					//}
				}

			}
			sp=sp->next;
		}
	}

	//Now we can print it!!!


	fprintf(fp,"digraph simple_hierarchy {\noverlap = scale\n");

	//Do the timeline first:	//TIME LINE GRAPH:
	fprintf(fp,"/* TIME LINE */\n{node [shape=plaintext, fontsize=16];\n");
	fprintf(fp,"0");
	sp = spl->species;
	while(sp!=NULL){
		if(sp->tspp){
			if(sp->pf > minno || sp->anc){
				fprintf(fp," -> %d",sp->tspp);
			}
		}
		sp=sp->next;
	}
	//Have a line for NOW:
	fprintf(fp," -> %d",time);
	fprintf(fp,";\n}\n");


	fprintf(fp,"\n\n/* ANCESTORS: */");
	sp=spl->species;
	while(sp!=NULL){

		//if(sp->pf<0 || sp->pf > minno){//It's an ancestor
		if(sp->anc || sp->pf > minno){//Its an ancestor, NOT a minno
			if(!sp->pf)
				fprintf(fp,"\t\t\t {rank=same; %d;  anc%03d [style=filled,color=\"grey\",label=\"%03d\"]}\n",sp->tspp,sp->spp,sp->spp);
			else
				fprintf(fp,"\t\t\t {rank=same; %d;  anc%03d [label=\"%03d\"]}\n",sp->tspp,sp->spp,sp->spp);

			fprintf(fp,"\t\t\t cl%03d [shape=box,width=0.15,height=0.15,style=filled,color=",sp->spp);
			switch(sp->sptype){
			case 'C':
				fprintf(fp,"green");
				break;
			case 'A':
				fprintf(fp,"red");
				break;
			case 'P':
				fprintf(fp,"black");
				break;
			default:
				fprintf(fp,"grey");
			}
			fprintf(fp,",label=\"\"]\n");
			if(sp->pp->pa!=NULL){
				fprintf(fp,"\t\t\t anc%03d -> cl%03d\n",sp->pp->pa->spp,sp->spp);
				sp2=get_spp(sp->pp->pa->spp);
				if(sp2->pf>0 && sp2->pf <=minno){
					fprintf(fp,"\t\t\t {rank=same; %d;  anc%03d [label=\"%03d\"]}\n",sp2->tspp,sp2->spp,sp2->spp);
				}
			}
			if(sp->pp->pp!=NULL){
				fprintf(fp,"\t\t\t anc%03d -> cl%03d [color=grey]\n",sp->pp->pp->spp,sp->spp);
				sp2=get_spp(sp->pp->pp->spp);
				if(sp2->pf>0 && sp2->pf <=minno){
					fprintf(fp,"\t\t\t {rank=same; %d;  anc%03d [label=\"%03d\"]}\n",sp2->tspp,sp2->spp,sp2->spp);
				}
			}
			if(sp->pp->pa!=NULL || sp->pp->pp!=NULL)
				fprintf(fp,"\t\t\t cl%03d -> anc%03d\n",sp->spp,sp->spp);
		}



		sp=sp->next;
	}


	fprintf(fp,"\n\n/* CURRENT ECOSYSTEM: */");
	sp=spl->species;
	int minnoct=0;
	while(sp!=NULL){
		if(sp->pf>minno){
			fprintf(fp,"\t\t\t {rank=same; %d;  spp%03d [style=filled,color=black,fontcolor=white,label=\"%03d\"]}\n",time,sp->spp,sp->spp);
			fprintf(fp,"\t\t\t anc%03d -> spp%03d [color=red,headlabel=\"%d\"]\n\n",sp->spp,sp->spp,sp->pf);
		}
		//if(sp->pf>0 && sp->pf <= minno){
		//	minnoct += sp->pf;
		if(sp->pf && sp->pf <=minno){
			minnoct+= sp->pf;
		}
		sp=sp->next;
	}
	if(minnoct){
		fprintf(fp,"\t\t\t {rank=same; %d;  minnos [style=filled,color=black,fontcolor=red,label=\"rare\"]}\n",time);
		fprintf(fp,"ancminno [style=invis]");
		fprintf(fp,"\t\t\t ancminno -> minnos [color=red,headlabel=\"%d\"]\n\n",minnoct);
	}


	fprintf(fp,";\n}\n");
	fflush(fp);
	fclose(fp);

}


void stringPM::print_spp_strings(FILE *fp){


	if(fp==NULL)
		fp=fopen("species.txt","w");

	l_spp *sp;


	fprintf(fp,"SppNo\tTIME\tNumber\tP1\tP2\tbtype\tlen\tSTRING\n");
	for(sp=spl->species;sp!=NULL;sp=sp->next){
		int len = strlen(sp->S);
		//TODO: record all parent pairs.
		s_parent *ppp;
		for(ppp = sp->pp; ppp!= NULL; ppp=ppp->next){
			if(sp->pp->pa!=NULL || sp->pp->pp!=NULL){
				fprintf(fp,"%05d\t%06d\t%04d\t%d\t%d\t%d\t%d\t%s\t%p\n",sp->spp, sp->tspp, sp->pf, sp->pp->pp->spp, sp->pp->pa->spp, sp->sptype, len,  sp->S,  &(sp->S));
			}
		}
	}

	l_spp *ss;
	s_parent *pp;

	//For debugging - print the master list of species
	fprintf(fp,"\n\nMASTER SPECIES LIST");
	fprintf(fp,"SppNo\ttP1\tP2\tbtype\tlen\tSTRING\n");
	for(ss=spl->species ;ss!=NULL;ss=ss->next){
		int len = strlen(ss->S);
		for(pp=ss->pp; pp!=NULL; pp=pp->next){
			if(ss->pp->pa != NULL)
				fprintf(fp,"%05d\t%d\t%d\t%d\t%d\t%s\t%p\n",ss->spp, ss->pp->pp->spp, ss->pp->pa->spp, ss->sptype, len,  ss->S,  &(ss->S));
			else
				fprintf(fp,"%05d\t_\t_\t%d\t%d\t%s\t%p\n",ss->spp, ss->sptype, len,  ss->S,  &(ss->S));
		}
	}

	fflush(fp);
	fclose(fp);
}




int stringPM::Network_cleave(s_ag *act){

	int dac = 0,cpy;
	s_ag *c,*pass,*csite;

	pass = act->pass;

	//pick the mol containing the cleave site:
	csite = act->ft?act:pass;

	if(act->f[act->ft]-csite->S < csite->len){

		//1: MAKE THE NEW MOLECULE FROM THE CLEAVE POINT

		//Can't really say what the label is easily - for ECAL, it's always pass
		c = make_ag(pass->label,1);


		//if(c->idx>2375){
		//		fflush(stdout);
		//}

		//Copy the cleaved string to the agent

		char *cs;

		c->S =(char *) malloc(maxl0*sizeof(char));
		memset(c->S,0,maxl0*sizeof(char));

		cs = csite->S;
		cpy = strlen(cs);
		cpy -= act->f[act->ft]-cs;
		if(!cpy){
			printf("WARNING - zero length copy being created\n");
		}

		strncpy(c->S,act->f[act->ft],cpy);
		c->len = strlen(c->S);
#ifdef VERBOSE
		printf("String %d created:\n%s\n",c->idx,c->S);
#endif

		//Filling in the birth certificate is now done in update_lineage
		//c->pp = splist->get_parents(c->spp,act->spp,pass->spp);
		//c->paspp = act->spp;
		//c->ppspp = pass->spp;

		//append the agent to nexthead
		update_lineage(c,'C',0,act->spp,pass->spp,act->biomass);

		free_ag(c);
		//append_ag(&nexthead,c);

		//2: HEAL THE PARENT

		memset(act->f[act->ft],0,cpy*sizeof(char));

		csite->len = strlen(csite->S);

		if((dac = check_ptrs(act))){
			switch(dac){
			case 1://Destroy active - only append passive
				unbind_ag(pass,'P',0,act->spp, pass-> spp);
				//append_ag(&nexthead,pass);
				//free_ag(act);
				break;
			case 2://Destroy passive - only append active
				unbind_ag(act,'A',0,act->spp, pass-> spp);
				//append_ag(&nexthead,act);
				//free_ag(pass);
				break;
			case 3://Destroy both
				printf("This should never happen\n");
				unbind_ag(act,'A',0,act->spp, pass-> spp);
				unbind_ag(pass,'P',0,act->spp, pass-> spp);
				//free_ag(act);
				//free_ag(pass);
				break;
			default://This can't be rightcan it?
				if(act->ft == act->it){
					act->i[act->it]--;
				}
				break;
			}
		}
	}
	else{
		printf("Problem with cleave!\n");
		//getchar();
	}
	if(!dac){
		act->i[act->it]++;
		//dac=-1;
	}
	return dac;
}




int stringPM::Network_exec_step(s_ag *act, s_ag *pass){

	int finished=0;
	char *tmp; //Holds the modifier, where needed
	//s_ag *c;
	//int cpy;
	int dac=0;
	int safe_append=1;

	switch(*(act->i[act->it])){

	case '$'://h-search/n/staff/sjh/current/ewASP
		//act->ft = act->it;
				char *cs;
				if(act->ft)
					cs = act->S;
				else
					cs = act->pass->S;
				act->f[act->ft] = HSearch(act->i[act->it],cs,blosum,&(act->it),&(act->ft),maxl);
				act->i[act->it]++;
				break;

				/*
				//put the active flow head on the string where the active ip is:
				act->ft = act->it;
				char *cs;
				if(act->it)
					cs = act->S;
				else
					cs = act->pass->S;
				act->f[act->ft] = HSearch(act->i[act->it],cs,blosum);
				act->i[act->it]++;
				break;
				*/

	case '>'://p-move: ?iptr? to fptr
		//make sure we toggle ->it to ->ft!!!
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
				//iptr[it]++;
				break;
			}
			//iptr[it]++;
			break;


	/************
	 *   HCOPY  *
	 ************/
	case '='://h-copy
		if(hcopy(act)<0){
			unbind_ag(act,'A',0,act->spp, pass-> spp);
			unbind_ag(pass,'P',0,act->spp, pass-> spp);
			finished = 1;
		}
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
			act->i[act->it]=IfLabel(act->i[act->it],act->r[act->rt],act->S,blosum,maxl);
			break;


	/************
	 *  CLEAVE  *
	 ************/
	case '%':
			//cleave(p->f[p->ft],);
			//Create new agent with state unbound

			if((dac = Network_cleave(act))){
				//unbind_ag(act);
				//unbind_ag(pass);
				finished = 1;
				safe_append=0;	//extract_ag(&nowhead,p);
			}
			break;

			/*

			*/
	case 0:
	case '}'://ex-end - finish execution

//#ifdef V_VERBOSE
	//		printf("Unbinding...\n");
//#endif
			unbind_ag(act,'A',0,act->spp, pass-> spp);
			unbind_ag(pass,'P',0,act->spp, pass-> spp);

			finished = 1;
			break;

	default://Just increment the i-pointer
		act->i[act->it]++;
		break;
	}
//#ifdef V_VERBOSE
//	printf("Exec step - looks like:\n");
//	print_exec(stdout,act,pass);
//#endif

	return 1;
	//return finished;
}

void stringPM::get_spp_network(char *fn){

	float bprob;
	char *comp;
	align sw;
	l_spp *sp1,*sp2;
	s_ag *pag,*bag,*act,*pas,tmp;
	tmp.S=(char *) malloc(maxl0*sizeof(char));

	//Check for errors on spp_network.
	//s_ag *nerr;
	//nerr = nowhead;
	//while(nerr->next !=NULL){
	//	nerr = nerr->next;
	//}


	domut = 0;
	float oldindel = indelrate;
	float oldsub = subrate;

	FILE *fp;
	fp=fopen(fn,"w");
	int rno=0;
	int newmolc=0;
	int newmold=0;

	pag=make_ag('A',0);
	pag->S =(char *) malloc(maxl0*sizeof(char));
	bag=make_ag('B',0);
	bag->S =(char *) malloc(maxl0*sizeof(char));

	FILE *fpd;
	char fnd[128];
	sprintf(fnd,"%s.dot",fn);
	fpd = fopen(fnd,"w");
	fprintf(fpd,"digraph simple_hierarchy {\noverlap = scale\n");


	sp1 = spl->species;
	while (sp1!=NULL){
		sp2=spl->species;
		comp = string_comp(sp1->S);
		while (sp2!=NULL){
			if(sp1->pf && sp2->pf){
				//bprob = SmithWaterman(comp,sp2->sp->S,&sw,blosum,0);
				bprob = SmithWatermanV2(comp,sp2->S,&sw,blosum,0);
				if(sw.match){
					bprob = get_bprob(&sw);
					//printf("Match found between spp %d and %d\n",sp1->spp,sp2->spp);
					fprintf(fp,"\nBind details: sp%d: %d .. %d; sp%d: %d .. %d\n",sp1->spp,sw.s1,sw.e1,sp2->spp,sw.s2,sw.e2);
					fprintf(fp,"sp%d + sp%d -%f> r%d\n",sp1->spp,sp2->spp,bprob,++rno);

					fprintf(fpd,"\t\t\t r%d [shape=box,width=0.15,height=0.15,style=filled,label=\"\",color=\"green\"]\n",rno);
					fprintf(fpd,"\t\t\t sp%d -> r%d\n",  sp1->spp,rno);
					fprintf(fpd,"\t\t\t sp%d -> r%d\n",  sp2->spp,rno);


					//ok, now we have to get reaction data!
					//unbind_ag(&pag);
					//unbing_ag(&bag);
					memset(pag->S,0,maxl0*sizeof(char));
					memset(bag->S,0,maxl0*sizeof(char));
					strncpy(pag->S,sp1->S,maxl);
					strncpy(bag->S,sp2->S,maxl);
					pag->spp=sp1;//->spp;
					bag->spp=sp2;//->spp;
					pag->status=bag->status=B_UNBOUND;


					//figure out which is the executing string:
					set_exec(pag,bag,&sw);
					if(pag->status==B_ACTIVE){
						act=pag;
						pas=bag;
					}
					else{
						act=bag;
						pas=pag;
					}

					int nsteps=0;
					while(act->status==B_ACTIVE && nsteps < 100000){
						int found;
						//Commenting out to debug:

						if( *(act->i[act->it])=='%'){//A cleave
							memset(tmp.S,0,maxl0*sizeof(char));
							//strcpy(tmp.S,act->i[act->it]);
							strcpy(tmp.S,act->f[act->ft]);
							found = update_lineage(&tmp,'X',0,act->spp, pas-> spp, 0);
							fprintf(fp,"r%d -%f>  ",rno,1./(float) nsteps);

							if(!found){
								//we have to add a new molecule to the network...
								fprintf(fp," r%d + nc%d\n",++rno, ++newmolc);

								fprintf(fpd,"\t\t\t r%d [shape=box,width=0.15,height=0.15,style=filled,label=\"\",color=\"red\"]\n",rno);
								fprintf(fpd,"\t\t\t r%d -> r%d \n", rno-1, rno);
								fprintf(fpd,"\t\t\t r%d -> nc%d \n", rno-1, newmolc);
							}
							else{
								fprintf(fp," r%d + sp%d\n",++rno, found);

								fprintf(fpd,"\t\t\t r%d [shape=box,width=0.15,height=0.15,style=filled,label=\"\",color=\"red\"]\n",rno);
								fprintf(fpd,"\t\t\t r%d -> r%d \n", rno-1, rno);
								fprintf(fpd,"\t\t\t r%d -> sp%d \n", rno-1, found);
							}
							fflush(fp);
							nsteps = 0;
						}

						Network_exec_step(act,pas);
						nsteps++;
					}
					//*/



					//Report how long it took, and check if anything has changed...
					int out1,out2;
					out1 = update_lineage(act,'X',0,act->spp, pas-> spp, 0 );
					out2 = update_lineage(pas,'X',0,act->spp, pas-> spp, 0);
					fprintf(fp,"r%d -%f> ",rno,1./(float) nsteps);
					if(!out1){
						fprintf(fp,"nd%d + ",++newmold);
						fprintf(fpd,"\t\t\t r%d -> nd%d \n", rno, newmold);
					}
					else{
						fprintf(fp,"sp%d + ",sp1->spp);
						fprintf(fpd,"\t\t\t r%d -> sp%d \n", rno, sp1->spp);
					}

					if(!out2){
						fprintf(fp,"nd%d\n\n",++newmold);
						fprintf(fpd,"\t\t\t r%d -> nd%d\n\n", rno, newmold);
					}
					else{
						fprintf(fp,"sp%d\n\n",sp2->spp);
						fprintf(fpd,"\t\t\t r%d -> sp%d\n\n", rno, sp2->spp);
					}


				}
			}
			sp2 = sp2->next;
		}
		free(comp);
		sp1 = sp1->next;
	}
	fflush(fp);
	fclose(fp);

	fprintf(fpd,"\n}\n");
	fflush(fpd);
	fclose(fpd);

	//Free the strings
	free(tmp.S);
	free_ag(pag);
	free_ag(bag);

	//if(nerr->next != NULL)
	//	printf("Test agent assigned to real set!\n");

	domut = 1;
	indelrate = oldindel;
	subrate = oldsub;
}


void stringPM::set_epochs(){

	lastepoch=get_ecosystem();
	thisepoch=lastepoch;
	nepochs=1;
}


int stringPM::share_agents(s_ag **hp){

	s_ag *pa,**head,*tmp,*aa,*tmphead;
	int ntot=0,herect=0,therect=0,*pct,*nact;
	float rno;

	tmphead = *hp;
	pa = tmphead;
	*hp = NULL;

	while(pa !=NULL){
		ntot++;

		//Work out where it's going:
		if((rno = rand0to1()) < 0.5){
			head = hp;
			therect++;
			nact = &therect;
		}
		else{
			head = &nowhead;
			herect++;
			nact = &herect;
		}

		//Move it and any partners.
		tmp=pa->next;
		extract_ag(&tmphead,pa);
		append_ag(head,pa);

		switch(pa->status){
		case B_ACTIVE:
			aa = pa->pass;
			break;
		case B_PASSIVE:
			aa = pa->exec;
			break;
		case B_UNBOUND:
			aa = NULL;
			break;
		}

		if(aa){
			*pct++;
			ntot++;
			if(tmp==aa)
				tmp=tmp->next;
			extract_ag(&tmphead,aa);
			append_ag(head,aa);
		}

		pa=tmp;
	}

	printf("Total %d. Moved %d. Left %d\n",ntot,herect,therect);
	return 1;
}



//A[c]->copy_agents(A[c2]->nowhead);
//int stringPM::copy_agents(s_ag **head){

//	s_ag *aa,*pa,*tmp;
//	float rno;int ntot=0, ncop=0;

//	pa=*head;
//	while(pa!=NULL){

		/*
		aa = make_ag(pa->label,-324);
		aa->S =(char *) malloc(maxl0*sizeof(char));
		memset(aa->S,0,maxl0*sizeof(char));
		strncpy(aa->S,pa->S,strlen(pa->S));
		aa->len = strlen(aa->S);
		append_ag(&nowhead,aa);
		*/

	/*
		//It is actually EASIER to copy 1/2 the cell contents, since then we can just move stuff rather than
		//restructure all the reactive environments. More realistic too!
		ntot++;
		if((rno = rand0to1()) < 0.5){
			ncop++;
			tmp=pa->next;
			switch(pa->status){
			case B_ACTIVE:

				aa = pa->pass;
				//make up for moving two across
				//tmp=tmp->next;
				//make sure tmp isn't about to be extracted
				if(tmp==aa)
					tmp=tmp->next;
				extract_ag(head,pa);
				append_ag(&nowhead,pa);
				extract_ag(head,aa);
				append_ag(&nowhead,aa);
				ncop++;

				break;
			case B_PASSIVE:

				aa = pa->exec;
				//make up for moving two across
				//tmp=tmp->next;
				//make sure tmp isn't about to be extracted
				if(tmp==aa)
					tmp=tmp->next;
				extract_ag(head,pa);
				append_ag(&nowhead,pa);
				extract_ag(head,aa);
				append_ag(&nowhead,aa);
				ncop++;

				break;
			case B_UNBOUND:
				extract_ag(head,pa);
				append_ag(&nowhead,pa);
				break;
			}
			pa=tmp;
		}
		else{

			pa = pa->next;
		}

	}
	printf("Copied %d agents out of %d\n",ncop,ntot);
	return 1;
}
*/


void print_agent_cfg(FILE *fp, s_ag *pa){

	if(pa->status == B_ACTIVE){
		fprintf(fp,"############### ,%d,%s\n",pa->spp->spp,pa->spp->S);
		fprintf(fp,"REACTION_ACTIVE ,%d,%s,irwf:",pa->spp->spp,pa->S);
		//IRWF
		fprintf(fp,",%d,%d,%d",pa->it , (int) (pa->i[0]-&(pa->pass->S[0])) , (int) (pa->i[1]-&(pa->S[1])) );
		fprintf(fp,",%d,%d,%d",pa->rt , (int) (pa->r[0]-&(pa->pass->S[0])) , (int) (pa->r[1]-&(pa->S[1])) );
		fprintf(fp,",%d,%d,%d",pa->wt , (int) (pa->w[0]-&(pa->pass->S[0])) , (int) (pa->w[1]-&(pa->S[1])) );
		fprintf(fp,",%d,%d,%d",pa->ft , (int) (pa->f[0]-&(pa->pass->S[0])) , (int) (pa->f[1]-&(pa->S[1])) );
		fprintf(fp,"\n");
	}
	else{
		fprintf(fp,"REACTION_PASSIVE,%d,%s\n",pa->spp->spp,pa->S);
		fprintf(fp,"############### ,%d,%s\n\n",pa->spp->spp,pa->spp->S);
	}

}


int stringPM::print_conf(FILE *fp){


	if(fp==NULL)
		fp=stdout;



/* TODO: record RANDSEED and NCONTAINERS values */


/*
	%%%CELL PARAMETERS
	CELLRAD		2500
	AGRAD		10
	ENERGY		0
	NSTEPS		1200000000
	USING	/n/staff/sjh/current/ewSTRING/SMconfigs/instr_set1.mis
*/
	fprintf(fp,"%%%%%%AUTOMATICALLY GENERATED Stringmol CONFIG FILE\n");
	fprintf(fp,"%%%%%%Generated at time: %d\n\n",extit);
	fprintf(fp,"%%%%%%CELL PARAMETERS\n");
	fprintf(fp,"CELLRAD		%d\n",(int) cellrad);
	fprintf(fp,"AGRAD		%d\n",(int) agrad);
	fprintf(fp,"ENERGY		%d\n",(int) energy);
	fprintf(fp,"NSTEPS		%d\n",(int) nsteps);//1200000000\n");

	/*TODO: Mutation rate is complicated. In the original paper, there were two rates:
	 * 		Indelrate: The rate of insertion and deletion
	 * 		Subrate: The rate of substitution.
	 * 	Later papers set these to be equal, read in by the MUTATE option. This means there are three possible
	 * 	settings:
	 * 		1: original rates - these rates are hard-coded, and set if there is no MUTATE parameter
	 * 		2: mutation rate - happens if a MUTATE value is set
	 * 		3: no mutation - happens with the line MUTATE 0
	 * 	There is also a parameter called 'domut' which can be set to zero at anytime to turn mutation off.
	 *
	 */


	if((indelrate - subrate)<FLT_EPSILON){
		if(indelrate<FLT_EPSILON)
			fprintf(fp,"MUTATE		0\n");
		else
			fprintf(fp,"MUTATE		%f\n",indelrate);
	}
	else{
		//No mutation rate needs to be set: the hard-wired alife values will be loaded.
	}


	/*TODO: Need to distinguish between 'USING' and 'SUBMAT' configurations.
	 * USING: files have the '.mis' extension
	 * SUBMAT: files have the '.mtx' extension
	 * Also need to know what to do if these are missing - do we create a separate file?
	 */
	if(strstr(swt_fn,".mtx")!=NULL){
		fprintf(fp,"SUBMAT	%s\n\n",swt_fn);
	}
	else{
		if(strstr(swt_fn,".mis")!=NULL){
			fprintf(fp,"USING	%s\n\n",swt_fn);
		}
		else{
			printf("Warning! no SWT substitution matrix specified\n");
		}
	}

	s_ag *spp,*pa;
	int spc,count;
	int finished = 0;
	int nag,*done;
	int i,found;

	nag = nagents(nowhead,-1);

	done = (int *) malloc(nag*sizeof(int));
	memset(done,0,nag*sizeof(int));

	do{
		i = 0;
		found=0;
		finished = 1;
		for(i=0,pa=nowhead;i<nag;i++,pa=pa->next){
			if(!done[i]){
				if(!found){
					done[i]=1;
					count=1;
					finished=0;
					found=1;
					spp = pa;
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

		/*
			% Basic seed replicase should look something like:
			AGENT OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B>C$=?>$$BLUBO%}OYHOB 150 Q
		*/
		//Write to file
		if(!finished){
			fprintf(fp,"%%%%%% SPECIES %d\n",spp->spp->spp);
			fprintf(fp,"AGENT\t%s\t%d\tQ\n\n",spp->spp->S,count);
		}

	}while(!finished);

	//Now write the reacting molecules (could be a bit tricky this)

	int nreactions = 0;
	for(i=0,pa=nowhead;i<nag;i++,pa=pa->next){
		if(pa->status == B_ACTIVE){
			fprintf(fp,"%%%%%% REACTION %d\nREACTION\n",++nreactions);
			print_agent_cfg(fp, pa);
			print_agent_cfg(fp,pa->pass);
		}
	}



	free(done);

	return 0;

}
