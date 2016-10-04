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
#include "SMspp.h"


//extern const int  maxl0;

SMspp::SMspp() {
	spp_count=1;  // Why is this 1??
	species = NULL;
}

SMspp::~SMspp() {
	/* nothing needs doing here at present */
}


s_parent *SMspp::get_parents(l_spp * c, l_spp *paspp, l_spp  *ppspp){

	s_parent *pp;
	for(pp = c->pp;pp!=NULL;pp=pp->next){
		if(paspp == pp->pa)
			if(ppspp == pp->pp){
				pp->n++;
				return pp;
			}
	}

	//We'll only make it to here if a parent is not found
	pp = make_parents(paspp,ppspp);
	append_parents(c,pp);
	return pp;
}


void SMspp::append_parents(l_spp *c, s_parent *pp){
	s_parent *oo;

	if(c->pp==NULL){
		c->pp = pp;
		return;
	}
	else{
		oo = c->pp;
		while(oo->next != NULL){
			oo=oo->next;
		}
		oo->next = pp;
	}
}

//NB: We have to know what the parents are to do this, and to have checked that
// the parent pair does not already exist.
s_parent * SMspp::make_parents(l_spp * paspp, l_spp * ppspp){

	s_parent *pp;
	pp = (s_parent *) malloc(sizeof(s_parent));

	//if((pp->pa = getspp(paspp))==NULL)
	//		printf("OOPS! bad active parent\n");
	//if((pp->pp = getspp(ppspp))==NULL)
	//		printf("OOPS! bad passive parent\n");
	pp->pa = paspp;
	pp->pp = ppspp;
	pp->n=1;
	pp->next=NULL;

	return pp;
}


l_spp * SMspp::make_spp_from_string(char *S, int extit, const int maxl0, const int spno){

	int l;
	l_spp *sp;

	//TODO: check for malloc fails here...
	sp = (l_spp *) malloc(sizeof(l_spp));
	sp->S=(char *) malloc(maxl0+1*sizeof(char));

	memset(sp->S,0,maxl0*sizeof(char));

	//TODO: what if l > maxl0 ??
	l = strlen(S);
	strncpy(sp->S,S,l);

	sp->pp = NULL;
	if(spno>-1){
		//Todo - check that spno doesn't already exist:
		l_spp *lsp;
		for(lsp = species; lsp != NULL; lsp = lsp->next){
			if(lsp->spp == spno){
				printf("ERROR: species %d already exists! string is: %s\n",lsp->spp,lsp->S);
				return NULL; //todo handle the error more gracefully!
			}
		}
		sp->spp = spno;
		spp_count = spp_count>spno?spp_count:(spno+1);
	}
	else{
		sp->spp=spp_count++;
	}
	sp->sptype=0;
	//These are from make_lspp in the stringPM class:
	sp->count=0;
	sp->next=NULL;
	sp->tspp=extit;
	sp->biomass=0;

	return sp;

}


l_spp * SMspp::make_spp_from_agent(s_ag *a, int extit, const int maxl0){

	l_spp *sp;

	sp = make_spp_from_string(a->S,extit,maxl0,-1);
	//Here we assume that parents of the agent are already set!
	//sp->pp = make_parents(a->pp->pa,a->pp->pp);
	//a->pp = sp->pp;

	a->pp = NULL;
	return sp;
}


void SMspp::prepend_spp(l_spp *sp){

	sp->next = species;
	species = sp;

}

l_spp * SMspp::find_spp(char *S, const int maxl0){

	l_spp *p;
	//First, check if it's in the list:
	for(p=species;p!=NULL;p=p->next){
		if(!strcmp(p->S,S))
			return p;
	}

	return NULL;

}


/*gets a spp no. and makes a new one if needed */
l_spp * SMspp::getspp(s_ag *a, int extit,const int maxl0){//s_spp * paspp, s_spp * ppspp){

	l_spp *p;

	p = find_spp(a->S,maxl0);

	if(p==NULL){
		//If not found, make a new one, append it and return the address of the new.
		p = make_spp_from_agent(a,extit,maxl0);
		prepend_spp(p);
	}

	return p;
}

/*gets a spp no. and makes a new one if needed */
l_spp * SMspp::getspp_from_string(char *S, int extit,const int maxl0, const int spno){//s_spp * paspp, s_spp * ppspp){

	l_spp *p;

	p = find_spp(S,maxl0);

	if(p==NULL){
		//If not found, make a new one, append it and return the address of the new.
		p = make_spp_from_string(S,extit,maxl0,spno);
		prepend_spp(p);
	}

	return p;
}


/*
s_spp * SMspp::getspp(int spno){

	s_spp *p;
	//First, check if it's in the list:
	for(p=species;p!=NULL;p=p->next){
		if(spno == p->spp)
			return p;
	}

	//If not, make a new one, append it and return the address of the new.
	printf("ERROR - species %d not found\n",spno);
	fflush(stdout);
	return NULL;
}
*/

void SMspp::free_parent_list(s_parent *pp){

	s_parent *tmp;

	while(pp !=NULL){
		tmp=pp;
		pp=pp->next;
		free(tmp);
	}
}


void SMspp::free_spp(l_spp *sp){
	free(sp->S);
	free_parent_list(sp->pp);
	free(sp);
	//NEED ALSO TO FREE PARENT LIST!
}


int SMspp::clear_list(){
	l_spp *ps,*n;

	for(ps=species;ps!=NULL;){
		n = ps->next;
		free_spp(ps);
		ps = n;
	}
	species = NULL;
	spp_count=1;  // Why is this 1??
	return 0;
}


int SMspp::print_spp_list(FILE *fp){

	l_spp *ps;
	s_parent *pp;
	int count=0;

	if(fp==NULL)
		fp=stdout;

	for(ps=species;ps!=NULL;ps=ps->next){
			for(pp=ps->pp;pp!=NULL;pp=pp->next){
				count++;
				if(pp->pa && pp->pp)//Need to find a better way....
					fprintf(fp,"%d,%d,%d,%f,%d,%d,%d,%s\n",ps->spp, pp->pa->spp, pp->pp->spp ,ps->sig_sc,pp->n,ps->tspp,ps->biomass,ps->S);
				else{
					fprintf(fp,"%d,",ps->spp);
					if(pp->pa)
						fprintf(fp,"%d,",pp->pa->spp);
					else
						fprintf(fp,"-1,");
					if(pp->pp)
						fprintf(fp,"%d,",pp->pp->spp);
					else
						fprintf(fp,"-1,");
					fprintf(fp,"%d,%s\n",pp->n,ps->S);
				}
			}
	}
	return count;
}
