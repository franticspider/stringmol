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


#ifndef SMSPP_H_
#define SMSPP_H_


enum s_bind{B_UNBOUND,B_PASSIVE,B_ACTIVE};

//Forward declare the parent structure
struct s_parent;
struct l_spp;

struct s_ag{//THIS DEFINES AN INDIVIDUAL AGENT IN A STRINGMOL SYSTEM
	char 	*S;		//The executing string;
	//char 	*comp;	//the complement of the string
	int		label;
	int 	idx;	//the index of the agent.
	int		len;	//the length of the agent string.
	int 	ect;	//count of the number of instructions executed on bind
	int		nbind;	//number of binds
	s_bind 	status;
	char *i[2];		//instruction pointer
	char *r[2]; 	//read pointer
	char *w[2]; 	//write pointer
	char *f[2];		//flow (loop) pointer
	//Toggles: passive == 0; active == 1.
	int	it;			//intruction pointer toggle
	int rt;			//read pointer toggle
	int wt;			//write pointer toggle
	int ft;			//flow (loop) pointer toggle
	s_ag		*exec;
	s_ag		*pass;
	s_ag		*next;
	s_ag		*prev;

	//Lineages:
	l_spp 		*spp;
	s_parent 	*pp;

	//biomass
	int			biomass; //the biomass created during a reaction

	//spatial
	bool 		set;
	int 		x;
	int			y;
};

/*
struct s_spp{//THIS DEFINES WHAT A MOLECULAR SPECIES IS - can be pointed to by an agent OR just be in a list
	char 		*S;		//The string of the species
	int  		spp;	//The index number
	//int		paspp;	//The index of the active parent
	//int		ppspp;	//The index of the passive parent.
	s_parent	*pp;	//Pointer to the set of parent reactions
	//s_spp 		*next;
	char  		sptype;	//'A' = new from an active bind; 'P' = new from a passive bind; 'C' = new from a cleave bind
	int 		anc;	//A processing flag, specifically to flag an ancestor - reset when beginning an analysis!
};
*/

struct l_spp{//THIS IS THE LIST OF SPECIES THAT EVOLVE
	//s_spp 	*sp;	//Pointer to details of the species
	//int 		first;	//The first reported generation that the species appeared in...

	int 		count;  //The number of individuals in the lineage
	int			pf;		//processing flag = reset this when beginning an analysis!
	l_spp 		*next;	//Pointer to next element in the list
	int			tspp;	//Time of speciation - the first time the species was produced...

	//The following were from the s_spp typedef above
	char 		*S;		//The string of the species
	int  		spp;	//The index number
	s_parent	*pp;	//Pointer to the set of parent reactions
	char  		sptype;	//'A' = new from an active bind; 'P' = new from a passive bind; 'C' = new from a cleave bind
	int 		anc;	//A processing flag, specifically to flag an ancestor - reset when beginning an analysis!

	//A Score for the signal molecule
	float 		sig_sc;

	//Total biomass produced in this reaction over a trial
	int			biomass;

};

struct s_parent{
	l_spp 		*pa;	//The active parent
	l_spp 		*pp;	//The passive parent
	//int 		run;	//The trial in which the reaction was first seen
	int			n;		//The number of times the reaction has been seen
	s_parent 	*next;	//Pointer to the next data structure
};

class SMspp {
public:
	l_spp *species; 		//The list of species
	long int spp_count; 	//The current number of species

	SMspp();
	virtual ~SMspp();

	l_spp * make_spp(s_ag *a, int extit, const int maxl0);
	void 	prepend_spp(l_spp *p);
	l_spp * getspp(s_ag *a, int extit,const int maxl0);//, s_spp * paspp, s_spp * ppspp);
	//s_spp * getspp(s_ag *a);
	void free_spp(l_spp *sp);

	s_parent * make_parents(l_spp * paspp, l_spp * ppspp);
	s_parent * get_parents(l_spp * c, l_spp *paspp, l_spp  *ppspp);
	void 		append_parents(l_spp *c, s_parent *pp);
	void 		free_parent_list(s_parent *pp);

	int print_spp_list(FILE *fp);

	int clear_list();

	//TODO: write these - save memory!
	//int		print_list();
	//int		list_to_file(const char *fn, const char *mode);


};

#endif /* SMSPP_H_ */
