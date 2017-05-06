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


#ifndef STRINGPM_H_
#define STRINGPM_H_

enum e_mut{M_NONE,M_INCREMENT,M_DECREMENT,M_INSERT,M_DELETE};

/* Cellular Automata parameters */
enum s_gstatus{G_EMPTY,G_NOW,G_DOING,G_NEXT};

/* Reaction Load params */
//L_EVOEVO 		was only used for 2016 EvoEvo experiments
//L_REPLICABLE	is designed to have runs that can be restarted
enum s_loadtype{
	L_EVOEVO,
	L_REPLICABLE};


typedef struct td_smsprun{
	int 	 gridx;
	int		 gridy;
	s_ag *** grid; //2d grid of pointers to agents.
	enum s_gstatus	 **  status; //current status of cell
} smsprun;



class stringPM: public agents_base
{

private:

public:
	s_ag * make_ag(int alab);//, int randpos);

	s_ag *nowhead;
	s_ag *nexthead;

	//l_spp *species;
	SMspp *spl;

	int	spp_count;

	swt	*blosum;
	s_sw *swlist;

	long agct;
	unsigned int extit; 	//record of the cell iteration count

	unsigned long randseed;
	s_loadtype loadtype;		//Type of load we are doing (for backwards compatibility)

	long biomass; 	//used as a measure of fitness
	long bstart;    //time of biomass reset

	int domut;
	int dodecay;

	//Epoch recording
	int lastepoch;
	int thisepoch;
	int nepochs;

	//Conservation of Mass structures:
	int *mass; //This can be built and populated after "blosum" has been set...

	int run_number; //The run number
	char swt_fn[256];

	//Mutation
	float subrate;
	float indelrate;

	//Decay
	float decayrate;

	//Reporting flags
	int verbose_bind;
	int verbose_load;
	unsigned int splprint;		//the number of timesteps before printing out the specieslist

	//Reporting timings
	unsigned int report_every;			//How often to write splists and configs
	unsigned int image_every;			//How often to generate an image (spatial stringmol only)

	//Signal string - for getting a signal score:
	char *signal;

	//Max molecular lengths
	unsigned int maxl;
	unsigned int maxl0;

	//energy step
	unsigned int estep;

	//toggle for granularity
	int granular_1;

	//File name for the popdy file...
	char popdyfn[128];

	//linecount variable for loading
	int linecount;


	//constructor
	stringPM(SMspp * pSP);
	//destructor
	~stringPM();

	void preset();
	void clearout(int verbose = 0);

	//return a string describing the error
	char * parse_error(int errno);

	//Loading
	int load_splist(const char *fn,int verbose);

	//int load_agents(char *fn,                int test, int verbose);
	// VJH version for youShare
	int load_agents(const char *fn, char *fninput, int test, int verbose);

	float load_mut(const char *fn, int verbose); //load the mutation rate
	float load_decay(const char *fn, int verbose); //load the decay rate

	int load_reactions(const char *fn, char *fntab, int test, int verbose);

	int load_table_matrix(const char *fn);
	//int load(char *fn);


	//Aspatial testing:
	int aspatial_coverage(int n);
	int aspbind(rules *rset, s_ag *pag);

	//Iteration
	void make_next();
	int testbind(s_ag *pag);
	int testdecay(s_ag *pag);
	int hasdied();

	void replenish_operons();
	void divide();

	//list stuff
	int 	append_ag(s_ag **list, s_ag *ag);
	int 	extract_ag(s_ag **list, s_ag *ag);
	int 	nagents(s_ag *head, int state);
	s_ag * 	rand_ag(s_ag *head, int state);
	int 	free_ag(s_ag *pag);
	bool 	ag_in_list(s_ag *list, s_ag *tag);

	//First version works fine, but no species analysis...
	//void 	unbind_ag(s_ag * pag,char sptype);
	int 	unbind_ag(s_ag * pag, char sptype, int update, l_spp *pa, l_spp *pp);

	void update();

	//INSTRUCTION SET:
	int hcopy(s_ag *act); 	//	=	HCOPY
	int cleave(s_ag *act);  //	=	CLEAVE

    //Influx
	void influx_special(int t);
	void influx(int i);


    //list / array interface:
	void update_aac();

	//Diagnostics
	void print_agents(FILE *fp, const char *spec, int verbose);
	int print_agent_idx(FILE *fp, int det, int idx);
	void testprop();
	int countcomp(char E, char P);
	void sanity_check();


	//Distance calcs
	float close_dist(s_ag *a1, s_ag *a2);
	void move_ag(s_ag *a1);

	//Cellular Automaton stuff
	smsprun *grid;
	smsprun * init_smprun(const int gridx, const int gridy);
	void free_grid();
	void print_grid(FILE *fp);

	//Loading from config file (see also agents_base)
	int 	load_table(const char *fn);
	int 	load_replicable(const char *fn);

	s_ag * 	read_unbound_agent(FILE **fp, char line[], const int llen);
	s_ag *  read_active_agent(FILE **fp, char line[], const int llen, int &pidx);
	s_ag *  read_passive_agent(FILE **fp, char line[], const int llen);




	//String & alignment stuff
	int 	h_pos(s_ag *pag, char head); 	//Find the position of a particular head.
	void 	get_string_comp(s_ag *pag);
	float 	get_sw(s_ag *a1, s_ag *a2, align *sw);
	float 	get_bprob(align *sw);
	void 	set_exec(s_ag *A, s_ag *B, align *sw);
	int 	exec_step(s_ag *act, s_ag *pass);
			//print string stuff
	void 	print_ptr_offset(FILE *fp, char *S, char *p,int F, char c);
	void 	print_exec(FILE *fp, s_ag *act, s_ag *pas);
	void 	free_swt(swt *pSWT, int verbose);
	int 	check_ptrs(s_ag* act);
	int 	rewind_bad_ptrs(s_ag* pag);


	//Checking the energy model: (THIS RESULTS IN AN UNSTABLE SYSTEM)
	void 	energetic_make_next();
	int 	energetic_exec_step(s_ag *act, s_ag *pass);
	int 	energetic_testbind(s_ag *pag);



	//Trying conservation of mass
	int 	load_comass(const char *fn, int verbose); //load single value from a file
	int 	set_mass(int *param);  //load a set of values from an array
	void 	comass_make_next();
	int 	comass_testdecay(s_ag *pag);
	int 	comass_free_ag(s_ag *pag);
	int 	update_mass(char *S, int len, int val, const int doconcat);
	int 	comass_exec_step(s_ag *act, s_ag *pass);
	int 	comass_hcopy(s_ag *act);

	//Speigelman's monster

	int 	speig_hcopy(s_ag *act);


	//Molecular species analysis:
	//void 		update_lineage(s_ag *p,char sptype);
	int 		get_ecosystem();
	int 		update_lineage(s_ag *p, char sptype, int add, l_spp *paspp, l_spp * ppspp, int mass);
	void 		print_lineage_dot(FILE *fp, int time,int step); //traces everything descending from the initial set.
	void 		print_ancestry_dot(FILE *fp, int time,int step); //takes all current agents and traces them back
	void 		print_spp_strings(FILE *fp);
	//s_spp * 	make_spp(s_ag *a);
	//l_spp * 	make_lspp(s_ag *a);
	void 		free_lspp(l_spp *sp);
	l_spp * 	get_spp(int n);
	//int 		append_spp(s_spp *sp);
	int 		append_lspp(l_spp *sp);
	int 		count_spp();
	void 		print_spp_count(FILE *fp,int style, int state);
	void 		get_spp_count(int state);//Count the number of individuals of each species present in the system
	//find a species that a molecule belongs to
	int 		id_spp(l_spp *sp, s_ag *pag, int  aspno, char *spp_string);

	//Network analysis
	void 		get_spp_network(char *fn);
	int 		Network_cleave(s_ag *act);
	int 		Network_exec_step(s_ag *act, s_ag *pass);

	//Epoch recording
	void 		set_epochs();

	//Container functions
	int 		share_agents(s_ag **head); //Copies a set of agents onto nowhead;

	//Signal functions
	float 		sigalign(char *str);

	//Output a config file at any point in the trial:
	int 		print_conf(FILE *fp);
	void 		print_agent_cfg(FILE *fp, s_ag *pa, const int pass_index);

	//Write a species list to a config file
	void 		write_extant_spp(FILE *fp);

};


#endif /* STRINGPM_H_ */
