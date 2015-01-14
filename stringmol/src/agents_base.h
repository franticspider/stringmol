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



#ifndef AGENT_H_
#define AGENT_H_


struct s_ix{
	int n;
	float prob;
	int start;
	int stop;
	int	label;
	s_ix *next;
};


class agents_base{

	private:

		//s_ix *ifxhead;

		//virtual s_ag * make_ag(int alab, int randpos)=0;


	public:

		//for testing propensity
		int  *bct; //count of times we've tried for a particular reagent count
		int  *bpp; //count of times we've been closenough to run a bind...
		int  bmax; //max reagents we are going to bother with.
		void 	note_propensity(int N,int X);
		void 	print_propensity(FILE *fp);
		int 	proper_prop(const int n);
		int 	eqn_prop(const int n);

		s_ix *ifxhead;
		//s_ag *nowhead;
		//s_ag *nexthead;

		float cellrad;  //The radius of the "cell"
		float agrad;    //The active radius of the agent
		float vcellrad;	//The radius as a function of the number of div count things
		float move;		//The amount an agent can move. Set to 1.1 * agrad...
		long energy;	//The initial energy present in the system
		float nsteps;	//The number of steps (SHOULD BE AN INT!)

		//division
		//float divtime;	//The time to divide (SHOULD BE AN INT!)
		//int	*adc;		//Agent division count  - if all these conditions are met, divide.
		//int *aro;		//list of operons to replenish at division

		//agent state tables
		//int ntt;  		//number of agent types we are tracking.
		//int *aat; 		//array of agent types
		//int *aac; 		//array of agent counts
		//int **btab; 	//table of binding pairs
		//int **com; 		//table of agent complexes for replenishment at division
		//int **dcom;		//table of agent complexes for division counts

		//diagnostics
		int ict;		//count of influx rules
		int *irt;		//influx rule tested
		int *irf;   	//influx rule fired
		int *irl;   	//influx rule label
		int *tr; 		//count of the number times a rule has been tried
		int *fr;		//count of the number times a rule has fired.

		//creators and destructors
		agents_base();
		~agents_base();
		void preset();
		void clearout(int verbose);


		//fromfile stuff
		void load(char *fn, char *fninput, int test, int verbose);
		void load(char *fn, int test, int verbose);
		int test(char *fn);
		void test2();

		int load_params(char *fn, int test, int verbose);
		int load_influx(char *fn);
		int load_division(char *fn);
		int load_replenish(char *fn);
		//declaring this as virtual and calling it from load caused problems...
		//the `=0' is key!
		//virtual int load_agents(char *fn, int test, int verbose) = 0;
		virtual int load_agents(char *fn, char *fninput, int test=0, int verbose=0) = 0;// VJH - added this function



		//influx info
		s_ix * make_influx(int lab, int n, float prob, int start, int stop);
		//void influx(int i);
		int append_ix(s_ix **list, s_ix *ax);
		//void influx_special(int t);

		//Diagnoistics
		//void print_agents(const char *spec);
		void print_agents_header(FILE *fp);
		void print_agents_count(FILE *fp);
		//void update_aac();

		//Counting and indexing
		int aac_count(int lab);
		int nttindex(const int label);
		void make_btab(rules *rset);
		void make_com();
		void make_dcom();
		void makefr(int nr);
		void printfr(FILE *fp, rules *rset);

		//Iteration stuff
		//virtual void make_next(rules *rset)=0;
		virtual void make_next()=0;

		//Cell division
		int divide_conditions(int time);
		//void divide();
		//void replenish_operons();

		//Death!
		int hasdied();

};

#endif /*AGENT_H_*/
