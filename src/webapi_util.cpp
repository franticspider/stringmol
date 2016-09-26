/* Copyright (C) 2009-2015 Simon Hickinbotham, Adam Nellis, Ed Clark    */
/* When you use this, send an email to: sjh436@gmail.com                */
/* with an appropriate reference to your work.                          */

/* This file is part of STRINGMOL					*/

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

//setup
// Writing PNGs
#include "lodepng.h"
#include <iostream>
#include "setupSM.h"

#include "webapi_util.h"

//extern const int maxl;//150;// 512;
//extern const int maxl0; //allow room for a terminating 0


//#define DEBUG

const int MAXSTR = 5000;


//TODO: This belongs elsewhere...possibly in stringPM itself
/** construct a stringmol from a 0-terminated string */
s_ag * make_mol(stringPM *A, char* string){

	l_spp *species;
	s_ag *pag;


	//PRE-PROCESS THE STRINGS...
	/**
	 * Convert '1' to '>'....
	 */
	char *p;
	for(p=&(string[0]);*p!='\0';p++){
		if(*p == '1')
			*p='>';
	}


	pag = A->make_ag('X');

	pag->S =(char *) malloc(A->maxl0*sizeof(char));

	memset(pag->S,0,A->maxl0*sizeof(char));
	strncpy(pag->S,string,strlen(string));
	pag->len = strlen(pag->S);

	//No parents for these initial agents!
	pag->pp = A->spl->make_parents(NULL,NULL);


	//add the mols to the species list
	species = A->spl->getspp(pag,0,A->maxl0);
	species->tspp = 0;

	pag->spp=species;

	return pag;
}



char * generate_bind_data(char * string1, char * string2){

	int 	bindImpossible=0,
			bindLength=0,
			start1=0,
			end1=0,
			start2=0,
			end2=0,
			//pointers:
			instr1=0,
			instr2=0,
			flow1=0,
			flow2=0,
			read1=0,
			read2=0,
			write1=0,
			write2=0,
			instrToggle=0,
			flowToggle=0,
			readToggle=0,
			writeToggle=0;
	float	bindProb=1.0;
	char * output;

	output = (char *)malloc(MAXSTR*sizeof(char));
	memset(output,0,MAXSTR*sizeof(char));


#ifdef DEBUG
	/* Show that the conversion worked */
	printf("<p>after conversion </p><p> string1 = \"%s\"</p><p> string2 = \"%s\"</p>\n",string1,string2);
#endif

	/*Run the bind*/
	SMspp		SP;
	stringPM 	A(&SP);
	align sw;
	s_ag* a1,*a2,*pact;
	 //Load default table
	A.blosum = default_table();
	 //Construct the molecules
	a1 = make_mol(&A,string1);
	a2 = make_mol(&A,string2);
	 //Get the Smith Waterman alignment
	A.get_sw(a1,a2,&sw);
	 //Figure out the active molecule
	A.set_exec(a1,a2,&sw);
	if(a1->status == B_ACTIVE)
		pact = a1;
	else
		pact = a2;
	 //Get the bind probability
	float bprob = A.get_bprob(&sw);


	/*TODO: Questions for Adam:
	 * 1: Below what value of bindProb do we say that bindImpossible is TRUE?
	 * 2: How is bind length calculated? is it min(end1-start1,end2-start2), or something else?
	 */


	if(bprob<0.00001)
		bindImpossible=1;
	else{
		bindProb = bprob;//sw.prob;
		bindLength = (sw.e1-sw.s1)>(sw.e2-sw.s2)?(sw.e1-sw.s1):(sw.e2-sw.s2);
		start1 = sw.s1;
		end1 = sw.e1;
		start2 = sw.s2;
		end2 = sw.e2;
		//string1 is already set
		//string2 is already set
		instr1 = pact->i[1] - pact->S;
		instr2 = pact->i[0] - pact->pass->S;
		flow1 =  pact->f[1] - pact->S;
		flow2 =  pact->f[0] - pact->pass->S;
		read1 =  pact->r[1] - pact->S;
		read2 =  pact->r[0] - pact->pass->S;
		write1 = pact->w[1] - pact->S;
		write2 = pact->w[0] - pact->pass->S;
		instrToggle = pact->it;
		flowToggle	= pact->ft;		//19
		readToggle	= pact->rt;		//21
		writeToggle	= pact->wt;		//22
	}


	//Create the output
	//                      1             2           3            4       5         6       7          8          9         10        11       12       13       14       15        16        17             18            19            20             21
	sprintf(output,
			"bindImpossible=%d&bindLength=%d&bindProb=%0.3f&start1=%d&end1=%d&start2=%d&end2=%d&string1=%s&string2=%s&instr1=%d&instr2=%d&flow1=%d&flow2=%d&read1=%d&read2=%d&write1=%d&write2=%d&instrToggle=%d&flowToggle=%d&readToggle=%d&writeToggle=%d",
			bindImpossible,	//1
			bindLength,		//2
			bindProb,		//3
			start1,			//4
			end1,			//5
			start2,			//6
			end2,			//7
			string1,		//8
			string2,		//9
			instr1,			//10
			instr2,			//11
			flow1,			//12
			flow2,			//13
			read1,			//14
			read2,			//15
			write1,			//16
			write2,			//17
			instrToggle,	//18
			flowToggle,		//19
			readToggle,		//21
			writeToggle		//22
			);
	return output;

}




/*TODO: This is what the output string data needs to be:
'string1=('+$(".stringActive", currentState).html().replace(/&gt;/g, "1")+')&'+
'string2=('+$(".stringPassive", currentState).html().replace(/&gt;/g, "1")+')&'+
'instr1=('+$(".instr1", currentState).html()+')&'+
'instr2=('+$(".instr2", currentState).html()+')&'+
'flow1=('+$(".flow1", currentState).html()+')&'+
'flow2=('+$(".flow2", currentState).html()+')&'+
'read1=('+$(".read1", currentState).html()+')&'+
'read2=('+$(".read2", currentState).html()+')&'+
'write1=('+$(".write1", currentState).html()+')&'+
'write2=('+$(".write2", currentState).html()+')&'+
'instrToggle=('+$(".instrToggle", currentState).html()+')&'+
'flowToggle=('+$(".flowToggle", currentState).html()+')&'+
'readToggle=('+$(".readToggle", currentState).html()+')&'+
'writeToggle=('+$(".writeToggle", currentState).html()+')',
*/
char * generate_step_data(char * string1, char * string2, int instr1, int instr2, int flow1, int flow2, int read1, int read2, int write1, int write2, int instrToggle, int flowToggle, int readToggle, int writeToggle){

	char * output;

	output = NULL;

	//Run the bind
	SMspp		SP;
	stringPM 	A(&SP);
	s_ag *a1,*a2;


	//Load default table
	A.blosum = default_table();

	//Construct the molecules
	a1 = make_mol(&A,string1);
	a2 = make_mol(&A,string2);

	//TODO: We are assuming that mol1 is active - this may not be the case...
	a1->status = B_ACTIVE;
	a2->status = B_PASSIVE;

	a1->i[1] = &(a1->S[instr1]);
	a1->i[0] = &(a2->S[instr2]);

	a1->f[1] = &(a1->S[flow1]);
	a1->f[0] = &(a2->S[flow2]);

	a1->r[1] = &(a1->S[read1]);
	a1->r[0] = &(a2->S[read2]);

	a1->w[1] = &(a1->S[write1]);
	a1->w[0] = &(a2->S[write2]);

	a1->it = instrToggle;
	a1->ft = flowToggle;
	a1->rt = readToggle;
	a1->wt = writeToggle;

	//Execute the next step:
	A.exec_step(a1,a2);

	//Now write the string

/*TODO: Finish this function...



		instr1 = pact->i[1] - pact->S;
		instr2 = pact->i[0] - pact->pass->S;
		flow1 =  pact->f[1] - pact->S;
		flow2 =  pact->f[0] - pact->pass->S;
		read1 =  pact->r[1] - pact->S;
		read2 =  pact->r[0] - pact->pass->S;
		write1 = pact->w[1] - pact->S;
		write2 = pact->w[0] - pact->pass->S;
		instrToggle = pact->it;
		flowToggle	= pact->ft;		//19
		readToggle	= pact->rt;		//21
		writeToggle	= pact->wt;		//22
*/
	
	return output;
}





