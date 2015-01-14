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
//#include <string.h>
//#include <math.h>

//#include "memoryutil.h"
//#include "randutil.h"
//#include "params.h"

//string stuff
//#include "stringmanip.h"
#include "alignment.h"
//#include "instructions.h"

//metabolism stuff
#include "rules.h"
#include "agents_base.h"
#include "SMspp.h"
#include "stringPM.h"




void set_energy(const int NCON, int epercon, stringPM **A, char *signal, int *earray, float *score){


	//The idea is that we distribute the energy according to the signal score of each container.
	//float score[NCON];


	float sumsc=0;
	for(int i=0;i<NCON;i++){
		score[i]=0;

		s_ag *p;

		p = A[i]->nowhead;
		while(p!=NULL){
			if(p->status==B_UNBOUND){
				score[i] += p->spp->sig_sc;
			}
			p = p->next;
		}
		//score[i] /= A[i]->nagents(A[i]->nowhead,-1);
		sumsc+=score[i];
	}

	float toten = epercon * NCON;
	for(int i=0;i<NCON;i++){
		//score[i] /= sumsc;
		A[i]->energy += (float) toten * score[i] / sumsc;
		if(earray!=NULL){
			earray[i] = (float) toten * score[i] / sumsc;
		}
	}


}

