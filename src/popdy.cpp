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

#include "popdy.h"





int popdy::load_popdy(char *fn){


	FILE *fp;
	const int ml = 256;
	char line[256];
	int lcount = 0;

	long time;
	char *spp;
	int count;

	pd_time * pt;
	pd_spp * psp;


	if((fp = fopen(fn,"r"))!=NULL){


		spp = (char *) malloc(ml*sizeof(char));
		while((fgets(line,ml,fp))!=NULL){
			sscanf(line,"%l,%s,%d",&time,spp,&count);

			pt  = add_time(time);
			psp = add_spp(pt,spp);
			lcount++;
		}

		fclose(fp);
		free(spp);
	}
	else
		printf("Unable to open file %s\n",fn);

	return 0;
}



struct pd_time * popdy::add_time(long time){
	pd_time * pt;

	if(time_head == NULL){
		time_head = (pd_time *) malloc(sizeof(pd_time));

		time_head->time = time;
		return time_head;

	}
	else{
		for(pt = time_head; pt!=NULL; pt = pt->next){
			if(pt->time == time)
				return pt;
		}
		//if we've made it to here, we need a new time value
		pt->next = (pd_time *) malloc(sizeof(pd_time));
	}

	return pt;
}


struct pd_spp * popdy::add_spp(pd_time *pt, char *spp){
	pd_spp *psp;

	psp = NULL;

	return psp;
}
