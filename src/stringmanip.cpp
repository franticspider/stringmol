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
#include <string.h>

int AlphaComp(int t){

	int a;
	if(t>=65 && t<=90){ //n-ops
		a = t-65;
		a+=13;
		a = a%26;
		a = a+65;
	}
	else{//x-ops:
		a=t;
		/*
		switch(t){
		//flow
		case '$':
			a = '?';
			break;
		case '?':
			a = '$';
			break;

		//move
		case '>':
			a = '^';
			break;
		case '^':
			a = '>';
			break;

			//modify
		case '=':
			a = '%';
			break;
		case '%':
			a = '=';
			break;

		//end
		case '}':
			a = '}';
			break;
		}
		*/
	}
	return a;

}


char * string_comp(char *S){
	char *comp;
	int len = strlen(S)+1;
	unsigned int i;

	comp = (char *) malloc(len*sizeof(char));
	memset(comp,0,len*sizeof(char));

	for(i=0;i<strlen(S);i++)
		if(S[i]!=0)
			comp[i] = AlphaComp(S[i]);

	return comp;
}
