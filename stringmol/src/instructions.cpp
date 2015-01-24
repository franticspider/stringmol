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

#include "memoryutil.h"
#include "randutil.h"

//string stuff
#include "stringmanip.h"
#include "alignment.h"


//extern const int  maxl;

//Use this to control whether h-search is stochastic or not
#define SOFT_SEARCH

int LabLength(char *ip, const int maxl){

	int len=0;
	ip = ip+1;
	while(*ip > 64 && *ip <91){
		len++;
		ip++;
	}
	if(len>maxl){
		printf("Label = %d, longer than maxl (= %d!!\n",len,maxl);
	}
	return len;
}


////////////////////
// $: H-Search    //
////////////////////
char * HSearch(char *iptr, char *sp, swt *T, int *itog, int *ftog,const int maxl){

	char *ip,*tp,tmp[maxl];
	ip = iptr;
	int i,len=0;
	align A;

	memset(tmp,0,maxl*sizeof(char));

	/*NOTE: We are currently searching from the start of the string with the active flow pointer.
	Perhaps we should start AT the flow pointer, and "loop around" to the beginning of the string if
	there is no match in the first part.

	So the string:

	ABCDEFGHIJKLMNOPQRSTUVWXYZ
	          f

	would be serached as if it was:

	JKLMNOPQRSTUVWXYZABCDEFGHI

	..the best match along this line would be returned. Note it is possible that we could position <f> at the
	start of the line using this technique with a little modification. But probably better to implement a "decrement"
	operator */

	/*
	//first get the length of the string
	while(*ip > 64 && *ip <91){
		len++;
		ip++;
	}
	*/
	len = LabLength(ip, maxl);
	tp = iptr+len;
	ip=iptr+1;

	if(!len){
		//Ensure that the toggles are set:
		*ftog = *itog;
		return iptr;
	}

	memset(tmp,0,128*sizeof(char));
	strncpy(tmp,iptr+1,len);
	//generate the complement:
	for(i=0;i<len;i++)
		tmp[i] = AlphaComp(tmp[i]);

	//SmithWaterman(tmp,sp,&A,T,0);
	float bprob = SmithWatermanV2(tmp,sp,&A,T,0);

	//TODO: this will always match if any symbols match. There is no stochastic element..
#ifndef SOFT_SEARCH
	if(A.match)
#else
	int l = A.e1-A.s1 < A.e2-A.s2 ? A.e1-A.s1 : A.e2-A.s2;
	if(l<=2)
		bprob=0;
	else
		bprob = pow(A.score,l)/pow(l,l);


	float s = A.score<l-1.124? A.score : l-1.124;
	bprob = s/(l-1.124);

	float rno = rand0to1();
	if(rno<bprob)//search success!
#endif
		return tp + A.e2 - (tp-sp);

	//Ensure that the toggles are set - if no match found, we currenty move *F to *I - might leave it on the opposite string if it was there in the 1st place...:
	*ftog = *itog;
	return tp;
}



////////////////////
// ?: If-Label    //
////////////////////
char * IfLabel(char *ip, char *rp, char *sp, swt *T, const int maxl){

	char tmp[maxl],tmp2[maxl];
	int i,len = LabLength(ip, maxl);
	ip++;
	align A;

	switch(len){

	case 0:
		if(!*rp)
			return ip+len+1;
		return ip+len;

		break;


	case 1: //Possibly switch to look at different heads here....not implemented yet....
		if(!*rp)
			return ip+len+1;
		return ip+len;
		break;

	default:

		memset(tmp ,0,maxl*sizeof(char));
		memset(tmp2,0,maxl*sizeof(char));
		strncpy(tmp,ip,len);
		//generate the complement:
		for(i=0;i<len;i++)
			tmp[i] = AlphaComp(tmp[i]);

		strncpy(tmp2,rp,len);

		//SmithWaterman(tmp,tmp2,&A,T,0);
		SmithWatermanV2(tmp,tmp2,&A,T,0);

		if(align_event(&A,len))
			return ip+len+1;
		return ip+len;
		break;
	}

	//TODO: We should never reach this point, but check:
	return 0;
}

