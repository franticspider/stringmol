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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "memoryutil.h"

#include "rules.h"


rules::rules(char *fn){

	const int maxl = 128;
	FILE *fp;

	if((fp=fopen(fn,"r"))!=NULL){
	    char line[maxl];
	    char label[maxl];
	    int rulecount = 0;
		while((fgets(line,maxl,fp))!=NULL){
			memset(label,0,maxl);
			sscanf(line,"%2000s",label);
			printf("line = %s",line);
			if(!strncmp(label,"RULE",4))
				rulecount++;
		}
		printf("\nfinished scanning - %d rules - rewinding\n",rulecount);
		//allocate memory for the rules now
		nr = rulecount;
		rset = (int **) arr2alloc(nr,4,sizeof(int));
		rval = (float *) mymalloc(nr,sizeof(float));
		re = (float *) mymalloc(nr,sizeof(float));

		//fill the rule array
		rewind(fp);
		rulecount = 0;
		while((fgets(line,maxl,fp))!=NULL){
			memset(label,0,maxl);
			sscanf(line,"%2000s",label);
			if(!strncmp(label,"RULE",4)){//fill the rule data
				sscanf(line,"%*s %c %c %c %c %f %f",(char * )&(rset[rulecount][0]),(char * )&(rset[rulecount][1]),(char * )&(rset[rulecount][2]),(char * )&(rset[rulecount][3]),&(rval[rulecount]),&(re[rulecount]));
				printf("RULE: %c %c %c %c %f %f\n",(rset[rulecount][0]),(rset[rulecount][1]),(rset[rulecount][2]),(rset[rulecount][3]),(rval[rulecount]),(re[rulecount]));
					rulecount++;
			}
		}

		//close
		fclose(fp);

	}
	else{
		printf("Unable to open file %s\n",fn);
		fflush(stdout);
	}
}

int rules::findinrange(int row, int min, int max, int label){
	int i;
	for(i=min;i<=max;i++){
		if(rset[row][i]==label)
			return 1;
	}
	return 0;
}

//says if a label is present in any rule at a specified side of the equations
//found("lhs",bag->label)
int rules::found(const char *side, int label){
	int i,sidecheck=0;
	if(!strncmp(side,"lhs",3)){
		sidecheck = 1;
		for(i=0;i<nr;i++){
			if(findinrange(i,0,1,label))
				return 1;
		}
	}
	if(!strncmp(side,"rhs",3)){
		sidecheck = 1;
		for(i=0;i<nr;i++){
			if(findinrange(i,2,3,label))
				return 1;
		}
	}

	if(!sidecheck){
		printf("\"%s\" is not a known specifer string. Use \"lhs\" or \"rhs\" instead\n",side);
	}
	return 0;
}



int rules::getrule(const char *side,int l1, int l2){
	int sidecheck = 0;
	if(!strncmp(side,"lhs",3)){
		sidecheck = 1;
		for(int i=0;i<nr;i++){
			if(rset[i][0]==l1 && rset[i][1]==l2)
				return i;
			else if (rset[i][1]==l1 && rset[i][0]==l2)
				return i;
		}
	}

	if(!sidecheck){
		printf("\"%s\" is not a known specifer string. Use \"lhs\" or \"rhs\" instead\n",side);
	}
	return -1;
}





rules::~rules(){

}


