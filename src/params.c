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
#include "params.h"
/* This holds parameters */


int read_flag(FILE *fp,const char *label, int verbose){

	const int maxl = 128;
	char line[maxl];

	rewind(fp);
	while((fgets(line,maxl,fp))!=NULL){
		if(!strncmp(label,line,strlen(label))){
			if(verbose)
				printf("%s set\n",label);
			return 0;
		}
	}

	return 1;
}


void report_param_error(int error, int doexit){
	switch(error){
	case 0://Do nothing - no error
		break;
	case 1:
		printf("PARAMETER ERROR - Parameter not found\n");
		break;
	case 2:
		printf("PARAMETER ERROR - Duplicate parameter found\n");
		break;
	default:
		printf("PARAMETER ERROR - unspecified problem\n");
		break;
	}
	if(error && doexit){
		exit(error);
	}
}



int read_param_int(FILE *fp,const char *label, unsigned int *val, int verbose){

	const int maxl = 128;
	char line[maxl];
	int found=0;
	int err = 1;

	rewind(fp);
	while((fgets(line,maxl,fp))!=NULL){
		if(!strncmp(label,line,strlen(label))){
			if(!found){
				sscanf(line,"%*s %d",val);
				if(verbose)
					printf("%s = %d\n",label,*val);
				found = 1;
				err = 0;
			}
			else{ //found 2 labels - bad error!
				printf("found a second label!\n");
				fflush(stdout);
				getchar();
				err = 2;
			}
		}
	}
	if(!found)err=1;
	return err;
}


int read_param_float(FILE *fp,const char *label, float *val, int verbose){

	const int maxl = 128;
	char line[maxl];
	int found=0;
	int err = 0;

	rewind(fp);
	while((fgets(line,maxl,fp))!=NULL){
		if(!strncmp(label,line,strlen(label))){
			if(!found){
				sscanf(line,"%*s %f",val);
				if(verbose)
					printf("%s = %f\n",label,*val);
				found = 1;
			}
			else{ //found 2 labels - bad error!
				printf("found a second label!\n");
				fflush(stdout);
				getchar();
				err = 2;
			}
		}
	}
	if(!found)err=1;
	return err;
}



char * read_param_string(FILE **pfp,const char *label, int verbose){

	const int maxl = 128;
	char line[maxl];
	char st[maxl];
	char *s;
	s = NULL;
	int found=0;
	FILE *fp;

	fp = *pfp;

	memset(st,0,maxl*sizeof(char));

	while(((fgets(line,maxl,fp))!=NULL)&&!found){
		if(!strncmp(label,line,strlen(label))){
			if(!found){
				sscanf(line,"%*s %s",(char *) &st);
				if(verbose)
					printf("%s = %s\n",label,st);
				found = 1;
			}
		}
	}
	*pfp = fp;
	if(!found)
		return NULL;
	else{
		//TODO: checking needed on size of string!
		s = (char *) malloc(sizeof(char)*256);//(strlen(st)+1));
		memset(s,0,256*sizeof(char));
		strncpy(s,st,strlen(st));
		return s;
	}
}


int readordef_param_int(const char *fn, const char *label, unsigned int *val, const int defaultvalue, const int verbose){
	int rerr=1;
	int errcode = 3;
	FILE *fp;

	if((fp=fopen(fn,"r"))!=NULL){
		rerr = read_param_int(fp,label,val,0);
		switch(rerr){
		case 2:
			printf("Multiple instances of %s specified. Check config file\n",label);
			getchar();
			exit(0);
		case 0:
			if(verbose)printf("Setting %s to %d\n",label,*val);
			errcode = 0;
			break;
		default:
			if(defaultvalue < 0){
				errcode = 1;
			}
			else{
				*val = defaultvalue;
				errcode = 0;
				if(verbose)printf("%s not specified;\nSetting %s to %d\n",label,label,*val);
			}
			break;
		}
		fclose(fp);
	}
	if(verbose)printf("%s = %d\n",label,*val);

	return errcode;
}
