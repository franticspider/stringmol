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

#include "memoryutil.h"

void memerror(){

	printf("Memory allocation error\n");
	fflush(stdout);
	getchar();
	exit(6);
}


void * mymalloc(const int number, const int size){
	void *mem;
	if((mem = malloc(number*size))==NULL)
		memerror();
	return mem;
}


int ** intarr2alloc(const int n1, const int n2){

	int i;
	int **aa;

	aa=(int **) malloc(n1*sizeof(int *));
	for(i=0;i<n1;i++){
		aa[i]=(int *) malloc(n2*sizeof(int ));
	}

	return aa;
}

void  intarr2free(int **aa, const int n1, const int n2){

	int i;

	for(i=0;i<n1;i++){
		free(aa[i]);
	}
	free(aa);
}



void ***old_arr3alloc(const int size, const int n1, const int n2, const int n3){
	int i,j;
	int ***aa;

	//here we assume that pointers are always the same size regardless of type:

	aa=(int ***) malloc(n1*sizeof(int **));
	for(i=0;i<n1;i++){
		aa[i]=(int **) malloc(n2*sizeof(int *));
		for(j=0;j<n2;j++)
			aa[i][j] = (int *) malloc(n3*sizeof(size));
	}

	return (void ***) aa;

}




int **** intarr4alloc(const int n1,const int n2,const int n3,const int n4){
	int i,j,k;
	int ****aa;


	aa=(int ****) malloc(n1*sizeof(int ***));
	for(i=0;i<n1;i++){
		aa[i] = (int ***) malloc(n2*sizeof(int **));
		for(j=0;j<n2;j++){
			aa[i][j] = (int **) malloc(n3*sizeof(int *));
			for(k=0;k<n3;k++){
				aa[i][j][k] = (int *) malloc(n4*sizeof(int));
			}
		}
	}


	return aa;
}

void intarr4free(int ****aa, const int n1,const int n2,const int n3,const int n4){
	int i,j,k;

	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++){
			for(k=0;k<n3;k++){
				free(aa[i][j][k]);
			}
			free(aa[i][j]);
		}
		free(aa[i]);
	}
	free(aa);


}


int *** intarr3alloc(const int n1, const int n2, const int n3){

	int i,j;
	int ***aa;

	aa=(int ***) malloc(n1*sizeof(int **));
	for(i=0;i<n1;i++){
		aa[i]=(int **) malloc(n2*sizeof(int *));
		for(j=0;j<n2;j++)
			aa[i][j] = (int *) malloc(n3*sizeof(int));
	}

	return aa;
}


void  intarr3free(int ***aaa, const int n1, const int n2, const int n3){

	int i,j;

	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++)
			free(aaa[i][j]);
		free(aaa[i]);
	}
	free(aaa);
}





float *** fltarr3alloc(const int n1, const int n2, const int n3){

	int i,j;
	float ***aa;

	aa=(float ***) malloc(n1*sizeof(float **));
	for(i=0;i<n1;i++){
		aa[i]=(float **) malloc(n2*sizeof(float *));
		for(j=0;j<n2;j++)
			aa[i][j] = (float *) malloc(n3*sizeof(float));
	}

	return aa;
}


void  fltarr3free(float ***aaa, const int n1, const int n2, const int n3){

	int i,j;

	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++)
			free(aaa[i][j]);
		free(aaa[i]);
	}
	free(aaa);
}

void  arr3free(void ***aaa, const int n1, const int n2, const int n3){

	int i,j;

	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++)
			free(aaa[i][j]);
		free(aaa[i]);
	}
	free(aaa);
}


void	*** arr3alloc(const int n1, const int n2, const int n3, const int size){

	int i,j;
	void ***aa;

	aa=(void  ***) malloc(n1*sizeof(void **));
	for(i=0;i<n1;i++){
		aa[i]=(void **) malloc(n2*sizeof(void *));
		for(j=0;j<n2;j++)
			aa[i][j] = (void *) malloc(n3*size);
	}

	return aa;
}




void	** arr2alloc(const int n1, const int n2, const int size){

	int j;
	void **aa;

	aa=(void **) malloc(n1*sizeof(void *));
	for(j=0;j<n1;j++)
		aa[j] = (void *) malloc(n2*size);

	return aa;
}


void  arr2free(void **aa, const int n1, const int n2){

	int i;

	for(i=0;i<n1;i++){
		free(aa[i]);
	}
	free(aa);
}



double *** dblarr3alloc(const int n1, const int n2, const int n3){

	int i,j;
	double ***aa;

	aa=(double ***) malloc(n1*sizeof(double **));
	for(i=0;i<n1;i++){
		aa[i]=(double **) malloc(n2*sizeof(double *));
		for(j=0;j<n2;j++)
			aa[i][j] = (double *) malloc(n3*sizeof(double));
	}

	return aa;
}


void  dblarr3free(double ***aaa, const int n1, const int n2, const int n3){

	int i,j;

	for(i=0;i<n1;i++){
		for(j=0;j<n2;j++)
			free(aaa[i][j]);
		free(aaa[i]);
	}
	free(aaa);
}



double maxdinarr(double *arr, const int len, int *idx){
	int i,ii=0;
	double max = arr[0];
	for(i=0;i<len;i++){
		if(max<arr[i]){
			max = arr[i];
			ii=i;
		}
	}
	if(idx!=NULL)
		*idx = ii;
	return max;
}


float maxfinarr(float *arr, const int len){
	int i;
	float max = arr[0];
	for(i=0;i<len;i++){
		max = max<arr[i] ? arr[i] : max;
	}
	return max;
}
