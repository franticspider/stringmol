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


#ifdef __cplusplus
extern "C" {
#endif

#ifndef MEMUTIL_H_
#define MEMUTIL_H_

void memerror();
void * mymalloc(const int number, const int size);

//void ***arr3alloc(const int size, const int n1, const int n2, const int n3);
void	arr3free(void ***aaa, const int n1, const int n2, const int n3);
void	*** arr3alloc(const int n1, const int n2, const int n3, const int size);

void	** arr2alloc(const int n1, const int n2, const int size);
void	arr2free(void **aa, const int n1, const int n2);

int ** intarr2alloc(const int n1, const int n2);
void  intarr2free(int **aa, const int n1, const int n2);


int ***	intarr3alloc(const int n1, const int n2, const int n3);
void	intarr3free(int ***aaa, const int n1, const int n2, const int n3);


float ***	fltarr3alloc(const int n1, const int n2, const int n3);
void  fltarr3free(float ***aaa, const int n1, const int n2, const int n3);
//void	fltarr3free(double ***aaa, const int n1, const int n2, const int n3);

double ***	dblarr3alloc(const int n1, const int n2, const int n3);
void	dblarr3free(double ***aaa, const int n1, const int n2, const int n3);


int **** intarr4alloc(const int n1,const int n2,const int n3,const int n4);

void intarr4free(int ****aa, const int n1,const int n2,const int n3,const int n4);



//not strictly memutil - but util nevertheless...

float maxfinarr(float *arr, const int len);
double maxdinarr(double *arr, const int len, int *idx);



#endif /*MEMUTIL_H_*/

#ifdef __cplusplus
}
#endif


