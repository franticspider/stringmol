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


/* This is a header file for mt19937-2.c, which allows the functions in */
/* it to be called from other files. It can be used in conjunction with */
/* randutil.c and randutil.h which should appear along with this file   */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MT199372_H_
#define MT199372_H_

//FUNCTION HEADERS FOR RANDOM NUMBER GENERATOR
void	sgenrand(unsigned long seed);
double 	genrand(); /* generating reals */
unsigned long genrandint();/* unsigned long */ /* for integer generation */

int mt_get_mti();
void mt_set_mti(int val);

/*read/write mt state */
void print_mt(FILE *fp);
int load_mt(char *fn);

/* ERROR CODES FOR LOADING THE MT STATE */
enum load_mt_errcode{
	load_mt_nofile,
	load_mt_bad_mti,
	load_mt_bad_mt,
	load_mt_success
};

#endif /* MT199372_H_ */

#ifdef __cplusplus
}
#endif
