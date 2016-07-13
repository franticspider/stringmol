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

#ifndef PARAMS_H_
#define PARAMS_H_

int read_param_float(FILE *fp, const char *label, float *val, int verbose);
int read_param_int(FILE *fp,const char *label, unsigned int *val, int verbose);
char * read_param_string(FILE **pfp,const char *label, int verbose);

int readordef_param_int(char *fn, const char *label, unsigned int *val, const int defaultvalue, const int verbose);

void report_param_error(int error, int doexit);

#endif /*PARAMS_H_*/

#ifdef __cplusplus
	}
#endif
