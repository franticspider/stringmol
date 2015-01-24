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

#ifndef RANDUTIL_H_
#define RANDUTIL_H_

		double raisin();

		int initmyrand(int seed);
		unsigned long longinitmyrand(unsigned long *inseed);

		double rand0to1();
		int rand_in_rad(const float rad, float *x, float *y);

		unsigned long randint();
		int * randintarray(const int size,const int Min,const int max);
		int * randboolarray(const int size);

#endif /*RANDUTIL_H_*/

#ifdef __cplusplus
	}
#endif
