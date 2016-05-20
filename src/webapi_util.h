/* Copyright (C) 2009-2015 Simon Hickinbotham                           */
/* When you use this, send an email to: sjh436@gmail.com                */
/* with an appropriate reference to your work.                          */

/* This file is part of STRINGMOL										*/

/* STRINGMOL is free software: you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */

/* STRINGMOL is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */

/* You should have received a copy of the GNU General Public License    */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#ifdef __cplusplus
	extern "C" {
#endif

#ifndef SRC_WEBAPI_UTIL_H_
#define SRC_WEBABI_UTIL_H_


/* FUNCTION DECLARATIONS */

s_ag * make_mol(stringPM *A, char* string);
char * generate_bind_data(char * string1, char * string2);




#endif /* SRC_WEBAPI_UTIL_H_ */
#ifdef __cplusplus
	}
#endif
