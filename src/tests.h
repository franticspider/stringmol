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

#ifndef SRC_TESTS_H_
#define SRC_TESTS_H_


enum failmodes
{
	FAIL_RAND,		//general problem with the random number generator
	FAIL_RAND_INIT	//problem with random initialisation
};





/* FUNCTION DECLARATIONS */

int test_all(int argc, char *argv[]);

/* Premiminary check that a config is sane */
stringPM * test_config_settings( int argc, char *argv[], int return_SM);

/* Test that we can store the state of the rng */
void test_rand_config(int argc, char *argv[]);




#endif /* SRC_TESTS_H_ */
#ifdef __cplusplus
	}
#endif
