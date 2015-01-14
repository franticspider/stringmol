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

#ifndef RULES_H_
#define RULES_H_

class rules{


public:

	//variables.
	int nr; 	//The number of rules
	int **rset;	//The set of rules
	float *rval; //The values attached to the rules. (for various purposes)
	float *re;

	//creators and destructors
	rules(char *fn);
	~rules();

	//Finding etc
	int found(const char *side, int label);
	int maxlhs(int label);
	int getrule(const char *side,int l1, int l2);

	//diagnostics:
	void printasint();
	void printrulesfor(const char * side, int rule);

private:

	int findinrange(int row, int min, int max, int label);

};


#endif /*RULES_H_*/
