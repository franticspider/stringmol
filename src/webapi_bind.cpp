/* Copyright (C) 2009-2015 Simon Hickinbotham, Adam Nellis, Ed Clark    */
/* When you use this, send an email to: sjh436@gmail.com                */
/* with an appropriate reference to your work.                          */

/* This file is part of STRINGMOL					*/

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


/* TODO: Description */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


//utils
//extern "C" {
#include "randutil.h"
#include "params.h"
#include "hsort.h"
//}

//stringmol
#include "alignment.h"

//metabolism
#include "rules.h"
#include "agents_base.h"
#include "SMspp.h"
#include "stringPM.h"

//signal
#include "signalSM.h"

//setup
#include "setupSM.h"

#define DEBUG

int main(int argc, char*argv[]){

#ifdef DEBUG
	printf ("Hello web!\n");
#endif
	//Create the container
	SMspp		SP;
	stringPM 	A(&SP);

	//Create two molecules from strings


}
