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
//extern const int maxl;//150;// 512;
//extern const int maxl0; //allow room for a terminating 0


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
// Writing PNGs
#include "lodepng.h"
#include <iostream>
#include "setupSM.h"


#include "webapi_util.h"


//TODO: this is defined in webapi_util.cpp
const int MAXSTR = 5000;




int main(int argc, char*argv[]){

	char 	*output;
	char 	*string1,	//[] = "ADFAERQEDAGADGAERQEWRQERDFSDFASDFERQEW",
			*string2;	//[] = "DSFASDFAFQEREWASDFASDFAREQRQWEFDSFAEQWR";
	char 	*data;

	string1 = (char *) malloc(MAXSTR*sizeof(char));
	string2 = (char *) malloc(MAXSTR*sizeof(char));


#ifdef DEBUG
	printf("DEBUG IS #DEFINED\n");
	fflush(stdout);
#endif

	data = getenv("QUERY_STRING");

#ifdef DEBUG
	printf("%s%c%c\n",
	"Content-Type:text/html;charset=iso-8859-1",13,10);
	printf("<TITLE>Multiplication results</TITLE>\n");
	printf("<H3>Multiplication results</H3>\n");

	if(data==NULL){
		sprintf(string1,"OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B1C$=?1$$BLUBO%%}OYHOB");//"THISISATESTSTRING1ASDF1%%$^=?}ASDF1ASDF1ASDFASDFASDF");
		sprintf(string2,"OOGEOLHHHRLUEUOBBBRBXUUUDYGRHBLROORE$BLUBO^B1C$=?1$$BLUBO%%}OYHOB");//"QWERQW1ERQWER1QWERQ1WERQ1WERQWERQWER");
		printf("<p> Debugging with %s and %s</p>\n",string1,string2);
	}
	else{
		printf("<p>Query string is '%s'</p>\n",data);
	}
#else

	printf("Content-type: text/html\n\n");
#endif

	/*Read the query string into string1 and string2, checking for data errors */
	if(data!=NULL){
		if(sscanf(data,"string1=(%[^)])&string2=(%[^)])",string1,string2)!=2){
				printf("<P>Error! Invalid data. Data must be numeric.</p>");
				sprintf(string1,"ASDFADFASDFASDFASDFASDFASDFASDF");
				sprintf(string2,"QWERQWERQWERQWERQWERQWERQWERQWER");
		}
	}
	else{
		printf("<p> Problem with environment variable QUERY_STRING (== NULL), using debug values instead: </p>\n");
		printf("<p> string1 = \"%s\"</p><p> string2 = \"%s\"</p>",string1,string2);
	}

	/*Now try the bind*/
	printf("%s\n",output = generate_bind_data(string1,string2));

	/*Clean up*/
	free(output);
	free(string1);
	free(string2);

	return 0;
}
