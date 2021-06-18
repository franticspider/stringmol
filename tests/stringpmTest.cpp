//#include <stdlib.h>
//#include <stdio.h>
#include <string.h>

#include "catch.hpp"

/*NB: The following object files were needed to compile for stringPM: 
../release/mt19937-2.o 
../release/randutil.o 
../release/stringPM.o 
../release/agents_base.o 
../release/SMspp.o 
../release/rules.o 
../release/alignment.o 
../release/params.o 
../release/hsort.o 
../release/stringmanip.o 
../release/memoryutil.o 
../release/instructions.o
*/

//#include "randutil.h"
//#include "../src/params.h"      // reading and loading parameters
//ls #include "../src/hsort.h"       // hsort - ?
//}

//stringmol
#include "../src/alignment.h"   // the align object

//metabolism
#include "../src/rules.h"       // the rules object
#include "../src/SMspp.h"       // this defines s_ag among other things
#include "../src/agents_base.h" // parent class of stringPM
#include "../src/stringPM.h"


TEST_CASE("s_ag lifecycle execules correctly (without firing any opcodes)"){
	SMspp		SP;
	stringPM	spm(&SP);
    char        sequence[] = "BLUBO";
    s_ag        *pag;

    pag = NULL;
    REQUIRE(pag == NULL);
	pag = spm.make_ag('A');

    REQUIRE(pag != NULL);

	pag->S =(char *) malloc(spm.maxl0*sizeof(char));
	memset(pag->S,0,spm.maxl0*sizeof(char));
    sprintf(pag->S,sequence,strlen(sequence));

    spm.free_ag(&pag);
    REQUIRE(pag == NULL);
}

