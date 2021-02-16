#include "catch.hpp"

#include "../src/mt19937-2.h"
#include "../src/randutil.h"

/*Test the random number generator. It must:
 *
 * 1: give the same sequence from a seed
 * 2: give the current position so that we can re-seed
 */
TEST_CASE("seeding with initmyrand(436) produces 436 as output"){
    int rin = 436;
    int rout = initmyrand(rin);
    REQUIRE(rout == rin);
}


TEST_CASE("reseeding with initmyrand(436) resets the rng and the same sequence is output"){
    int rin = 436;
    int rout = initmyrand(rin);
    REQUIRE(rout == rin);
    
    const int nrand = 10;
    float rnos1[nrand], rnos2[nrand];
    for(int rr=0;rr<nrand;rr++)
        rnos1[rr] = rand0to1();

    rout = initmyrand(rin);
    REQUIRE(rout == rin);


    for(int rr=0;rr<nrand;rr++){
        rnos2[rr] = rand0to1();
        REQUIRE(rnos1[rr] == rnos2[rr]);
    }
}


TEST_CASE("reseeding with initmyrand(-1) resets the rng and a different sequence is output"){
    int rin = -1;
    int rout2,rout1 = initmyrand(rin);
    REQUIRE(rout1 != rin);
    
    const int nrand = 10;
    float rnos1[nrand], rnos2[nrand];
    for(int rr=0;rr<nrand;rr++)
        rnos1[rr] = rand0to1();

    rout2 = initmyrand(rin);
    REQUIRE(rout2 != rin);

    REQUIRE(rout2 != rout1);

    for(int rr=0;rr<nrand;rr++){
        rnos2[rr] = rand0to1();
        REQUIRE(rnos1[rr] != rnos2[rr]);
    }
}
	
/*    
TEST_CASE("resetting Mersenne Twister index"){
    int pos2,pos = 22;
    int rout,rin = 444;
    const int ncalls = 345;
    const int nrand = 20;
    float rnos1[nrand], rnos2[nrand];

    rout = initmyrand(rin);
    //This repeats tests above, but hey- can't hurt can it?
    REQUIRE(rout == rin);
    for(int cc=0;cc<ncalls;cc++){
        rand0to1();
    }
    pos = get_mti();
    for(int cc=0;cc<nrand;cc++){
        rnos1[cc] = rand0to1();
    }
    set_mti(pos);
    pos2 = get_mti();
    
    REQUIRE(pos==pos2);

    for(int rr=0;rr<nrand;rr++){
        rnos2[rr] = rand0to1();
        REQUIRE(rnos1[rr] == rnos2[rr]);
    }    

}*/


TEST_CASE("save and load RNG to/from file"){
    FILE *mtf;
    int runf = 7919; //the 1000th prime
    const int nrand = 20;
    char fn[] = "rng.txt";
    float rnos1[nrand], rnos2[nrand];

    initmyrand(-1);
    //Run the RNG forward a couple of thousand times
    for(int cc=0;cc<runf;cc++){
        rand0to1();
    }

	if((mtf = fopen(fn,"w"))!=NULL){
        //Save the RNG
		print_mt(mtf);
		fclose(mtf);
	}
	else{
		printf("Failed to record RNG state to file %s\n",fn);
        FAIL();
	}

    //Run the RNG forward into array1
    for(int cc=0;cc<nrand;cc++){
        rnos1[cc] = rand0to1();
    }

    //Restore the RNG
    load_mt(fn);

    //Run the RNG forward into array2 and test
    for(int rr=0;rr<nrand;rr++){
        rnos2[rr] = rand0to1();
        REQUIRE(rnos1[rr] == rnos2[rr]);
    }    
    
}
