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

}




/*
TEST_CASE("Triangle Constructor Throws Exception When One Side Is Negative", "[Triangle]") {
	REQUIRE_THROWS(Triangle(-10.0f, 10.0f, 10.0f));
}

TEST_CASE("Triangle Constructor Throws Exception When Not A Valid Triangle", "[Triangle]") {
	REQUIRE_THROWS(Triangle(10.0f, 10.0f, 25.0f));
}

TEST_CASE("Triangle Constructor Works When Valid Data Is Used", "[Triangle]") {
	REQUIRE_NOTHROW(Triangle(10.0f, 15.0f, 20.0f));
}

TEST_CASE("Triangle Area Is Computed Correctly", "[Triangle]") {
	Triangle Triangle(10.0f, 15.0f, 20.0f);
	float RoundedArea = roundf(Triangle.GetArea() * 100) / 100;

	REQUIRE(RoundedArea == 72.62f);
}

TEST_CASE("Triangle Perimiter Is Computed Correctly", "[Triangle]") {
	Triangle Triangle(10.0f, 15.0f, 20.0f);
	float Perimiter = Triangle.GetPerimiter();

	REQUIRE(Perimiter == 45.0f);
}*/
