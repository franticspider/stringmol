#include "catch.hpp"
#include "../Source/Headers/Square.h"

TEST_CASE("Square Constructor Throws Exception When Side Length Is Zero", "[Square]") {
	REQUIRE_THROWS(Square(0.0f));
}

TEST_CASE("Square Constructor Throws Exception When Side Length Is Negative", "[Square]") {
	REQUIRE_THROWS(Square(-1.0f));
}

TEST_CASE("Square Constructor Works When Valid Data Is Used", "[Square]") {
	REQUIRE_NOTHROW(Square(10.0f));
}

TEST_CASE("Square GetPerimiter() Works Correctly", "[Square]") {
	Square Square(10.0f);

	REQUIRE(Square.GetPerimiter() == 40.0f);
}

TEST_CASE("Square GetArea() Works Correctly", "[Square]") {
	Square Square(10.0f);

	REQUIRE(Square.GetArea() == 100.0f);
}