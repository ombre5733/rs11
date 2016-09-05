#include "catch.hpp"

#include "../galois.h"

TEST_CASE("inversion", "[gf256]")
{
    for (unsigned value = 1; value < 256; ++value)
    {
        GF256Element a(value);
        REQUIRE(a * a.inverse() == GF256Element(1));
        REQUIRE(a.inverse() * a == GF256Element(1));
    }
}
