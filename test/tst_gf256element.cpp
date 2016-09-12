#include "catch.hpp"

#include "../GaloisField.hpp"

using namespace rs11::Galois;

TEST_CASE("addition", "[gf256]")
{
    for (unsigned value = 0; value < 256; ++value)
    {
        GF256Value a(value);
        REQUIRE(a + a == GF256Value(0));
    }
}

TEST_CASE("subtraction", "[gf256]")
{
    for (unsigned value = 0; value < 256; ++value)
    {
        GF256Value a(value);
        REQUIRE(a - a == GF256Value(0));
    }
}

TEST_CASE("inversion", "[gf256]")
{
    for (unsigned value = 1; value < 256; ++value)
    {
        GF256Value a(value);
        REQUIRE(a * a.inverse() == GF256Value(1));
        REQUIRE(a.inverse() * a == GF256Value(1));
    }
}
