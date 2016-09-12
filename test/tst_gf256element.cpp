#include "catch.hpp"

#include "../GaloisField.hpp"

using namespace rs11::Galois;

TEST_CASE("inversion", "[gf256]")
{
    for (unsigned value = 1; value < 256; ++value)
    {
        GF256Value a(value);
        REQUIRE(a * a.inverse() == GF256Value(1));
        REQUIRE(a.inverse() * a == GF256Value(1));
    }
}
