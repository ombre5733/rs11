#include "reed_solomon.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>


using namespace std;


#if 0
GF256Polynomial createReedSolomonPoly(std::size_t terms,
                                      const GF256Polynomial& previous)
{
    return terms == 0 ? previous : createReedSolomonPoly(terms - 1, previous * GF256Polynomial(1, GF256Element::pow(terms)));
}
#endif

constexpr GF256Polynomial<4> POLY_255_251 = { { 1, 30, 216, 231, 116 } };

#include "testvectors.h"

// Usage model for AT45DB FLASH (512/528 byte page)
// Shortened (zero padded) RS<171,175> build from RS<255,251>
// 512 bytes spitted in 3 segments:
// Option A:    171 171 170  +  4 4 4 parity = 524 bytes -> 4 bytes unused (wear leveling or possibly CRC32 for fast correction need check)
// Option B:    172 172 168  +  4 4 4 parity = 524 bytes -> 4 bytes unused (wear leveling or possibly CRC32 for fast correction need check)
// Only benefit might by that 172/4 is integer
// CRC32 might also ensure error detection for more than 2 symbols in error ???

int main()
{
	using namespace std;
    int i = 1;
    bool allpass = true;
    constexpr GF256Polynomial<4> ZERO;
    for (auto &r : testVectors)
    {
        ReedSolomon<255, 251> rs;
        rs.encode(&r[0], 251);
        auto S = rs.syndrome(r);
        auto LO = rs.RiBM(S);

        bool equal = memcmp(rs.begin(), &r[251], 4) == 0;
        allpass &= equal && S == ZERO;
        cout << "test Vector " << i++ << " " << (equal ? "pass" : "fail") << endl;
    }
    cout << "---------------" << endl << "Test " << (allpass ? "passed" : "failed") << endl;
}
