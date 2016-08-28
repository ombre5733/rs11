#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include "reed_solomon.h"

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

int main()
{
	using namespace std;
    int i = 1;
    bool allpass = true;
    for (auto &r : testVectors)
    {
        ReedSolomon<255, 251> rs;
        rs.encode(&r[0], 251);
        bool equal = memcmp(rs.begin(), &r[251], 4) == 0;
        allpass &= equal;
        cout << "test Vector " << i++ << " " << (equal ? "pass" : "fail") << endl;
        auto S = rs.syndrome(r);
    }
    cout << "---------------" << endl << "Test " << (allpass ? "passed" : "failed") << endl;
}
