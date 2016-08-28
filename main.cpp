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
        bool equal = memcmp(rs.begin(), &r[251], 4) == 0;
        allpass &= equal && S == ZERO;
        cout << "test Vector " << i++ << " " << (equal ? "pass" : "fail") << endl;
    }
    cout << "---------------" << endl << "Test " << (allpass ? "passed" : "failed") << endl;
}
