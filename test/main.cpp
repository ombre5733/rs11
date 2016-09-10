//#define CATCH_CONFIG_MAIN
//#include "catch.hpp"

#include "reed_solomon.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <utility>



using namespace std;

// ----=====================================================================----
// ----=====================================================================----

#include <iostream>
#include <iomanip>
#include <vector>


// Output a table.
ostream& operator<<(ostream& str, const detail::IncrArray& a)
{
    for (std::size_t idx = 0; idx < a.size(); ++idx)
        str << setw(3) << idx << "   " << setw(3) << unsigned(a[idx]) << endl;
    return str;
}

ostream& operator<<(ostream& str, GF256Element e)
{
    str << e.value();
    return str;
}

template <std::size_t TDeg>
ostream& operator<<(ostream& str, const GF256Polynomial<TDeg>& p)
{
    for (std::size_t i = 0; i <= TDeg; ++i)
    {
        if (i)
            str << " + ";
        str << p.coeff(i) << " x^" << i;
    }
    str << "   (" << p.degree() << ")";
    return str;
}




int main()
{
    using namespace std;

#if 1
    auto pwrTable = detail::createPowerTable(0, {});
    cout << "Power table: " << endl;
    cout << pwrTable << endl;

    auto logTable = detail::createLogTable(pwrTable, 0, {});
    cout << "Log table: " << endl;
    cout << logTable << endl;
#endif

    GF256Element a(2);
    GF256Element b(6);
    cout << (a * b).value() << endl;
    cout << (b + b).value() << endl;

    cout << "inv\n";
    cout << (a * a.inverse()).value() << endl;
    cout << (b * b.inverse()).value() << endl;

    for (int e : {255, 256, 511, 512, 513})
        cout << e << "   " << e % 256 + e / 256 << "   " << e % 255 << endl;

    cout << GF256Element(1) * GF256Element(3) << endl;

    GF256Polynomial<1> pa{1, 1};
    GF256Polynomial<2> pb{3, 4, 5};

    cout << GF256Polynomial<2>() << endl;
    cout << GF256Polynomial<2>({1,0,0}) << endl;
    cout << "pa(x) = " << pa << endl;
    cout << "pb(x) = " << pb << endl;
    cout << (pa + pb) << endl;
    cout << (pa * pb) << endl;

    cout << "pa'(2) = " << pa(2, derivative) << endl;
    cout << "pb(2) = " << pb(2) << endl;
    cout << "pb'(2) = " << pb(2, derivative) << endl;

    auto RSgen = detail::createReedSolomonGeneratorPolynomial<19>();
    cout << "RSgen:\n";
    cout << RSgen << endl;

    cout << "\n\nReed-Solomon encoding" << endl;
    ReedSolomonEncoder<255, 251> rsEnc;
    std::vector<std::uint8_t> data({0x48, 0x65, 0x6c, 0x6c, 0x6f, 0x57, 0x6f, 0x72, 0x6c, 0x64});
    for (auto b : data)
        rsEnc << b;
    rsEnc.finish();

    cout << "Parity: ";
    std::vector<std::uint8_t> encData = data;
    for (auto iter = rsEnc.begin(); iter != rsEnc.end(); ++iter)
    {
        cout << std::dec << int(*iter) << " ";
        encData.push_back(*iter);
    }
    cout << endl << endl;

    // Introduce errors.
    for (int i = 1; i <= 4; ++i)
        encData[i] = 'a';
    cout << "\n\nReed-Solomon decoding" << endl;
    ReedSolomonDecoder<255, 245> rsDec;
    for (auto d : encData)
        rsDec << d;
    rsDec.bm();
    cout << "Syndrome: ";
    for (auto iter : rsDec.m_syndromes)
    {
        cout << std::dec << iter.value() << " ";
    }
    cout << endl;

    cout << "Error locator: ";
    for (auto iter : rsDec.m_errorLocator._m_coefficients)
    {
        cout << std::dec << iter.value() << " ";
    }
    cout << endl << endl;


    cout << "Error evaluator: ";
    for (auto iter : rsDec.m_errorEvaluator._m_coefficients)
    {
        cout << std::dec << iter.value() << " ";
    }
    cout << endl << endl;


    GF256Polynomial<6> pm{3, 4, 5, 6, 7, 8, 9};
    cout << pm(2, derivative) << endl;
}
