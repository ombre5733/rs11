//#define CATCH_CONFIG_MAIN
//#include "catch.hpp"

#include "../ReedSolomon.hpp"

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

// Output a table.
ostream& operator<<(ostream& str, const rs11::Galois::Table& a)
{
    for (std::size_t idx = 0; idx < a.size(); ++idx)
        str << setw(3) << idx << "   " << setw(3) << unsigned(a[idx]) << endl;
    return str;
}

ostream& operator<<(ostream& str, rs11::Galois::GF256Value e)
{
    str << e.value();
    return str;
}

template <std::size_t TDeg>
ostream& operator<<(ostream& str, const rs11::Galois::GF256Polynomial<TDeg>& p)
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
    auto pwrTable = rs11::Galois::createPowerTable(0, {});
    cout << "Power table: " << endl;
    cout << pwrTable << endl;

    auto logTable = rs11::Galois::createLogTable(pwrTable, 0, {});
    cout << "Log table: " << endl;
    cout << logTable << endl;
#endif

    rs11::Galois::GF256Value a(2);
    rs11::Galois::GF256Value b(6);
    cout << (a * b).value() << endl;
    cout << (b + b).value() << endl;

    cout << "inv\n";
    cout << (a * a.inverse()).value() << endl;
    cout << (b * b.inverse()).value() << endl;

    for (int e : {255, 256, 511, 512, 513})
        cout << e << "   " << e % 256 + e / 256 << "   " << e % 255 << endl;

    cout << rs11::Galois::GF256Value(1) * rs11::Galois::GF256Value(3) << endl;

    rs11::Galois::GF256Polynomial<1> pa{1, 1};
    rs11::Galois::GF256Polynomial<2> pb{3, 4, 5};

    cout << rs11::Galois::GF256Polynomial<2>() << endl;
    cout << rs11::Galois::GF256Polynomial<2>({1,0,0}) << endl;
    cout << "pa(x) = " << pa << endl;
    cout << "pb(x) = " << pb << endl;
    cout << (pa + pb) << endl;
    cout << (pa * pb) << endl;

    cout << "pa'(2) = " << pa(2, rs11::Galois::derivative) << endl;
    cout << "pb(2) = " << pb(2) << endl;
    cout << "pb'(2) = " << pb(2, rs11::Galois::derivative) << endl;

    auto RSgen = rs11::detail::createReedSolomonGeneratorPolynomial<19>();
    cout << "RSgen:\n";
    cout << RSgen << endl;

    cout << "\n\nReed-Solomon encoding" << endl;
    using Encoder = rs11::ReedSolomonEncoder<255, 251>;
    Encoder rsEnc;
    std::vector<std::uint8_t> data({0x48, 0x65, 0x6c, 0x6c, 0x6f, 0x57, 0x6f, 0x72, 0x6c, 0x64});
    if (1)
    {
        for (auto b : data)
            rsEnc << b;
        rsEnc.finish();
    }
    else
    {
        rsEnc.encode(data.data(), data.size());
    }

    cout << "Parity: ";
    std::vector<std::uint8_t> encData = data;
    for (auto iter = rsEnc.begin(); iter != rsEnc.end(); ++iter)
    {
        cout << std::dec << int(*iter) << " ";
        encData.push_back(*iter);
    }
    cout << endl << endl;

    // Introduce errors.
    auto defData = encData;
    for (int i = 2; i <= 3; ++i)
        defData[i] = 'a';
    cout << "\n\nReed-Solomon decoding" << endl;
    Encoder::decoder_type rsDec;
    for (auto d : defData)
        rsDec << d;
    rsDec.finish();
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


    //rs11::Galois::GF256Polynomial<6> pm{3, 4, 5, 6, 7, 8, 9};
    //cout << pm(2, rs11::Galois::derivative) << endl;


    auto corrData = defData;
    if (rsDec.result() == rs11::DecoderResult::Correctable)
        for (auto e : rsDec.errors())
            corrData[e.index()] = e.correct(defData[e.index()]);




    auto printVector = [] (const char* name, const std::vector<std::uint8_t>& data)
    {
        cout << name;
        cout << " (" << data.size() << "): ";
        for (auto d : data)
            if (d > 32 && d < 128)
                cout << d;
            else
                cout << "\\" << std::hex << std::setw(2) << unsigned(d) << std::dec;
        cout << endl;
    };

    cout << "\n\n";

    cout << "Result: " << int(rsDec.result()) << endl;

    printVector("Original    ", data);
    printVector("Encoded     ", encData);
    printVector("Transmitted ", defData);
    printVector("Corrected   ", corrData);
}
