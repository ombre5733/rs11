//#define CATCH_CONFIG_MAIN
//#include "catch.hpp"


#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <utility>

using namespace std;

#define MAX   256


constexpr std::uint8_t reduction_poly = 0x1d;

constexpr
int nextPower(int prev)
{
    return (prev & 0x80) ? (((prev ^ 0x80) << 1) ^ reduction_poly) : (prev << 1);
}





template <typename... TArgs>
using index_sequence_for = make_index_sequence<sizeof...(TArgs)>;


namespace detail
{

class IncrArray
{
    template <std::size_t... TIndices>
    constexpr
    IncrArray doModify(std::size_t idx, std::uint8_t value,
                       index_sequence<TIndices...>) const
    {
        return { (TIndices == idx ? value : m_data[TIndices])... };
    }

public:
    constexpr
    IncrArray modified(std::size_t idx, int value) const
    {
        return doModify(idx, value, make_index_sequence<MAX>());
    }

    constexpr
    std::size_t size() const
    {
        return MAX;
    }

    constexpr
    std::uint8_t operator[](std::size_t idx) const
    {
        return m_data[idx];
    }

    std::uint8_t m_data[MAX];
};

constexpr
auto createPowerTable(std::size_t idx, const IncrArray& temp) -> IncrArray
{
    return idx == MAX ? temp
                      : createPowerTable(
                            idx + 1,
                            temp.modified(idx, idx == 0 ? 1 : nextPower(temp[idx-1])));
}

constexpr
IncrArray createLogTable(const IncrArray& powers, std::size_t idx, const IncrArray& temp)
{
    return idx == MAX ? temp
                      : createLogTable(
                            powers, idx + 1,
                            temp.modified(powers[idx], powers[idx] == 1 ? 0 : idx));
}

} // namespace detail





//! An element from GF(2^8).
//!
//! The GF256Element is an element from the Galois field GF(2^8). Every
//! element from GF(2^8) is a polynomial with degree strictly less than 8 of
//! the form \f$ p(x) = \sum_{i=0}^7 a_i x^i \f$, with \f$ a_i \in GF(2) \f$.
//! Since \f$ a_i \in {0,1} \f$, we represent an element of GF(2^8) as
//! 8-bit unsigned integer. For example, the polynomial \f$ x^5 + x^2 + x^1 \f$
//! has the binary representation <tt>0b00100110</tt>.
//!
//! Let \p alpha denote one primitive element of GF(2^8). The powers
//! \f$ \alpha^i, i = 0, ..., 254 \f$ generate all non-zero elements of
//! GF(2^8).
// \alpha is defined as root of an irreducible polynomial with degree 8, e.g.
// p(x) = x^8 + x^1 + 1
// template <unsigned TReductionPoly>
class GF256Element
{
public:
    //! Creates the zero element.
    constexpr
    GF256Element() noexcept
        : m_value(0)
    {
    }

    constexpr
    GF256Element(std::uint8_t value) noexcept
        : m_value(value)
    {
    }

    //! Computes the sum of two elements from GF(2^8).
    //!
    //! Adds \p b to this element and returns the sum.
    constexpr
    GF256Element operator+(GF256Element b) const noexcept
    {
        // Addition in GF(2) is the XOR function.
        return m_value ^ b.m_value;
    }

    //! Computes the difference between two elements from GF(2^8).
    //!
    //! Subtracts \p b from this element and returns the difference.
    constexpr
    GF256Element operator-(GF256Element b) const noexcept
    {
        // Subtraction in GF(2) is just an addition.
        return m_value ^ b.m_value;
    }

    //! Computes the product of two elements from GF(2^8).
    //!
    //! Multiplies this element with \p b and returns the product.
    constexpr
    GF256Element operator*(GF256Element b) const noexcept
    {
        // Exploit the relation a * b = g^{log_g(a) + log_g(b)}
        return (m_value == 0 || b.m_value == 0)
               ? 0
               : m_antiLogTable[mod(m_logTable[m_value] + m_logTable[b.m_value])];
    }

    GF256Element inverse() const noexcept
    {
        return m_antiLogTable[255 - m_logTable[m_value]];
    }

    // TODO: Remove
    constexpr
    unsigned value() const noexcept
    {
        return m_value;
    }

    constexpr
    explicit operator std::uint8_t() const noexcept
    {
        return m_value;
    }

    explicit
    operator bool() const noexcept
    {
        return m_value != 0;
    }

    //! Returns 2 raised to the power \p p.
    //!
    //! Returns the GF(2^8) element <tt>2^p</tt>.
    static constexpr
    GF256Element pow2(unsigned p)
    {
        return m_antiLogTable[p % 255];
    }

private:
    std::uint8_t m_value;

    static constexpr
    unsigned mod(unsigned x) noexcept
    {
        // Compute x % |g| = x % 255. |g| is the order of the generator
        // element, i.e. the number of non-zero elements of the field.
        return x < 255 ? x : x - 255;
    }

    // The entry at index i equals \alpha^i modulo the reduction polynomial.
    static constexpr detail::IncrArray m_antiLogTable = detail::createPowerTable(0, {});
    static constexpr detail::IncrArray m_logTable = detail::createLogTable(GF256Element::m_antiLogTable, 0, {});
};

constexpr detail::IncrArray GF256Element::m_antiLogTable;
constexpr detail::IncrArray GF256Element::m_logTable;


template <typename T>
constexpr
T const& cmax(T const& a, T const& b)
{
    return a > b ? a : b;
}

//! A polynomial with coefficients in GF(2^8).
//!
//! The GF256Polynomial is a polynomial whose coefficients are elements
//! of the Galois field GF(2^8).
template <std::size_t TDegree>
class GF256Polynomial
{
    template <std::size_t TDeg, std::size_t... TIndices>
    constexpr
    GF256Polynomial<sizeof...(TIndices) - 1>
    doAdd(const GF256Polynomial<TDeg>& b,
          index_sequence<TIndices...>) const noexcept
    {
        return { (coeff(TIndices) + b.coeff(TIndices))... };
    }

    template <std::size_t TDeg>
    constexpr
    GF256Element doConvPart(const GF256Polynomial<TDeg>& b,
                            std::size_t degA, std::size_t degB) const noexcept
    {
        return degA == 0 ? coeff(degA) * b.coeff(degB)
                         : coeff(degA) * b.coeff(degB)
                           + doConvPart(b, degA - 1, degB + 1);
    }

    template <std::size_t TDeg, std::size_t... TIndices>
    constexpr
    GF256Polynomial<sizeof...(TIndices) - 1>
    doMul(const GF256Polynomial<TDeg>& b,
          index_sequence<TIndices...>) const noexcept
    {
        return { doConvPart(b, TIndices, 0)... };
    }

public:
//    constexpr
//    GF256Polynomial() noexcept
//        : _m_coefficients{}
//    {
//    }

    template <std::size_t TDeg>
    constexpr
    auto operator+(const GF256Polynomial<TDeg>& b) const noexcept
        -> GF256Polynomial<cmax(TDegree ,TDeg)>
    {
        return doAdd(b, make_index_sequence<cmax(TDegree, TDeg) + 1>());
    }

    template <std::size_t TDeg>
    constexpr
    auto operator-(const GF256Polynomial<TDeg>& b) const noexcept
        -> GF256Polynomial<cmax(TDegree, TDeg)>
    {
        // Subtraction in GF(2^8) equals addition.
        return doAdd(b, make_index_sequence<cmax(TDegree, TDeg) + 1>());
    }

    //! Returns the product of two polynomials.
    template <std::size_t TDeg>
    constexpr
    auto operator*(const GF256Polynomial<TDeg>& b) const noexcept
        -> GF256Polynomial<TDegree + TDeg>
    {
        return doMul(b, make_index_sequence<TDegree + TDeg + 1>());
    }

    //! Evaluates the polynomial.
    GF256Element operator()(GF256Element x) const noexcept
    {
        GF256Element sum{_m_coefficients[TDegree]};
        for (std::size_t deg = TDegree; deg > 0; --deg)
            sum = sum * x + _m_coefficients[deg - 1];
        return sum;
    }

    //! Evaluates the polynomial of a certain degree.
    GF256Element operator()(GF256Element x, std::size_t degree) const noexcept
    {
        GF256Element sum{_m_coefficients[degree]};
        for (std::size_t deg = degree; deg > 0; --deg)
            sum = sum * x + _m_coefficients[deg - 1];
        return sum;
    }

    //! Returns the coefficient of a monomial.
    constexpr
    GF256Element coeff(std::size_t degree) const noexcept
    {
        return degree <= TDegree ? _m_coefficients[degree] : GF256Element();
    }

    //! Returns the coefficient of a monomial.
    GF256Element& operator[](std::size_t degree) noexcept
    {
        return _m_coefficients[degree];
    }

    //! Multiplies the polynomial inplace with the polynomial \f$ p(x) = x \f$.
    void timesX() noexcept
    {
        // Shift the coefficients [a, b, c, d] -> [d, a, b, c].
        std::rotate(&_m_coefficients[0], &_m_coefficients[TDegree], &_m_coefficients[TDegree + 1]);
        _m_coefficients[0] = 0;
    }

    //! Returns the degree of the polynomial.
    std::size_t degree() const noexcept
    {
        for (std::size_t deg = TDegree; deg > 0; --deg)
            if (_m_coefficients[deg])
                return deg;
        return 0;
    }



    // The coefficients of the polynomial.
    GF256Element _m_coefficients[TDegree + 1];
};


namespace detail
{

//! Returns the polynomial \f$ p(x) = prod_{i=0}^N (x - \alpha^i) \f$ with
//! \f$ alpha = 2 \f$.
template <std::size_t N>
constexpr
GF256Polynomial<N> createReedSolomonGeneratorPolynomial()
{
    return createReedSolomonGeneratorPolynomial<N - 1>()
            * GF256Polynomial<1>({1, GF256Element::pow2(N - 1)});
}

//! Returns the polynomial \f$ p(x) = 1 \f$.
template <>
constexpr
GF256Polynomial<0> createReedSolomonGeneratorPolynomial()
{
    return GF256Polynomial<0>({1});
}

} // namespace detail


//! A Reed-Solomon encoder.
//!
//! The Reed-Solomon code has a total number of \p M symbols (the block size),
//! \p N payload symbols and <tt>(M - N)</tt> parity symbols. The RS code makes
//! use of 8-bit symbols.
//!
//! The class holds internal state needed for encoding messages.
template <std::size_t M, std::size_t N>
class ReedSolomonEncoder
{
    static_assert(M <= 255, "Only 255 symbols are supported");
    static_assert(N < M, "The number of payload symbols must be smaller than the total number of symbols");

    static constexpr std::size_t numParitySyms = M - N;

public:
    constexpr
    ReedSolomonEncoder() noexcept
        : m_polynomial{detail::createReedSolomonGeneratorPolynomial<numParitySyms>()},
          m_scratchIter{&m_scratch[0]},
          m_finished{false}
    {
    }

    ReedSolomonEncoder(const ReedSolomonEncoder&) = delete;
    ReedSolomonEncoder& operator=(const ReedSolomonEncoder&) = delete;

    //! Encodes a message.
    //!
    //! Encodes the given \p message consisting of \p length bytes.
    void encode(std::uint8_t* message, std::size_t length) noexcept
    {
        encodePart(message, length);
        finish();
    }

    void encodePart(std::uint8_t* message, std::size_t length) noexcept
    {
        // TODO: If length exceeds the payload size, throw an exception.

        // Fill up the scratch space except the last slot.
        while (m_scratchIter != &m_scratch[numParitySyms] && length)
        {
            *m_scratchIter++ = *message++;
            --length;
        }

        // Perform a long division for every new element.
        while (length)
        {
            *m_scratchIter = *message++;
            --length;
            polyLongDiv();
        }
    }

    //! Finishes encoding.
    void finish() noexcept
    {
        while (m_scratchIter < &m_scratch[numParitySyms + 1])
            *m_scratchIter++ = 0;
        for (std::size_t count = 0; count < numParitySyms; ++count)
            polyLongDiv();
    }

    ReedSolomonEncoder& operator<<(std::uint8_t datum) noexcept
    {
        *m_scratchIter = datum;
        if (m_scratchIter == &m_scratch[numParitySyms])
            polyLongDiv();
        else
            ++m_scratchIter;
        return *this;
    }

    std::uint8_t* begin() noexcept
    {
        static_assert(sizeof(GF256Element) == sizeof(std::uint8_t), "");
        return reinterpret_cast<std::uint8_t*>(&m_scratch[0]);
    }

    std::uint8_t* end() noexcept
    {
        static_assert(sizeof(GF256Element) == sizeof(std::uint8_t), "");
        // Spare the last element.
        return reinterpret_cast<std::uint8_t*>(&m_scratch[numParitySyms]);
    }

    //! Returns the parity symbols.
    void paritySymbols() const noexcept;

    //! Returns the number of parity symbols.
    std::size_t numParitySymbols() const noexcept
    {
        return numParitySyms;
    }

private:
    // Performs a polynomial long division.
    //
    // Divides the polynomial in m_scratch by the generator.
    void polyLongDiv() noexcept
    {
        // We know that the highest order coefficient of the ERC polynomial
        // is one. So we can skip the division.
        GF256Element quotient = m_scratch[0];
        // When subtracting the intermediate result, the first term is
        // zeroed. Shift the lower order terms into place.
        for (size_t idx = 1; idx <= numParitySyms; ++idx)
            m_scratch[idx - 1] = m_scratch[idx] - m_polynomial.coeff(idx) * quotient;
    }

    //! The polynomial for Reed-Solomon encoding/decoding.
    GF256Polynomial<numParitySyms> m_polynomial;

    //! A scratch space to perform the polynomial long division.
    GF256Element m_scratch[numParitySyms + 1];
    //! An iterator over the scratch space.
    GF256Element* m_scratchIter;

    bool m_finished;
};



class Error
{
public:
    std::size_t location() const noexcept
    {
        return m_location;
    }

    std::uint8_t correct(std::uint8_t datum) const noexcept
    {
        return 0; // TODO
    }

private:
    std::size_t m_location;
    GF256Element m_magnitude;
};

enum class DecoderResult
{
    Correct,
    Correctable,
    NonCorrectable
};

template <std::size_t M, std::size_t N>
class ReedSolomonDecoder
{
    static_assert(M <= 255, "Only 255 symbols are supported");
    static_assert(N < M, "The number of payload symbols must be smaller than the total number of symbols");

    static constexpr std::size_t numParitySyms = M - N;

public:
    ReedSolomonDecoder()
        : m_size(0)
    {
    }

    ReedSolomonDecoder(const ReedSolomonDecoder&) = delete;
    ReedSolomonDecoder& operator=(const ReedSolomonDecoder&) = delete;

    ReedSolomonDecoder& operator<<(std::uint8_t datum) noexcept
    {
        ++m_size;

        // To compute the syndromes, the encoded message must be evaluated
        // at the roots of the generator polynomial, i.e. at
        // \alpha^i, i = 0 ... N, with \alpha = 2. If we had no errors in
        // the encoded message, the syndromes would be all zero.
        for (size_t idx = 0; idx < numParitySyms; ++idx)
        {
            // TODO: Speed up with
            // if m_syndrome[idx] == 0: m_syndrome[idx] = datum;
            m_syndromes[idx] = m_syndromes[idx] * GF256Element::pow2(idx) + GF256Element(datum);
        }

        return *this;
    }

    void bm()
    {
        // Berlekamp-Massey algorithm
        for (unsigned idx = 0; idx <= numParitySyms; ++idx)
            m_errorLocator[idx] = m_B[idx] = 0;
        m_errorLocator[0] = m_B[0] = GF256Element(1);

        std::size_t numErrors = 0;
        for (std::size_t iteration = 0; iteration < numParitySyms; ++iteration)
        {
            GF256Element discrepancy;
            for (std::size_t idx = 0; idx <= iteration; ++idx)
                discrepancy = discrepancy + m_errorLocator[idx] * m_syndromes[iteration - idx]; // TODO: +=

            cout << "iter = " << iteration << "   delta = " << discrepancy.value() << "\n";

            if (discrepancy.value() == 0) // TODO: operator==
            {
                // No discrepancy.
                // B(x) <- x * B(x)
                m_B.timesX();
            }
            else
            {
                // Lambda(x) <- Lambda(x) - Delta * x * B(x)
                m_temp[0] = m_errorLocator[0];
                for (size_t idx = 0; idx < numParitySyms; ++idx)
                   m_temp[idx + 1] = m_errorLocator[idx + 1] - discrepancy * m_B[idx];

                if (2 * numErrors <= iteration)
                {
                    numErrors = iteration + 1 - numErrors;
                    // B(x) <- Delta^{-1} Lambda(x)
                    discrepancy = discrepancy.inverse();
                    for (size_t idx = 0; idx <= numParitySyms; ++idx)
                        m_B[idx] = discrepancy * m_errorLocator[idx];
                }
                else
                {
                    // B(x) <- x * B(x)
                    m_B.timesX();
                }

                m_errorLocator = m_temp;
            }

            cout << "   m_errorLocator ";
            for (auto coeff : m_errorLocator._m_coefficients)
                cout << " " << coeff.value();
            cout << "\n";

            cout << "   m_errorLocator " << m_errorLocator << "\n";
        }


        // Find the roots of the error locator polynomial. The roots are the
        // inverse of the error locations.
        unsigned numRoots = 0;
        auto errorLocatorDegree = m_errorLocator.degree();
        for (unsigned idx = 1; idx <= 255; ++idx)
        {
            if (!m_errorLocator(GF256Element::pow2(idx), errorLocatorDegree))
            {
                cout << "   - Root: i = " << idx
                     << "   (" << 255 - idx << ")"
                     << "   " << GF256Element::pow2(idx).value()
                     << "   " << GF256Element::pow2(idx).inverse().value() << endl;
                cout << " DBG " << m_errorLocator(0).value() << endl;

                auto invIndex = 255 - idx;
                if (m_size - 1 < invIndex)
                    return; // TODO: unable to correct the error

                cout << "Error at " << m_size - 1 - invIndex << "\n";
                if (++numRoots == errorLocatorDegree)
                    break;
            }
        }

        if (numRoots != errorLocatorDegree)
            return; // TODO: unable to correct the error
    }

    bool hasError() const noexcept
    {
        for (auto syn : m_syndromes)
            if (syn.value() != 0)
                return true;
        return false;
    }

public:
    //! The syndromes.
    GF256Element m_syndromes[numParitySyms];

    //! The error locator polynomial.
    GF256Polynomial<numParitySyms> m_errorLocator;
    //! A scratch polynomial needed to update the error locator polynomial.
    GF256Polynomial<numParitySyms> m_temp;

    GF256Polynomial<numParitySyms> m_B;

    std::size_t m_size;
};


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
    cout << pa << endl;
    cout << pb << endl;
    cout << (pa + pb) << endl;
    cout << (pa * pb) << endl;

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
    for (int i = 1; i <= 1; ++i)
        encData[i] = 'a';
    cout << "\n\nReed-Solomon decoding" << endl;
    ReedSolomonDecoder<255, 251> rsDec;
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
}
