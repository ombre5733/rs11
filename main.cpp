#include <array>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <utility>

using namespace std;

#define MAX   256

constexpr
int nextPower(int prev)
{
    return (prev & 0x80) ? (((prev ^ 0x80) << 1) ^ 0x1d) : (prev << 1);
}





#ifndef _MSC_VER
template <typename T, T... TIndices>
struct integer_sequence
{
    typedef T value_type;

    static constexpr std::size_t size()
    {
        return sizeof...(TIndices);
    }
};

namespace weos_detail
{

template <typename T, typename TSequence1, typename TSequence2>
struct AppendIntegerSequence;

template <typename T, T... TIndices1, T... TIndices2>
struct AppendIntegerSequence<T, integer_sequence<T, TIndices1...>, integer_sequence<T, TIndices2...>>
{
    using type = integer_sequence<T, TIndices1..., (sizeof...(TIndices1) + TIndices2)...>;
};

template <typename T, std::size_t N>
struct MakeIntegerSequence
{
    using type = typename AppendIntegerSequence<T,
                                                typename MakeIntegerSequence<T, N/2>::type,
                                                typename MakeIntegerSequence<T, N - N/2>::type>::type;
};

template <typename T>
struct MakeIntegerSequence<T, 0>
{
    using type = integer_sequence<T>;
};

template <typename T>
struct MakeIntegerSequence<T, 1>
{
    using type = integer_sequence<T, 0>;
};

template <typename T, T N>
struct MakeSafeIntegerSequence
{
    // Break compilation if N is negative. Otherwise, the compiler tries to
    // create lots and lots of template instances.
    static_assert(N >= 0, "The sequence length must be positive.");
    using type = typename MakeIntegerSequence<T, N>::type;
};

} // namespace weos_detail

template <std::size_t... TValues>
using index_sequence = integer_sequence<std::size_t, TValues...>;

//! Creates the integer sequence 0, 1, 2, ... N - 1 of type T.
template <typename T, T N>
using make_integer_sequence = typename weos_detail::MakeSafeIntegerSequence<T, N>::type;

//! Creates the integer sequence 0, 1, 2, ... N - 1 with type std::size_t.
template <std::size_t N>
using make_index_sequence = make_integer_sequence<std::size_t, N>;
#endif

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

ostream& operator<< (ostream& str, const IncrArray& a)
{
    for (std::size_t idx = 0; idx < a.size(); ++idx)
        str << setw(3) << idx << "   " << setw(3) << unsigned(a[idx]) << endl;
    return str;
}

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
//! has the binary representation 0b00100110.
//!
//! Let \alpha denote one primitive element of GF(2^8). The powers
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
        return m_value ^ b.m_value;
    }

    constexpr
    GF256Element operator*(GF256Element b) const noexcept
    {
        // Exploit the relation a * b = g^{log_g(a) + log_g(b)}
        return (m_value == 0 || b.m_value == 0)
               ? 0
               : m_antiLogTable[mod(m_logTable[m_value] + m_logTable[b.m_value])];
    }

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

//    //! Returns the 2 raised to the power \p p.
//    //!
//    //! Returns the GF(2^8) element 2^p.
//    static GF256Element exp(unsigned p)
//    {
//        return m_antiLogTable[trueMod(p)];
//    }

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


template<typename T> constexpr
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
    constexpr GF256Polynomial<cmax(TDegree ,TDeg)> operator+(const GF256Polynomial<TDeg>& b) const noexcept
    {
        return doAdd(b, make_index_sequence<cmax(TDegree, TDeg) + 1>());
    }

    template <std::size_t TDeg>
    constexpr auto operator-(const GF256Polynomial<TDeg>& b) const noexcept -> GF256Polynomial<cmax(TDegree, TDeg)>
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

    //! Returns the coefficient of a monomial.
    constexpr
    GF256Element coeff(std::size_t degree) const noexcept
    {
        return degree <= TDegree ? _m_coefficients[degree] : GF256Element();
    }

    // The coefficients of the polynomial.
    GF256Element _m_coefficients[TDegree + 1];
};

template <std::size_t N>
constexpr
GF256Polynomial<N> prod()
{
    return GF256Polynomial<1>({1, N}) * prod<N - 1>();
}

template <>
constexpr
GF256Polynomial<0> prod<0>()
{
    return GF256Polynomial<0>({1});
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
    return str;
}


#if 0
GF256Polynomial createReedSolomonPoly(std::size_t terms,
                                      const GF256Polynomial& previous)
{
    return terms == 0 ? previous : createReedSolomonPoly(terms - 1, previous * GF256Polynomial(1, GF256Element::pow(terms)));
}
#endif


// Polynomial
// \prod_{i=0}^{?} (x - \alpha^i), \alpha \in GF(2^8)
template <std::size_t M, std::size_t N>
class ReedSolomon
{
    static_assert(N < M, "");
    static constexpr std::size_t numParity = M - N;

public:
    constexpr
    ReedSolomon() noexcept
        : m_scratchIter(&m_scratch[0]),
          m_finished(false)
    {
    }

    void finish() noexcept
    {
        while (m_scratchIter < &m_scratch[numParity + 1])
            *m_scratchIter++ = 0;
        for (std::size_t count = 0; count < numParity; ++count)
            polyLongDiv();
    }

    void encode(std::uint8_t* message, std::size_t length) noexcept
    {
        while (m_scratchIter != &m_scratch[numParity] && length)
        {
            *m_scratchIter++ = *message++;
            --length;
        }

        while (length)
        {
            *m_scratchIter = *message++;
            --length;
            polyLongDiv();
        }

        finish();
    }

    ReedSolomon& operator<< (std::uint8_t datum) noexcept
    {
        *m_scratchIter = datum;
        if (m_scratchIter == &m_scratch[numParity])
            polyLongDiv();
        else
            ++dest;
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
        return reinterpret_cast<std::uint8_t*>(&m_scratch[numParity]);
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
        for (size_t idx = 1; idx <= numParity; ++idx)
            m_scratch[idx - 1] = m_scratch[idx] - m_polynomial.coeff(idx) * quotient;
    }

    // createRootPolynomial(4, GF256Polynomial{1})
    GF256Polynomial<10> m_polynomial;
    GF256Element m_scratch[numParity + 1];
    GF256Element* m_scratchIter;
    bool m_finished;
};


int main()
{
    using namespace std;

#if 1
    auto pwrTable = detail::createPowerTable(0, {});
    cout << pwrTable << endl;

    auto logTable = detail::createLogTable(pwrTable, 0, {});
    cout << logTable << endl;
#endif

    GF256Element a(2);
    GF256Element b(6);
    cout << (a * b).value() << endl;
    cout << (b + b).value() << endl;

    for (int e : {255, 256, 511, 512, 513})
        cout << e << "   " << e % 256 + e / 256 << "   " << e % 255 << endl;

    cout << GF256Element(1) * GF256Element(3) << endl;

    GF256Polynomial<1> pa{1, 1};
    GF256Polynomial<2> pb{3, 4, 5};

    cout << pa << endl;
    cout << pb << endl;
    cout << (pa + pb) << endl;
    cout << (pa * pb) << endl;
}
