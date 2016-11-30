// -----------------------------------------------------------------------
// This file is distributed under a 2-clause BSD like license
// Alternatively this file can be used under GPLv2 license
// See LICENSE.TXT for details.

#ifndef RS11_GALOISFIELD_HPP
#define RS11_GALOISFIELD_HPP

#include "_config.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>


namespace rs11
{
namespace Galois
{


static constexpr std::size_t num_elements = 256;
static constexpr std::uint8_t reduction_poly = 0x1d;

constexpr
int nextPower(int prev)
{
    return (prev & 0x80) ? (((prev ^ 0x80) << 1) ^ reduction_poly) : (prev << 1);
}

//! A lookup table for logarithms and exponents over a Galois field.
class Table
{
    template <std::size_t... TIndices>
    constexpr
    Table doModify(std::size_t idx, std::uint8_t value,
                   std::index_sequence<TIndices...>) const noexcept
    {
        return { (TIndices == idx ? value : m_data[TIndices])... };
    }

public:
    constexpr
    Table modified(std::size_t idx, int value) const noexcept
    {
        return doModify(idx, value, std::make_index_sequence<num_elements>());
    }

    constexpr
    std::size_t size() const noexcept
    {
        return num_elements;
    }

    constexpr
    std::uint8_t operator[](std::size_t idx) const noexcept
    {
        return m_data[idx];
    }

    std::uint8_t m_data[num_elements];
};

constexpr
auto createPowerTable(std::size_t idx, const Table& temp) -> Table
{
    return idx == num_elements ? temp
                      : createPowerTable(
                            idx + 1,
                            temp.modified(idx, idx == 0 ? 1 : nextPower(temp[idx-1])));
}

constexpr
Table createLogTable(const Table& powers, std::size_t idx, const Table& temp)
{
    return idx == num_elements ? temp
                      : createLogTable(
                            powers, idx + 1,
                            temp.modified(powers[idx], powers[idx] == 1 ? 0 : idx));
}

// The entry at index i equals \alpha^i modulo the reduction polynomial.
static constexpr Table g_antiLogTable = createPowerTable(0, {});
// The entry at index i is j such that \alpha^j = i modulo the reduction polynomial.
static constexpr Table g_logTable = createLogTable(g_antiLogTable, 0, {});



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
class GF256Value
{
public:
    //! Creates the zero element.
    constexpr
    GF256Value() noexcept
        : m_value(0)
    {
    }

    constexpr
    GF256Value(std::uint8_t value) noexcept
        : m_value(value)
    {
    }

    //! Computes the sum of two elements from GF(2^8).
    //!
    //! Adds \p b to this element and returns the sum.
    constexpr
    GF256Value operator+(GF256Value b) const noexcept
    {
        // Addition in GF(2) is the XOR function.
        return m_value ^ b.m_value;
    }

    GF256Value& operator+=(GF256Value b) noexcept
    {
        m_value ^= b.m_value;
        return *this;
    }

    //! Computes the difference between two elements from GF(2^8).
    //!
    //! Subtracts \p b from this element and returns the difference.
    constexpr
    GF256Value operator-(GF256Value b) const noexcept
    {
        // Subtraction in GF(2) is just an addition.
        return m_value ^ b.m_value;
    }

    //! Computes the product of two elements from GF(2^8).
    //!
    //! Multiplies this element with \p b and returns the product.
    constexpr
    GF256Value operator*(GF256Value b) const noexcept
    {
        // Exploit the relation a * b = g^{log_g(a) + log_g(b)}
        return (m_value == 0 || b.m_value == 0)
               ? 0
               : g_antiLogTable[mod(g_logTable[m_value] + g_logTable[b.m_value])];
    }

    GF256Value inverse() const noexcept
    {
        return g_antiLogTable[255 - g_logTable[m_value]];
    }

    // TODO: Remove
    constexpr
    unsigned value() const noexcept
    {
        return m_value;
    }

    //! \brief Compares two elements from GF(2^8).
    //!
    //! Returns \p true, if this element is equal to \p b.
    constexpr
    bool operator==(GF256Value b) const noexcept
    {
        return m_value == b.m_value;
    }

    //! \brief Compares two elements from GF(2^8).
    //!
    //! Returns \p true, if this element is not equal to \p b.
    constexpr
    bool operator!=(GF256Value b) const noexcept
    {
        return m_value != b.m_value;
    }

    constexpr
    explicit operator std::uint8_t() const noexcept
    {
        return m_value;
    }

    //! \brief Compares the element against zero.
    //!
    //! Returns \p true, if this element is non-zero.
    explicit
    operator bool() const noexcept
    {
        return m_value != 0;
    }

    //! Returns 2 raised to the power \p p.
    //!
    //! Returns the GF(2^8) element <tt>2^p</tt>.
    static constexpr
    GF256Value pow2(unsigned p)
    {
        return g_antiLogTable[p % 255];
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
};



struct derivative_t
{
};

constexpr derivative_t derivative;
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
          std::index_sequence<TIndices...>) const noexcept
    {
        return { (coeff(TIndices) + b.coeff(TIndices))... };
    }

    template <std::size_t TDeg>
    constexpr
    GF256Value doConvPart(const GF256Polynomial<TDeg>& b,
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
          std::index_sequence<TIndices...>) const noexcept
    {
        return { doConvPart(b, TIndices, 0)... };
    }

public:
//    constexpr
//    GF256Polynomial() noexcept
//        : _m_coefficients{}
//    {
//    }
#ifndef _MSC_VER
    template <std::size_t TDeg>
    constexpr
    auto operator+(const GF256Polynomial<TDeg>& b) const noexcept
        -> GF256Polynomial<(TDegree > TDeg) ? TDegree : TDeg>
    {
        return doAdd(b, std::make_index_sequence<((TDegree > TDeg) ? TDegree : TDeg) + 1>());
    }


    template <std::size_t TDeg>
    constexpr
    auto operator-(const GF256Polynomial<TDeg>& b) const noexcept
        -> GF256Polynomial<(TDegree > TDeg) ? TDegree : TDeg>
    {
        // Subtraction in GF(2^8) equals addition.
        return doAdd(b, std::make_index_sequence<((TDegree > TDeg) ? TDegree : TDeg) + 1>());
    }
#endif
    //! Returns the product of two polynomials.
    template <std::size_t TDeg>
    constexpr
    auto operator*(const GF256Polynomial<TDeg>& b) const noexcept
        -> GF256Polynomial<TDegree + TDeg>
    {
        return doMul(b, std::make_index_sequence<TDegree + TDeg + 1>());
    }


    //! Evaluates the polynomial.
    GF256Value operator()(GF256Value x) const noexcept
    {
        GF256Value sum{_m_coefficients[TDegree]};
        for (std::size_t deg = TDegree; deg > 0; --deg)
            sum = sum * x + _m_coefficients[deg - 1];
        return sum;
    }

    //! Evaluates the polynomial of a certain degree.
    GF256Value operator()(GF256Value x, std::size_t degree) const noexcept
    {
        GF256Value sum{_m_coefficients[degree]};
        for (std::size_t deg = degree; deg > 0; --deg)
            sum = sum * x + _m_coefficients[deg - 1];
        return sum;
    }

    GF256Value operator()(GF256Value x, derivative_t) const noexcept
    {
        if (TDegree < 1)
            return GF256Value(0);

        // The even powers of the polynomial lead to the derivative
        // (p * x^n)' = n * p * x^(n-1), with n even. But we can write this as
        // n/2 * p * x^(n-1) + n/2 * p * x^(n-1) and so it is equal to
        // zero (because we are in GF(2) and addition is a XOR).
        //
        // The odd powers result in n * p * x^(n-1), with n odd. This reduces to
        // (n-1)/2 * p * x^(n-1) + (n-1)/2 * p * x^(n-1) + p * x^(n-1) =
        // p * x^(n-1). So only odd powers survive.
        GF256Value x2 = x * x;
        GF256Value sum(0);
        for (int even = (TDegree - 1) & ~1; even >= 0; even -= 2)
            sum = sum * x2 + _m_coefficients[even + 1];

        return sum;
    }

    //! Returns the coefficient of a monomial.
    //!
    //! Returns the coefficient of the monomial of given \p degree. Zero is
    //! returned if \p degree is higher than the degree of the polynomial.
    constexpr
    GF256Value coeff(std::size_t degree) const noexcept
    {
        return degree <= TDegree ? _m_coefficients[degree] : GF256Value();
    }

    //! Returns the coefficient of a monomial.
    GF256Value& operator[](std::size_t degree) noexcept
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


//#endif

    // The coefficients of the polynomial.
    GF256Value _m_coefficients[TDegree + 1];
};
} // namespace Galois
} // namespace rs11

#endif // RS11_GALOISFIELD_HPP
