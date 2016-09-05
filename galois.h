#pragma once

#include <algorithm>
#include <utility>
#include <ostream>

namespace detail
{

constexpr
int nextPower(int prev)
{
    return (prev & 0x80) ? (((prev ^ 0x80) << 1) ^ 0x1d) : (prev << 1);
}

template <typename... TArgs>
using index_sequence_for = std::make_index_sequence<sizeof...(TArgs)>;

static constexpr uint32_t MAX = 256;

class IncrArray
{
    template <std::size_t... TIndices>
    constexpr
        IncrArray doModify(std::size_t idx, std::uint8_t value,
            std::index_sequence<TIndices...>) const
    {
        return{ (TIndices == idx ? value : m_data[TIndices])... };
    }

public:
    constexpr
        IncrArray modified(std::size_t idx, int value) const
    {
        return doModify(idx, value, std::make_index_sequence<MAX>());
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

#if 0
    ostream& operator<< (ostream& str, const IncrArray& a)
    {
        for (std::size_t idx = 0; idx < a.size(); ++idx)
            str << setw(3) << idx << "   " << setw(3) << unsigned(a[idx]) << endl;
        return str;
    }
#endif

constexpr
auto createPowerTable(std::size_t idx, const IncrArray& temp) -> IncrArray
{
    return idx == MAX ? temp
        : createPowerTable(
            idx + 1,
            temp.modified(idx, idx == 0 ? 1 : nextPower(temp[idx - 1])));
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

    //! Computes the product between two elements from GF(2^8).
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

    constexpr
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

    //    //! Returns the 2 raised to the power \p p.
    //    //!
    //    //! Returns the GF(2^8) element 2^p.
    //    static GF256Element exp(unsigned p)
    //    {
    //        return m_antiLogTable[trueMod(p)];
    //    }

    //! Returns \p true, if this element is equal to \p b.
    constexpr
    bool operator==(GF256Element b) const noexcept
    {
        return m_value == b.m_value;
    }

    //! Returns \p true, if this element is non-zero.
    constexpr explicit
    operator bool() const noexcept
    {
        return m_value != 0;
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
    GF256Polynomial<sizeof...(TIndices) - 1> doAdd(
        const GF256Polynomial<TDeg>& b,
        std::index_sequence<TIndices...>) const noexcept
    {
        return{ (coeff(TIndices) + b.coeff(TIndices))... };
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
    GF256Polynomial<sizeof...(TIndices) - 1> doMul(
        const GF256Polynomial<TDeg>& b,
        std::index_sequence<TIndices...>) const noexcept
    {
        return{ doConvPart(b, TIndices, 0)... };
    }

public:
    //constexpr
    //GF256Polynomial() noexcept
    //    : _m_coefficients{}
    //{
    //}

    constexpr
    GF256Polynomial(const GF256Element(&poly)[TDegree + 1]) noexcept
    {
        using namespace std;
        copy(begin(poly), end(poly), begin(_m_coefficients));
    }

    template <std::size_t TDeg>
    constexpr
    GF256Polynomial<cmax(TDegree, TDeg)> operator+(const GF256Polynomial<TDeg>& b) const noexcept
    {
        return doAdd(b, std::make_index_sequence<cmax(TDegree, TDeg) + 1>());
    }

    template <std::size_t TDeg>
    constexpr
    auto operator-(const GF256Polynomial<TDeg>& b) const noexcept
        -> GF256Polynomial<cmax(TDegree, TDeg)>
    {
        // Subtraction in GF(2^8) equals addition.
        return doAdd(b, std::make_index_sequence<cmax(TDegree, TDeg) + 1>());
    }

    //! Returns the product of two polynomials.
    template <std::size_t TDeg>
    constexpr
    auto operator*(const GF256Polynomial<TDeg>& b) const noexcept
        -> GF256Polynomial<TDegree + TDeg>
    {
        return doMul(b, std::make_index_sequence<TDegree + TDeg + 1>());
    }

    //! Evaluates the polynomial.
    GF256Element operator()(GF256Element x) const noexcept
    {
        GF256Element sum{ _m_coefficients[TDegree] };
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

    GF256Element& coeff(std::size_t degree)
    {
        return _m_coefficients[degree];
    }

    //! Multiplies this polynomial inplace with the polynomial \f$ p(x) = x \f$.
    void timesX() noexcept
    {
        // Shift the coefficients: [a, b, c, d] -> [d, a, b, c].
        std::rotate(&_m_coefficients[0], &_m_coefficients[0] + TDegree,
                    &_m_coefficients[0] + TDegree + 1);
        _m_coefficients[0] = 0;
    }

    constexpr
    bool operator==(const GF256Polynomial& rhs) const noexcept
    {
        using namespace std;
        return equal(begin(_m_coefficients), end(_m_coefficients), begin(rhs._m_coefficients));
    }

protected:
    // The coefficients of the polynomial.
    GF256Element _m_coefficients[TDegree + 1] = { 0, };
};

template <std::size_t N>
constexpr
GF256Polynomial<N> prod()
{
    return GF256Polynomial<1>({ 1, N }) * prod<N - 1>();
}

template <>
constexpr
GF256Polynomial<0> prod<0>()
{
    return GF256Polynomial<0>({ 1 });
}



std::ostream& operator<<(std::ostream& str, GF256Element e)
{
    str << e.value();
    return str;
}

template <std::size_t TDeg>
std::ostream& operator<<(std::ostream& str, const GF256Polynomial<TDeg>& p)
{
    for (std::size_t i = 0; i <= TDeg; ++i)
    {
        if (i)
            str << " + ";
        str << p.coeff(i) << " x^" << i;
    }
    return str;
}

