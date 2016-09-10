#pragma once



#include <utility>
#include <cstdint>

#define MAX_IDX_SEQ   256

constexpr std::uint8_t reduction_poly = 0x1d;

constexpr
int nextPower(int prev)
{
    return (prev & 0x80) ? (((prev ^ 0x80) << 1) ^ reduction_poly) : (prev << 1);
}





template <typename... TArgs>
using index_sequence_for = std::make_index_sequence<sizeof...(TArgs)>;


namespace detail
{

    class IncrArray
    {
        template <std::size_t... TIndices>
        constexpr IncrArray doModify(std::size_t idx, std::uint8_t value, std::index_sequence<TIndices...>) const
        {
            return{ (TIndices == idx ? value : m_data[TIndices])... };
        }

    public:
        constexpr IncrArray modified(std::size_t idx, int value) const
        {
            return doModify(idx, value, std::make_index_sequence<MAX_IDX_SEQ>());
        }

        constexpr std::size_t size() const
        {
            return MAX_IDX_SEQ;
        }

        constexpr std::uint8_t operator[](std::size_t idx) const
        {
            return m_data[idx];
        }

        std::uint8_t m_data[MAX_IDX_SEQ];
    };

    constexpr auto createPowerTable(std::size_t idx, const IncrArray& temp) -> IncrArray
    {
        return idx == MAX_IDX_SEQ ? temp
            : createPowerTable(
                idx + 1,
                temp.modified(idx, idx == 0 ? 1 : nextPower(temp[idx - 1])));
    }

    constexpr IncrArray createLogTable(const IncrArray& powers, std::size_t idx, const IncrArray& temp)
    {
        return idx == MAX_IDX_SEQ ? temp
            : createLogTable(
                powers, idx + 1,
                temp.modified(powers[idx], powers[idx] == 1 ? 0 : int(idx)));
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
    constexpr GF256Element() noexcept
        : m_value(0)
    {
    }

    constexpr GF256Element(std::uint8_t value) noexcept
        : m_value(value)
    {
    }

    //! Computes the sum of two elements from GF(2^8).
    //!
    //! Adds \p b to this element and returns the sum.
    constexpr GF256Element operator+(GF256Element b) const noexcept
    {
        // Addition in GF(2) is the XOR function.
        return m_value ^ b.m_value;
    }

    GF256Element& operator+=(GF256Element b) noexcept
    {
        m_value ^= b.m_value;
        return *this;
    }

    //! Computes the difference between two elements from GF(2^8).
    //!
    //! Subtracts \p b from this element and returns the difference.
    constexpr GF256Element operator-(GF256Element b) const noexcept
    {
        // Subtraction in GF(2) is just an addition.
        return m_value ^ b.m_value;
    }

    //! Computes the product of two elements from GF(2^8).
    //!
    //! Multiplies this element with \p b and returns the product.
    constexpr GF256Element operator*(GF256Element b) const noexcept
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
    constexpr unsigned value() const noexcept
    {
        return m_value;
    }

    constexpr explicit operator std::uint8_t() const noexcept
    {
        return m_value;
    }

    explicit operator bool() const noexcept
    {
        return m_value != 0;
    }

    //! Returns 2 raised to the power \p p.
    //!
    //! Returns the GF(2^8) element <tt>2^p</tt>.
    static constexpr GF256Element pow2(unsigned p)
    {
        return m_antiLogTable[p % 255];
    }

private:
    std::uint8_t m_value;

    static constexpr unsigned mod(unsigned x) noexcept
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
constexpr T const& cmax(T const& a, T const& b)
{
    return a > b ? a : b;
}

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
    constexpr GF256Polynomial<sizeof...(TIndices) -1> doAdd(const GF256Polynomial<TDeg>& b, std::index_sequence<TIndices...>) const noexcept
    {
        return{ (coeff(TIndices) + b.coeff(TIndices))... };
    }

    template <std::size_t TDeg>
    constexpr GF256Element doConvPart(const GF256Polynomial<TDeg>& b, std::size_t degA, std::size_t degB) const noexcept
    {
        return degA == 0 ? coeff(degA) * b.coeff(degB)
            : coeff(degA) * b.coeff(degB)
            + doConvPart(b, degA - 1, degB + 1);
    }

    template <std::size_t TDeg, std::size_t... TIndices>
    constexpr GF256Polynomial<sizeof...(TIndices) -1> doMul(const GF256Polynomial<TDeg>& b, std::index_sequence<TIndices...>) const noexcept
    {
        return{ doConvPart(b, TIndices, 0)... };
    }

public:
    //    constexpr
    //    GF256Polynomial() noexcept
    //        : _m_coefficients{}
    //    {
    //    }

    template <std::size_t TDeg>
    constexpr auto operator+(const GF256Polynomial<TDeg>& b) const noexcept -> GF256Polynomial<cmax(TDegree, TDeg)>
    {
        return doAdd(b, std::make_index_sequence<cmax(TDegree, TDeg) + 1>());
    }

    template <std::size_t TDeg>
    constexpr auto operator-(const GF256Polynomial<TDeg>& b) const noexcept -> GF256Polynomial<cmax(TDegree, TDeg)>
    {
        // Subtraction in GF(2^8) equals addition.
        return doAdd(b, std::make_index_sequence<cmax(TDegree, TDeg) + 1>());
    }

    //! Returns the product of two polynomials.
    template <std::size_t TDeg>
    constexpr auto operator*(const GF256Polynomial<TDeg>& b) const noexcept -> GF256Polynomial<TDegree + TDeg>
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

    //! Evaluates the polynomial of a certain degree.
    GF256Element operator()(GF256Element x, std::size_t degree) const noexcept
    {
        GF256Element sum{ _m_coefficients[degree] };
        for (std::size_t deg = degree; deg > 0; --deg)
            sum = sum * x + _m_coefficients[deg - 1];
        return sum;
    }

    GF256Element operator()(GF256Element x, derivative_t) const noexcept
    {
        if (TDegree < 1)
            return GF256Element(0);

        {
            GF256Polynomial tmp;
            for (unsigned i = 0; i < TDegree; i += 2)
                tmp[i] = _m_coefficients[i + 1];
            //cout << "> > > > derivative should be " << tmp(x).value() << endl;
        }
        // 3 + 4 x + 5 x^2
        // 4 + 5 * 2 * x

        // The even powers of the polynomial lead to the derivative
        // (p * x^n)' = n * p * x^(n-1), with n even. But we can write this as
        // n/2 * p * x^(n-1) + n/2 * p * x^(n-1) and so it is equal to
        // zero (because we are in GF(2) and addition is a XOR).
        //
        // The odd powers result in n * p * x^(n-1), with n odd. This reduces to
        // (n-1)/2 * p * x^(n-1) + (n-1)/2 * p * x^(n-1) + p * x^(n-1) =
        // p * x^(n-1). So only odd powers survive.
        GF256Element x2 = x * x;
        GF256Element sum(0);
        for (int even = (TDegree - 1) & ~1; even >= 0; even -= 2)
            sum = sum * x2 + _m_coefficients[even + 1];
        //cout << "> > > > derivative is " << sum.value() << endl;

        sum = GF256Element(0);
        for (int idx = TDegree; idx >= 1; idx -= 1)
            sum = sum * x + ((idx % 2 == 0) ? GF256Element(0) : _m_coefficients[idx]);
        //cout << "> > > > derivative chk " << sum.value() << endl;

        return sum;
    }

    //! Returns the coefficient of a monomial.
    //!
    //! Returns the coefficient of the monomial of given \p degree. Zero is
    //! returned if \p degree is higher than the degree of the polynomial.
    constexpr GF256Element coeff(std::size_t degree) const noexcept
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
    constexpr GF256Polynomial<N> createReedSolomonGeneratorPolynomial()
    {
        return createReedSolomonGeneratorPolynomial<N - 1>()
            * GF256Polynomial<1>({ 1, GF256Element::pow2(N - 1) });
    }

    //! Returns the polynomial \f$ p(x) = 1 \f$.
    template <>
    constexpr GF256Polynomial<0> createReedSolomonGeneratorPolynomial()
    {
        return GF256Polynomial<0>({ 1 });
    }

} // namespace detail
