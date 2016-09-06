#pragma once

#include "galois.h"


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



template <std::size_t M, std::size_t N>
class Decoder
{
    static_assert(N < M, "");
    static constexpr std::size_t numParitySyms = M - N;

public:
    Decoder()
    {
    }

    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;

    Decoder& operator<<(std::uint8_t datum) noexcept
    {
        return *this;
    }

private:
    //! The syndromes of the encoded data.
    GF256Element m_syndromes[numParitySyms];
    //! The error locator polynomial.
    GF256Polynomial<numParitySyms> m_errorLocator;
    //! A scratch polynomial for the Berlekamp-Massey algorithm.
    GF256Polynomial<numParitySyms> m_temp;
};
