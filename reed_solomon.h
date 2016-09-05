#pragma once

#include "galois.h"

// Polynomial
// \prod_{i=0}^{?} (x - \alpha^i), \alpha \in GF(2^8)
template <std::size_t M, std::size_t N>
class ReedSolomon
{
    static_assert(N < M, "");
    static constexpr std::size_t numParitySyms = M - N;

public:
    constexpr
    ReedSolomon() noexcept
        : m_scratchIter(&m_scratch[0]),
          m_finished(false)
    {
        std::fill(std::begin(m_scratch), std::end(m_scratch), GF256Element());
    }

    ReedSolomon(const ReedSolomon&) = delete;
    ReedSolomon& operator=(const ReedSolomon&) = delete;

    void finish() noexcept
    {
        while (m_scratchIter != &m_scratch[numParitySyms + 1])
            *m_scratchIter++ = 0;
        for (std::size_t count = 0; count < numParitySyms; ++count)
            polyLongDiv();
    }

    void encode(const std::uint8_t* message, std::size_t length) noexcept
    {
        encodePart(message, length);
        finish();
    }

    void encodePart(const std::uint8_t* message, std::size_t length) noexcept
    {
        while (m_scratchIter != &m_scratch[numParitySyms] && length)
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
    }

    GF256Polynomial<4> syndrome(const std::uint8_t (&rx_message) [M]) const noexcept
    {
        GF256Polynomial<4> S;
        static constexpr uint8_t scale[] = {16,8,4,2};
        for (auto d : rx_message)
        {
            for (int i = 0; i < 4; i++)
                S.coeff(i) = S.coeff(i) * scale[i] + d;
        }
        return S;
    }

    ReedSolomon& operator<<(std::uint8_t datum) noexcept
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
        for (std::size_t idx = 1; idx <= numParitySyms; ++idx)
            m_scratch[idx - 1] = m_scratch[idx] - m_polynomial.coeff(idx) * quotient;
    }
protected:
    // RS(255,251) polynomial (285 dec)
    const GF256Polynomial<4> m_polynomial = { { 1, 30, 216, 231, 116 } };
    GF256Element m_scratch[numParitySyms + 1] = { 0, };
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
