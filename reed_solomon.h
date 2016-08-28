#pragma once

#include "galois.h"

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
        :
        //m_scratch(),
        m_scratchIter(&m_scratch[0]),
        m_finished(false)
    {
        std::fill(std::begin(m_scratch), std::end(m_scratch), GF256Element());
    }

    void finish() noexcept
    {
        while (m_scratchIter < &m_scratch[numParity + 1])
            *m_scratchIter++ = 0;
        for (std::size_t count = 0; count < numParity; ++count)
            polyLongDiv();
    }

    void encode(const std::uint8_t* message, std::size_t length) noexcept
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
protected:
    // RS(255,251) polynomial (285 dec)
    const GF256Polynomial<4> m_polynomial = { { 1, 30, 216, 231, 116 } };
    GF256Element m_scratch[numParity + 1] = { 0, };
    GF256Element* m_scratchIter;
    bool m_finished;
};

