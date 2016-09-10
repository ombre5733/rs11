#pragma once

#include "galois.h"

template <std::size_t M, std::size_t N>
class ReedSolomonEncoder;

template <std::size_t M, std::size_t N>
class ReedSolomonDecoder;


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
    using decoder_type = ReedSolomonDecoder<M, N>;

    #ifndef _MSC_VER
    // MSVC does not like the constexpr in here ...
    // error C2476: 'constexpr' constructor does not initialize all members
    // m_scratch not initialized 
    constexpr
    #endif
        ReedSolomonEncoder() noexcept
        : m_polynomial{ detail::createReedSolomonGeneratorPolynomial<numParitySyms>() },
        m_scratchIter{ &m_scratch[0] },
        m_finished{ false }
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

    std::uint8_t correct(std::uint8_t /*datum*/) const noexcept
    {
        return 0; // TODO
    }

    //private:
    std::uint8_t m_root;
    std::size_t m_location;
    GF256Element m_magnitude;
};

class ErrorRange
{
public:
    ErrorRange(Error* begin, Error* end)
        : m_begin(begin),
        m_end(end)
    {
    }

    ErrorRange(const ErrorRange&) = default;
    ErrorRange& operator=(const ErrorRange&) = default;

    const Error* begin() const noexcept
    {
        return m_begin;
    }

    const Error* end() const noexcept
    {
        return m_end;
    }

private:
    Error* m_begin;
    Error* m_end;
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
    using encoder_type = ReedSolomonEncoder<M, N>;

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
            m_syndromes[idx] = m_syndromes[idx] * GF256Element::pow2(unsigned(idx)) + GF256Element(datum);
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
                discrepancy += m_errorLocator[idx] * m_syndromes[iteration - idx];

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
        // inverses of the error locations.
        unsigned m_numErrors = 0;
        auto errorLocatorDegree = m_errorLocator.degree();
        for (unsigned idx = 1; idx <= 255; ++idx)
        {
            auto x = GF256Element::pow2(idx);
            if (!m_errorLocator(x, errorLocatorDegree))
            {
                cout << "   - Root: i = " << idx
                    << "   (" << 255 - idx << ")"
                    << "   " << x.value()
                    << "   " << x.inverse().value() << endl;
                cout << " DBG " << m_errorLocator(0).value() << endl;

                auto invIndex = 255 - idx;
                if (m_size - 1 < invIndex)
                    return; // TODO: unable to correct the error

                cout << "Error at " << m_size - 1 - invIndex << "\n";

                m_errors[m_numErrors].m_root = x.value();
                m_errors[m_numErrors].m_location = m_size - 1 - invIndex;
                if (++m_numErrors == errorLocatorDegree)
                    break;
            }
        }

        if (m_numErrors != errorLocatorDegree)
            return; // TODO: unable to correct the error

                    // Compute the error evaluator polynomial:
                    // omega(x) = s(x) * lambda(x) (mod x^numParitySyms).
                    // Use a specially crafted multiplication which omits the higher order
                    // terms as they will vanish after the modulo operation.
        for (unsigned idx = 0; idx < errorLocatorDegree; ++idx)
        {
            GF256Element sum(0);
            for (unsigned k = 0; k <= idx; ++k)
                sum += m_syndromes[k] * m_errorLocator.coeff(idx - k);
            m_errorEvaluator[idx] = sum;
        }

        // Apply Forney's algorithm to compute the error magnitudes.
        for (unsigned errorIdx = 0; errorIdx < m_numErrors; ++errorIdx)
        {
            cout << "\n";
            cout << "---- Forney " << errorIdx << endl;

            GF256Element Xi = m_errors[errorIdx].m_root;

            cout << "eval at " << Xi.value() << endl;
            GF256Element numerator = m_errorEvaluator(Xi);
            cout << "numerology " << numerator.value() << endl;

            // The denominator is the formal derivative of the error locator
            // polynomial evaluated at the root.
            GF256Element denominator = m_errorLocator(Xi, derivative);
            cout << "denominator " << denominator.value() << endl;

            if (!denominator)
                return; // TODO: unable to correct the error

            cout << "quot: " << (numerator * denominator.inverse()).value() << endl;
            cout << "??? " << (numerator * (denominator * Xi).inverse()).value() << endl;
            cout << "----\n";

            m_errors[errorIdx].m_magnitude = (numerator * (denominator * Xi).inverse()).value();
        }
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

    GF256Polynomial<numParitySyms> m_errorEvaluator;

    //! A scratch polynomial needed to update the error locator polynomial.
    GF256Polynomial<numParitySyms> m_temp;

    GF256Polynomial<numParitySyms> m_B;

    Error m_errors[numParitySyms];

    std::size_t m_size;
};


