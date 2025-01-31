// -----------------------------------------------------------------------
// This file is distributed under a 2-clause BSD like license
// Alternatively this file can be used under GPLv2 license
// See LICENSE.TXT for details.

#ifndef RS11_REEDSOLOMON_HPP
#define RS11_REEDSOLOMON_HPP

#include "_config.hpp"
#include "GaloisField.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

#if defined(_MSC_VER)
    #if _MSC_VER < 1911
        #define CONSTEXPR_WORKAROUND
    #endif
#endif

namespace rs11
{

namespace detail
{

//! Returns the polynomial \f$ p(x) = prod_{i=0}^N (x - \alpha^i) \f$ with
//! \f$ alpha = 2 \f$.
template <std::size_t N>
constexpr
Galois::GF256Polynomial<N> createReedSolomonGeneratorPolynomial()
{
    return createReedSolomonGeneratorPolynomial<N - 1>()
            * Galois::GF256Polynomial<1>({1, Galois::GF256Value::pow2(N - 1)});
}

//! Returns the polynomial \f$ p(x) = 1 \f$.
template <>
constexpr
Galois::GF256Polynomial<0> createReedSolomonGeneratorPolynomial()
{
    return Galois::GF256Polynomial<0>({1});
}

} // namespace detail


template <std::size_t M, std::size_t N>
class ReedSolomonEncoder;

template <std::size_t M, std::size_t N>
class ReedSolomonDecoder;


//! \brief A Reed-Solomon encoder.
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


    //! Creates a Reed-Solomon encoder.
#ifndef CONSTEXPR_WORKAROUND
    constexpr
#endif
    ReedSolomonEncoder() noexcept
        : m_polynomial(detail::createReedSolomonGeneratorPolynomial<numParitySyms>()),
          m_scratch(),
          m_scratchIter{&m_scratch[0]},
          m_size{0},
          m_finished{false}
    {
    }

    //! Creates a Reed-Solomon encoder, with agiven polynom
#ifndef CONSTEXPR_WORKAROUND
    constexpr
#endif
    ReedSolomonEncoder(const Galois::GF256Polynomial<numParitySyms>& poly) noexcept
        : m_polynomial(poly),
        m_scratch(),
        m_scratchIter{ &m_scratch[0] },
        m_size{ 0 },
        m_finished{ false }
    {
    }

    ReedSolomonEncoder(const ReedSolomonEncoder&) = delete;
    ReedSolomonEncoder& operator=(const ReedSolomonEncoder&) = delete;

    //! \brief Encodes a message.
    //!
    //! Encodes the given \p message consisting of \p length bytes. The
    //! encoding is finalized right away.
    void encode(const std::uint8_t* message, std::size_t length)
    {
        encodePart(message, length);
        finish();
    }

    //! \brief Encodes a part of a message.
    //!
    //! Encodes the given \p message consisting of \p length bytes. This method
    //! can be called multiple times until the whole message is complete.
    //! At the end, finish() must be called to finalize the encoding.
    void encodePart(const std::uint8_t* message, std::size_t length) RS11_ERR_NOEXCEPT
    {
        #ifndef RS11_NO_EXCEPTIONS        
            if (m_size + length > N)
                throw std::exception();
            if (m_finished)
                throw std::exception();
        #else
            if ((m_size + length > N) || m_finished)
                return;
        #endif
        
        m_size += length;

        // Fill up the scratch space except the last slot.
        while (m_scratchIter != &m_scratch[numParitySyms] && length)
        {
            *m_scratchIter++ = *message++;
            --length;
        }

        // Perform a long division for every new element.
        while (length--)
        {
            *m_scratchIter = *message++;
            polyLongDiv();
        }
    }

    //! \brief Adds one byte.
    //!
    //! Encodes the given \p datum. When all bytes have been added,
    //! finish() must be called to finalize the encoding.
    ReedSolomonEncoder& operator<<(std::uint8_t datum) RS11_ERR_NOEXCEPT
    {
        #ifndef RS11_NO_EXCEPTIONS
        if (m_size == N)
            throw std::exception();
        if (m_finished)
            throw std::exception();
        #else
            if((m_size == N) || m_finished)
                return 
        #endif

        ++m_size;
        *m_scratchIter = datum;
        if (m_scratchIter == &m_scratch[numParitySyms])
            polyLongDiv();
        else
            ++m_scratchIter;

        return *this;
    }

    #ifdef RS11_HAVE_GSL
    //! \brief Encodes a message.
    //!
    //! Encodes the given \p message. The encoding is automatically finalized
    //! with a call to finish().
    template<typename T, int S>
    void encode(gsl::span<const T, S> message)
    {
        static_assert(sizeof(T) != 8, "RS implemntation only works with bytes");
        encodePart(reinterpret_cast<const uint8_t*>(message.data()), message.size());
        finish();
    }

    //! \brief Encodes a part of a message.
    //!
    //! Encodes the message part given by \p message. When all parts have been
    //! encoded, finish() has to be called.
    template<typename T, int S>
    void encodePart(gsl::span<const T, S> message)
    {
        static_assert(sizeof(T) != 8, "RS implemntation only works with bytes");
        encodePart(message.data(), message.size());
    }
    #endif // RS11_HAVE_GSL

    //! \brief Finishes encoding.
    //!
    //! Finalizes the encoding. This method must be called after calling
    //! encodePart().
    void finish() noexcept
    {
        if (m_finished)
            return;

        m_finished = true;

        while (m_scratchIter < &m_scratch[numParitySyms + 1])
            *m_scratchIter++ = 0;
        for (std::size_t count = 0; count < numParitySyms; ++count)
            polyLongDiv();
    }

    //! \brief Resets the encoder.
    //!
    //! Resets the internal state and prepares the encoder for a new message.
    void reset()
    {
        m_scratchIter = &m_scratch[0];
        m_size = 0;
        m_finished = false;
    }

    //! \brief An iterator to the first parity symbol.
    //!
    //! Returns an iterator to the first parity symbol. This method can only
    //! be called when the encoding has been finalized (either with finish()
    //! or encode()).
    const std::uint8_t* begin() RS11_ERR_NOEXCEPT
    {
        static_assert(sizeof(Galois::GF256Value) == sizeof(std::uint8_t),
                      "Mismatch in the representation.");

        if (!m_finished)
            #ifndef RS11_NO_EXCEPTIONS
                throw std::exception();
            #else
                return nullptr;
            #endif

        return reinterpret_cast<std::uint8_t*>(&m_scratch[0]);
    }

    //! \brief An iterator past the last parity symbol.
    //!
    //! Returns an iterator past the last parity symbol. This method can only
    //! be called when the encoding has been finalized (either with finish()
    //! or encode()).
    const std::uint8_t* end() RS11_ERR_NOEXCEPT
    {
        static_assert(sizeof(Galois::GF256Value) == sizeof(std::uint8_t),
                      "Mismatch in the representation.");

        if (!m_finished)
            #ifndef RS11_NO_EXCEPTIONS
                throw std::exception();
            #else
                return nullptr;
            #endif

        // Spare the last element.
        return reinterpret_cast<std::uint8_t*>(&m_scratch[numParitySyms]);
    }

    //! \brief Returns the number of parity symbols.
    //!
    //! Returns the number of parity symbols which this encoder computes.
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
        Galois::GF256Value quotient = m_scratch[0];
        // When subtracting the intermediate result, the first term is
        // zeroed. Shift the lower order terms into place.
        for (size_t idx = 1; idx <= numParitySyms; ++idx)
            m_scratch[idx - 1] = m_scratch[idx] - m_polynomial[idx] * quotient;
    }

    //! The polynomial for Reed-Solomon encoding.
    Galois::GF256Polynomial<numParitySyms> m_polynomial;

    //! A scratch space to perform the polynomial long division.
    Galois::GF256Value m_scratch[numParitySyms + 1];
    //! An iterator over the scratch space.
    Galois::GF256Value* m_scratchIter;

    //! The number of symbols which have been encoded.
    std::size_t m_size;
    //! Set if the encoding has been finished.
    bool m_finished;
};



//! \brief An error in a message.
//!
//! The Error class encapsulates an error detected during Reed-Solomon decoding.
//! It stores the index of the erronous symbol in the message. It also stores
//! stores the magnitude of the error, which is needed to restore the
//! original value.
class Error
{
public:
    //! \brief The index of the defective symbols.
    //!
    //! Returns the index of the defective symbol.
    std::size_t index() const noexcept
    {
        return m_location;
    }

    //! \brief Corrects a defective symbol.
    //!
    //! Subtracts the error magnitude from the erronous \p datum and returns
    //! the corrected value.
    std::uint8_t correct(std::uint8_t datum) const noexcept
    {
        return std::uint8_t(Galois::GF256Value(datum) - m_magnitude);
    }

private:
    std::uint8_t m_root;
    std::size_t m_location;
    Galois::GF256Value m_magnitude;

    template <std::size_t M, std::size_t N>
    friend class ReedSolomonDecoder;
};

class ErrorRange
{
public:
    ErrorRange(const Error* begin, const Error* end) noexcept
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
    const Error* m_begin;
    const Error* m_end;
};

enum class DecoderResult
{
    Good,
    Correctable,
    Defective
};

//! \brief A Reed-Solomon decoder.
template <std::size_t M, std::size_t N>
class ReedSolomonDecoder
{
    static_assert(M <= 255, "Only 255 symbols are supported");
    static_assert(N < M, "The number of payload symbols must be smaller than the total number of symbols");

    static constexpr std::size_t numParitySyms = M - N;

public:
    using encoder_type = ReedSolomonEncoder<M, N>;

    //! Creates a Reed-Solomon decoder.
    ReedSolomonDecoder()
        : m_result{DecoderResult::Good},
          m_numErrors{0},
          m_size{0},
          m_finished{false}
    {
    }

    ReedSolomonDecoder(const ReedSolomonDecoder&) = delete;
    ReedSolomonDecoder& operator=(const ReedSolomonDecoder&) = delete;

    //! \brief Decodes a message.
    //!
    //! Decodes the given \p message consisting of \p length bytes. The
    //! decoding is finalized right away.
    void decode(const std::uint8_t* message, std::size_t length)
    {
        decodePart(message, length);
        finish();
    }

    //! \brief Decodes a part of a message.
    //!
    //! Decodes the given \p message consisting of \p length bytes. This method
    //! can be called multiple times until the whole message is complete.
    //! At the end, finish() must be called to finalize the decoding.
    void decodePart(const std::uint8_t* message, std::size_t length)
    {
        #ifndef RS11_NO_EXCEPTIONS
            if (m_size + length > M)
                throw std::exception();
            if (m_finished)
                throw std::exception();
        #endif

        m_size += length;
        while (length--)
        {
            Galois::GF256Value x = *message++;
            for (size_t idx = 0; idx < numParitySyms; ++idx)
            {
                // TODO: Speed up with
                // if m_syndrome[idx] == 0: m_syndrome[idx] = datum;
                m_syndromes[idx] = m_syndromes[idx] * Galois::GF256Value::pow2(idx) + x;
            }
        }
    }

    //! \brief Adds another byte to decoding.
    //!
    //! Adds the \p datum to the decoding. The decoding must be finalized
    //! with a call to finish().
    ReedSolomonDecoder& operator<<(std::uint8_t datum) //noexcept
    {
        if (m_size == N)
            throw std::exception();
        if (m_finished)
            throw std::exception();

        ++m_size;

        // To compute the syndromes, the encoded message must be evaluated
        // at the roots of the generator polynomial, i.e. at
        // \alpha^i, i = 0 ... N, with \alpha = 2. If we had no errors in
        // the encoded message, the syndromes would be all zero.
        for (size_t idx = 0; idx < numParitySyms; ++idx)
        {
            // TODO: Speed up with
            // if m_syndrome[idx] == 0: m_syndrome[idx] = datum;
            m_syndromes[idx] = m_syndromes[idx] * Galois::GF256Value::pow2(idx)
                               + Galois::GF256Value(datum);
        }

        return *this;
    }

#ifdef RS11_HAVE_GSL
    //! \brief Decodes a message.
    //!
    //! Decodes the given \p message. The decoding is automatically finalized
    //! with a call to finish().
    void decode(gsl::span<const std::uint8_t> message)
    {
        decodePart(message.data(), message.size());
        finish();
    }
    void decode(gsl::span<const gsl::byte> message)
    {
        decodePart(reinterpret_cast<const uint8_t*>(message.data()), message.size());
        finish();
    }

    //! \brief Decodes a part of a message.
    //!
    //! Decodes the message part given by \p message. When all parts have been
    //! decoded, finish() has to be called.
    void decodePart(gsl::span<const std::uint8_t> message)
    {
        decodePart(message.data(), message.size());
    }
    void decodePart(gsl::span<const gsl::byte> message)
    {
        decodePart(reinterpret_cast<const uint8_t*>(message.data()), message.size()) ;
    }
#endif // RS11_HAVE_GSL

    //! \brief Finishes the decoding.
    void finish() noexcept
    {
        using namespace Galois;

        if (m_finished)
            return;

        m_finished = true;

        // If all syndromes are zero, there is no error. In that case, we
        // can return immediately. If we find an erronous syndrome, assume
        // that the message is defective and cannot be corrected. Only set it
        // to correctable if we know for sure that correction is feasible.
        m_result = DecoderResult::Good;
        for (auto syn : m_syndromes)
            if (syn)
            {
                m_result = DecoderResult::Defective;
                break;
            }
        if (m_result == DecoderResult::Good)
            return;


        for (unsigned idx = 1; idx <= numParitySyms; ++idx)
            m_errorLocator[idx] = m_B[idx] = 0;
        m_errorLocator[0] = m_B[0] = GF256Value(1);

        // Berlekamp-Massey algorithm
        std::size_t numErrors = 0;
        for (std::size_t iteration = 0; iteration < numParitySyms; ++iteration)
        {
            GF256Value discrepancy;
            for (std::size_t idx = 0; idx <= iteration; ++idx)
                discrepancy += m_errorLocator[idx] * m_syndromes[iteration - idx];

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
        }

        // Find the roots of the error locator polynomial. The roots are the
        // inverses of the error locations.
        // TODO: Chien search was faster, if we stored the error locator
        // polynomial in logarithmic form. We would spare a lot of table
        // look-ups for the multiplications.
        auto errorLocatorDegree = m_errorLocator.degree();
        for (unsigned idx = 1; idx <= 255; ++idx)
        {
            auto x = GF256Value::pow2(idx);
            if (!m_errorLocator(x, errorLocatorDegree))
            {
                auto invIndex = 255 - idx;
                if (m_size - 1 < invIndex)
                    return; // TODO: unable to correct the error

                m_errors[m_numErrors].m_root = x.value();
                m_errors[m_numErrors].m_location = m_size - 1 - invIndex;
                if (++m_numErrors == errorLocatorDegree)
                    break;
            }
        }

        // If not all roots of the error locator polynomial were found (the
        // number of roots is smaller than the degree), the message cannot
        // be corrected.
        if (m_numErrors != errorLocatorDegree)
            return; // TODO: unable to correct the error

        for (unsigned idx = 0; idx <= numParitySyms; ++idx)
            m_errorEvaluator[idx] = 0;

        // Compute the error evaluator polynomial:
        // omega(x) = s(x) * lambda(x) (mod x^numParitySyms).
        // Use a specially crafted multiplication which omits the higher order
        // terms as they will vanish after the modulo operation.
        for (unsigned idx = 0; idx < errorLocatorDegree; ++idx)
        {
            GF256Value sum(0);
            for (unsigned k = 0; k <= idx; ++k)
                sum += m_syndromes[k] * m_errorLocator.coeff(idx - k);
            m_errorEvaluator[idx] = sum;
        }

        // Apply Forney's algorithm to compute the error magnitudes.
        for (unsigned errorIdx = 0; errorIdx < m_numErrors; ++errorIdx)
        {
            GF256Value x = m_errors[errorIdx].m_root;
            GF256Value numerator = m_errorEvaluator(x);
            // The denominator is the formal derivative of the error locator
            // polynomial evaluated at the root.
            GF256Value denominator = m_errorLocator(x, derivative);
            if (!denominator)
                return; // TODO: unable to correct the error

            m_errors[errorIdx].m_magnitude
                = (numerator * (denominator * x).inverse()).value();
        }

        m_result = DecoderResult::Correctable;
    }

    //! \brief Resets the decoder.
    //!
    //! Resets the internal state and prepares the decoder for a new message.
    void reset() noexcept
    {
        for (auto& syn : m_syndromes)
            syn = 0;

        m_result = DecoderResult::Good;
        m_numErrors = 0;
        m_size = 0;
        m_finished = false;
    }

    //! \brief Returns the errors.
    //!
    //! Returns a range of errors. This method must only be called after the
    //! decoding has been finalized with a call to finish() and if the
    //! result is not DecoderResult::Defective. If the message did not
    //! contain errors, an empty range is returned.
    ErrorRange errors() const
    {
            if (!m_finished)
                throw std::exception();
            if (m_result == DecoderResult::Defective)
                throw std::exception();

        return ErrorRange(&m_errors[0], &m_errors[0] + m_numErrors);
    }

    //! \brief Returns the result of the decoding.
    DecoderResult result() const
    {
        if (!m_finished)
        {
            #ifdef RS11_NO_EXCEPTIONS
                return DecoderResult::Defective;
            #else
                throw std::exception();
            #endif
        }
        return m_result;
    }

    //! \brief Returns the number of parity symbols.
    //!
    //! Returns the number of parity symbols.
    std::size_t numParitySymbols() const noexcept
    {
        return numParitySyms;
    }

public:
    //! The syndromes.
    Galois::GF256Value m_syndromes[numParitySyms];

    //! The error locator polynomial.
    Galois::GF256Polynomial<numParitySyms> m_errorLocator;
    //! The error evaluator polynomial.
    Galois::GF256Polynomial<numParitySyms> m_errorEvaluator;

    //! A scratch polynomial needed to update the error locator polynomial.
    Galois::GF256Polynomial<numParitySyms> m_temp;

    Galois::GF256Polynomial<numParitySyms> m_B;

    Error m_errors[numParitySyms];

    //! The result of the decoding.
    DecoderResult m_result;
    //! The number of errors, which have been detected during decoding.
    std::size_t m_numErrors;

    //! The number of symbols which have been decoded.
    std::size_t m_size;
    //! Set if the decoding has been finished.
    bool m_finished;
};

} // namespace rs11

#endif // RS11_REEDSOLOMON_HPP
