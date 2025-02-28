// -----------------------------------------------------------------------
// This file is distributed under a 2-clause BSD like license
// Alternatively this file can be used under GPLv2 license
// See LICENSE.TXT for details.

#include "catch.hpp"

#include "../ReedSolomon.hpp"
#include "testvectors.h"

#include <algorithm>

TEST_CASE("RS encoding", "[rs]")
{
    std::vector<std::uint8_t> data({0x48, 0x65, 0x6c, 0x6c, 0x6f, 0x57, 0x6f, 0x72, 0x6c, 0x64});

    using Encoder = rs11::ReedSolomonEncoder<255, 251>;
    Encoder encoder;
    REQUIRE(encoder.numParitySymbols() == 4);

    encoder.encode(&data[0], data.size());
    std::vector<std::uint8_t> parity1;
    std::copy(encoder.begin(), encoder.end(), std::back_inserter(parity1));

    encoder.reset();
    encoder.encodePart(&data[0], 5);
    encoder.encodePart(&data[5], data.size() - 5);
    encoder.finish();
    std::vector<std::uint8_t> parity2;
    std::copy(encoder.begin(), encoder.end(), std::back_inserter(parity2));

    encoder.reset();
    for (auto d : data)
        encoder << d;
    encoder.finish();
    std::vector<std::uint8_t> parity3;
    std::copy(encoder.begin(), encoder.end(), std::back_inserter(parity3));

    REQUIRE(parity1 == parity2);
    REQUIRE(parity1 == parity3);
}

TEST_CASE("RS syndromes", "[rs]")
{
    constexpr rs11::Galois::GF256Polynomial<4> POLY_255_251 = { { 1, 30, 216, 231, 116 } };
    for (auto &r : testVectors)
    {
        rs11::ReedSolomonEncoder<255, 251> rs(POLY_255_251);
        gsl::span<const std::uint8_t,251> d = { &r[0], 251 };
        rs.encode(d);
        REQUIRE(memcmp(rs.begin(), &r[251], 4) == 0);

        rs11::ReedSolomonDecoder<255, 251> rs_d;
        rs_d.decode({ r });
        auto res = rs_d.result();
        REQUIRE((res == rs11::DecoderResult::Correctable || res == rs11::DecoderResult::Good));
        REQUIRE(rs_d.errors().begin() == rs_d.errors().end());
    }

}

TEST_CASE("RS round-trip", "[rs]")
{
    // Encode a message.
    using Encoder = rs11::ReedSolomonEncoder<255, 251>;
    Encoder encoder;
    std::vector<std::uint8_t> data({0x48, 0x65, 0x6c, 0x6c, 0x6f, 0x57, 0x6f, 0x72, 0x6c, 0x64});
    encoder.encode(data.data(), data.size());

    // Append the parity symbols.
    std::vector<std::uint8_t> encodedData = data;
    std::copy(encoder.begin(), encoder.end(), std::back_inserter(encodedData));

    REQUIRE(encodedData.size() == data.size() + 4);

    // Introduce some errors.
    auto defectiveData = encodedData;
    for (int i = 2; i <= 3; ++i)
        defectiveData[i] = 'a';

    REQUIRE(defectiveData != encodedData);

    // Decode the message.
    Encoder::decoder_type decoder;
    for (auto d : defectiveData)
        decoder << d;
    decoder.finish();
    REQUIRE(decoder.result() == rs11::DecoderResult::Correctable);

    // Correct the message.
    auto correctedData = defectiveData;
    for (auto e : decoder.errors())
        correctedData[e.index()] = e.correct(defectiveData[e.index()]);

    REQUIRE(correctedData.size() == data.size() + 4);

    for (std::size_t index = 0; index < data.size(); ++index)
    {
        REQUIRE(correctedData[index] == data[index]);
    }
}

TEST_CASE("RS round-trip - gsl", "[rs]")
{
    // Encode a message.
    using Encoder = rs11::ReedSolomonEncoder<255, 251>;
    Encoder encoder;
    std::vector<std::uint8_t> data({ 0x48, 0x65, 0x6c, 0x6c, 0x6f, 0x57, 0x6f, 0x72, 0x6c, 0x64 });
    gsl::span<const std::uint8_t> d = { data };
    encoder.encode(d);

    // Append the parity symbols.
    std::vector<std::uint8_t> encodedData = data;
    std::copy(encoder.begin(), encoder.end(), std::back_inserter(encodedData));

    REQUIRE(encodedData.size() == data.size() + 4);

    // Introduce some errors.
    auto defectiveData = encodedData;
    for (int i = 2; i <= 3; ++i)
        defectiveData[i] = 'a';

    REQUIRE(defectiveData != encodedData);

    // Decode the message.
    Encoder::decoder_type decoder;
    gsl::span<const std::uint8_t> dd = { defectiveData };
    decoder.decodePart(dd);
    decoder.finish();
    REQUIRE(decoder.result() == rs11::DecoderResult::Correctable);

    // Correct the message.
    auto correctedData = defectiveData;
    for (auto e : decoder.errors())
        correctedData[e.index()] = e.correct(defectiveData[e.index()]);

    REQUIRE(correctedData.size() == data.size() + 4);

    for (std::size_t index = 0; index < data.size(); ++index)
    {
        REQUIRE(correctedData[index] == data[index]);
    }
}
