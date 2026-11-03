// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2010-2015 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


/**
 *  \brief Test cases for bcfhelpers
 *
 *
 * \file test_bcfhelpers.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <limits>

#include "helpers/BCFHelpers.hh"
#include <htslib/bgzf.h>

using namespace bcfhelpers;

BOOST_AUTO_TEST_CASE(testClassifyAlleleString)
{
    BOOST_CHECK_EQUAL(classifyAlleleString("").first, AlleleType::MISSING);
    BOOST_CHECK_EQUAL(classifyAlleleString("ACgT").first, AlleleType::NUC);
    BOOST_CHECK_EQUAL(classifyAlleleString("*").first, AlleleType::SYMBOLIC_DEL);
    BOOST_CHECK_EQUAL(classifyAlleleString("<DEL>").first, AlleleType::SYMBOLIC_DEL);
    BOOST_CHECK_EQUAL(classifyAlleleString("<M>").first, AlleleType::SYMBOLIC_OTHER);
    BOOST_CHECK_EQUAL(classifyAlleleString("aj0213").first, AlleleType::UNKNOWN);
    BOOST_CHECK_EQUAL(classifyAlleleString(".").first, AlleleType::MISSING);
}

#include "helpers/Roc.hh"
#include "Variant.hh"

BOOST_AUTO_TEST_CASE(rocDenominatorBoundaries)
{
    roc::Level level;
    BOOST_CHECK_EQUAL(level.precision(), 0);
    BOOST_CHECK_EQUAL(level.recall(), 0);
    BOOST_CHECK_EQUAL(level.fScore(), 0);
    level.unk(8);
    BOOST_CHECK_EQUAL(level.precision(), 0);
    BOOST_CHECK_EQUAL(level.fScore(), 0);
    BOOST_CHECK_EQUAL(level.na(), 1);
    level.fn(3);
    level.fp(2);
    BOOST_CHECK_EQUAL(level.fScore(), 0);
    level.tp(3);
    level.tp2(6);
    BOOST_CHECK_CLOSE(level.recall(), 0.5, 0.0001);
    BOOST_CHECK_CLOSE(level.precision(), 0.75, 0.0001);
    BOOST_CHECK_CLOSE(level.fScore(), 0.6, 0.0001);
}

BOOST_AUTO_TEST_CASE(referenceIdentifierSuffixes)
{
    const char *tmp = std::getenv("TMPDIR");
    const boost::filesystem::path base = tmp ? tmp : ".";
    const auto dir = base / boost::filesystem::unique_path("reference-name-%%%%-%%%%");
    boost::filesystem::create_directories(dir);
    for (const std::string suffix : {".fa", ".fasta", ".fa.gz", ".fasta.gz", "", ".unknown"})
    {
        const auto fasta = (dir / ("T2T-CHM13v2.0" + suffix)).string();
        const std::string sequence = ">custom\nACGT\n";
        if (suffix.size() >= 3 && suffix.compare(suffix.size() - 3, 3, ".gz") == 0)
        {
            BGZF *out = bgzf_open(fasta.c_str(), "w");
            BOOST_REQUIRE(out);
            BOOST_REQUIRE_EQUAL(bgzf_write(out, sequence.data(), sequence.size()), sequence.size());
            BOOST_REQUIRE_EQUAL(bgzf_close(out), 0);
        }
        else { std::ofstream out(fasta); out << sequence; }
        auto header = ph(bcf_hdr_init("w"));
        bcfHeaderFromReference(header.get(), fasta.c_str());
        auto record = bcf_hdr_get_hrec(header.get(), BCF_HL_GEN, "reference", nullptr, nullptr);
        BOOST_REQUIRE(record);
        const std::string expected = "T2T-CHM13v2.0" + (suffix == ".unknown" ? suffix : "");
        BOOST_CHECK_EQUAL(record->value, expected);
        BOOST_CHECK_EQUAL(bcf_hdr_name2id(header.get(), "custom"), 0);
    }
    boost::filesystem::remove_all(dir);
}

BOOST_AUTO_TEST_CASE(referenceDictionaryRoundtrip)
{
    const char *tmp = std::getenv("TMPDIR");
    const boost::filesystem::path base = tmp ? tmp : ".";
    const auto dir = base / boost::filesystem::unique_path("reference-test-%%%%-%%%%");
    boost::filesystem::create_directories(dir);
    struct Cleanup {
        boost::filesystem::path path;
        ~Cleanup() {
            // NFS can retain .nfs files while the process-wide FASTA cache maps them.
            boost::system::error_code ec;
            boost::filesystem::remove_all(path, ec);
            if (ec) std::cerr << "Deferred test cleanup " << path << ": " << ec.message() << "\n";
        }
    } cleanup{dir};
    const std::string fasta = (dir / "custom.fa").string();
    { std::ofstream out(fasta); out << ">zeta\nACGT\n>alpha\nACGTAC\n"; }
    for (const std::string extension : {".vcf", ".bcf"})
    {
        const auto output = (dir / ("out" + extension)).string();
        {
            variant::VariantWriter writer(output.c_str(), fasta.c_str());
            writer.addHeader("##fileformat=VCFv4.1");
            writer.addHeader("##fileformat=VCFv4.2");
            writer.addHeader("##reference=hg19");
            writer.addHeader("##contig=<ID=alpha,length=999>");
            writer.addHeader("##contig=<ID=not_in_reference,length=10>");
            writer.addSample("haploid");
            writer.addSample("diploid");
            writer.addSample("missing");
            writer.addSample("haploid_missing");
            variant::Variants record;
            record.chr = "alpha";
            record.pos = 0;
            record.len = 1;
            record.calls.resize(4);
            record.calls[3].ngt = 1;
            record.calls[3].gt[0] = -1;
            record.calls[0].ngt = 1;
            record.calls[0].gt[0] = 0;
            record.calls[1].ngt = 2;
            record.calls[1].gt[0] = record.calls[1].gt[1] = 0;
            record.calls[1].phased = true;
            writer.put(record);
            record.calls[0].ngt = MAX_GT + 1;
            BOOST_CHECK_THROW(writer.put(record), importexception);
            record.calls[0].ngt = 1;
            record.chr = "not_in_reference";
            BOOST_CHECK_THROW(writer.put(record), std::exception);
        }
        htsFile *input = bcf_open(output.c_str(), "r");
        BOOST_REQUIRE(input);
        auto hdr = ph(bcf_hdr_read(input));
        BOOST_REQUIRE(hdr);
        BOOST_CHECK_EQUAL(bcf_hdr_get_version(hdr.get()), "VCFv4.1");
        int versions = 0;
        for (int i = 0; i < hdr->nhrec; ++i)
            versions += std::string(hdr->hrec[i]->key) == "fileformat";
        BOOST_CHECK_EQUAL(versions, 1);
        BOOST_CHECK_EQUAL(bcf_hdr_nsamples(hdr.get()), 4);
        BOOST_REQUIRE(bcf_hdr_get_hrec(hdr.get(), BCF_HL_FMT, "ID", "GT", nullptr));
        BOOST_CHECK_EQUAL(bcf_hdr_name2id(hdr.get(), "zeta"), 0);
        BOOST_CHECK_EQUAL(bcf_hdr_name2id(hdr.get(), "alpha"), 1);
        BOOST_CHECK_EQUAL(bcf_hdr_name2id(hdr.get(), "chr1"), -1);
        BOOST_CHECK_EQUAL(bcf_hdr_name2id(hdr.get(), "not_in_reference"), -1);
        auto contig = bcf_hdr_get_hrec(hdr.get(), BCF_HL_CTG, "ID", "alpha", nullptr);
        BOOST_REQUIRE(contig);
        const int key = bcf_hrec_find_key(contig, "length");
        BOOST_REQUIRE(key >= 0);
        BOOST_CHECK_EQUAL(contig->vals[key], "6");
        auto rec = pb(bcf_init());
        BOOST_REQUIRE_EQUAL(bcf_read(input, hdr.get(), rec.get()), 0);
        BOOST_CHECK_EQUAL(getChrom(hdr.get(), rec.get()), "alpha");
        BOOST_CHECK_EQUAL(rec->rid, 1);
        int gt[MAX_GT], count;
        bool phased;
        getGT(hdr.get(), rec.get(), 0, gt, count, phased);
        BOOST_CHECK_EQUAL(count, 1);
        BOOST_CHECK_EQUAL(gt[0], 0);
        BOOST_CHECK_EQUAL(gt[1], -1);
        getGT(hdr.get(), rec.get(), 1, gt, count, phased);
        BOOST_CHECK_EQUAL(count, 2);
        BOOST_CHECK(phased);
        getGT(hdr.get(), rec.get(), 2, gt, count, phased);
        BOOST_CHECK_EQUAL(count, 2);
        BOOST_CHECK_EQUAL(gt[0], -1);
        BOOST_CHECK(!phased);
        BOOST_CHECK_EQUAL(gt[1], -1);
        getGT(hdr.get(), rec.get(), 3, gt, count, phased);
        BOOST_CHECK_EQUAL(count, 1);
        BOOST_CHECK_EQUAL(gt[0], -1);
        int32_t *encoded = nullptr;
        int capacity = 0;
        BOOST_REQUIRE_EQUAL(bcf_get_genotypes(hdr.get(), rec.get(), &encoded, &capacity), 8);
        BOOST_CHECK_EQUAL(encoded[1], bcf_int32_vector_end);
        BOOST_CHECK(bcf_gt_is_phased(encoded[3]));
        BOOST_CHECK(bcf_gt_is_missing(encoded[4]));
        BOOST_CHECK(bcf_gt_is_missing(encoded[5]));
        BOOST_CHECK(bcf_gt_is_missing(encoded[6]));
        BOOST_CHECK_EQUAL(encoded[7], bcf_int32_vector_end);
        free(encoded);
        bcf_close(input);
    }
    auto hdr = ph(bcf_hdr_init("w"));
    BOOST_CHECK_THROW(bcfHeaderFromReference(hdr.get(), ""), importexception);
    BOOST_CHECK_THROW(bcfHeaderFromReference(hdr.get(), (dir / "absent.fa").string().c_str()), importexception);
}

BOOST_AUTO_TEST_CASE(rocEqualQualityInclusiveThresholds)
{
    using roc::DecisionType;
    const std::vector<roc::Observation> observations = {
        {0, DecisionType::FN, 1, roc::OBS_FLAG_GM},
        {10, DecisionType::FP, 1, roc::OBS_FLAG_GM},
        {2, DecisionType::UNK, 1, roc::OBS_FLAG_GM},
        {40, DecisionType::TP, 1, roc::OBS_FLAG_GM},
        {40, DecisionType::TP2, 2, roc::OBS_FLAG_GM},
        {40, DecisionType::FP, 3, roc::OBS_FLAG_AM}
    };
    for (bool reverse : {false, true})
    {
        roc::Roc curve;
        for (size_t i = 0; i < observations.size(); ++i)
            curve.add(observations[reverse ? observations.size() - i - 1 : i]);
        for (uint64_t mask : {uint64_t(0), roc::OBS_FLAG_GM})
        {
            std::vector<roc::Level> levels;
            curve.getLevels(levels, 0, mask);
            BOOST_REQUIRE_EQUAL(levels.size(), 4);
            BOOST_CHECK_EQUAL(levels[0].level, 0);
            BOOST_CHECK_EQUAL(levels[1].level, 2);
            BOOST_CHECK_EQUAL(levels[2].level, 10);
            BOOST_CHECK_EQUAL(levels[3].level, 40);
            for (const auto &level : levels)
            {
                BOOST_CHECK_EQUAL(level.tp(), 1);
                BOOST_CHECK_EQUAL(level.tp2(), 2); // Representations remain independent.
                BOOST_CHECK_EQUAL(level.fn(), 1);
                BOOST_CHECK_EQUAL(level.fp(), (mask ? 0 : 3) + (level.level <= 10 ? 1 : 0));
                BOOST_CHECK_EQUAL(level.unk(), level.level <= 2 ? 1 : 0);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(infoFloatWithIntegerStorage)
{
    bcf_hdr_t *header = bcf_hdr_init("w");
    BOOST_REQUIRE(header);
    BOOST_REQUIRE_EQUAL(bcf_hdr_append(header, "##contig=<ID=chrT,length=100>"), 0);
    BOOST_REQUIRE_EQUAL(bcf_hdr_append(header, "##INFO=<ID=SB,Number=.,Type=Float,Description=\"Strand bias\">"), 0);
    BOOST_REQUIRE_EQUAL(bcf_hdr_sync(header), 0);
    bcf1_t *record = bcf_init();
    BOOST_REQUIRE(record);
    record->rid = 0;
    record->pos = 9;
    const int32_t values[] = {-958, -3, 40000};
    BOOST_REQUIRE_EQUAL(bcf_update_info_int32(header, record, "SB", values, 3), 0);
    const auto floats = getInfoFloats(header, record, "SB");
    BOOST_REQUIRE_EQUAL(floats.size(), 3);
    for (size_t i = 0; i < floats.size(); ++i)
        BOOST_CHECK_EQUAL(floats[i], values[i]);
    const auto integers = getInfoInts(header, record, "SB");
    BOOST_REQUIRE_EQUAL(integers.size(), 3);
    for (size_t i = 0; i < integers.size(); ++i)
        BOOST_CHECK_EQUAL(integers[i], values[i]);
    for (const int32_t value : {-3, -958, 40000})
    {
        BOOST_REQUIRE_EQUAL(bcf_update_info_int32(header, record, "SB", &value, 1), 0);
        const bcf_info_t *stored = bcf_get_info(header, record, "SB");
        BOOST_REQUIRE(stored);
        const int expected_type = value == -3 ? BCF_BT_INT8 : value == -958 ? BCF_BT_INT16 : BCF_BT_INT32;
        BOOST_CHECK_EQUAL(stored->type, expected_type);
        BOOST_REQUIRE_EQUAL(getInfoFloats(header, record, "SB").size(), 1);
        BOOST_CHECK_EQUAL(getInfoFloats(header, record, "SB")[0], value);
    }
    const int32_t special[] = {7, bcf_int32_missing, bcf_int32_vector_end};
    BOOST_REQUIRE_EQUAL(bcf_update_info_int32(header, record, "SB", special, 3), 0);
    const auto missing = getInfoFloats(header, record, "SB");
    BOOST_REQUIRE_EQUAL(missing.size(), 2);
    BOOST_CHECK_EQUAL(missing[0], 7);
    BOOST_CHECK(std::isnan(missing[1]));
    BOOST_CHECK_EQUAL(getInfoInts(header, record, "SB")[1], bcf_int32_missing);
    float native_missing;
    bcf_float_set_missing(native_missing);
    BOOST_REQUIRE_EQUAL(bcf_update_info_float(header, record, "SB", &native_missing, 1), 0);
    BOOST_CHECK_THROW(getInfoInts(header, record, "SB"), importexception);
    const float invalid[] = {std::numeric_limits<float>::infinity(),
                             -std::numeric_limits<float>::infinity(),
                             std::numeric_limits<float>::quiet_NaN(), 1e30f};
    for (const float value : invalid)
    {
        BOOST_REQUIRE_EQUAL(bcf_update_info_float(header, record, "SB", &value, 1), 0);
        BOOST_CHECK_THROW(getInfoInts(header, record, "SB"), importexception);
    }
    bcf_destroy(record);
    bcf_hdr_destroy(header);
}
