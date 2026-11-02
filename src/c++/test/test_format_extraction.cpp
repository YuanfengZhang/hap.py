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
 * \brief
 *
 * \file test_format_extraction.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#define BOOST_TEST_NO_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/filesystem/path.hpp>

#include <iostream>
#include <sstream>
#include <cmath>

#include "Variant.hh"
#include "helpers/BCFHelpers.hh"

using namespace variant;


BOOST_AUTO_TEST_CASE(variantReadingFormats)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("test_formats.vcf.gz");

    std::cout << "Reading " << tp << std::endl;

    VariantReader r;
    r.addSample(tp.string().c_str(), "NA12877");

    size_t count = 0;

    const char * expected[] = {
        "chr1:826157-826160 826157-826160:C 0/0 PASS",
        "chr1:826160-826160 826160-826160:CTTTGAT 0/1 153 PASS",
        "chr1:826161-826167 826161-826167:T 0/0 PASS",
        "chr1:826168-826180 826168-826180:T 0/0 PASS",
        "chr1:826180-826182 826180-826182:C 826180-826182:T 0/1 191 PASS"
    };

    float expected_qual[] = {
    	-1,
    	153.0,
    	-1,
    	-1,
		191,
    };

    float expected_gq[] = {
    	138.0,
		193.0,
		102.0,
		96.0,
		231.0,
    };

    float expected_ad[] = {
    	0, 0,
		34, 5,
    	0, 0,
    	0, 0,
		38, 7,
    };

    float expected_ad_r[] = {
    	-1,
	 	34,
    	-1,
    	-1,
		 38,
    };

    float expected_ad_o[] = {
    	0,
	 	0,
    	0,
    	0,
		0,
    };


    int expected_dp[] = {
    	47,
		55,
		35,
		33,
		48,
    };

    while(r.advance())
    {
        Variants & v = r.current();
        std::ostringstream ss;
        ss << v;
        // std::cerr << ss.str() << "\n";
        BOOST_CHECK(count < sizeof(expected)/sizeof(const char *));
        BOOST_CHECK_EQUAL(ss.str(), expected[count]);

        if(expected_qual[count] >= 0)
        {
        	BOOST_CHECK_CLOSE(v.calls[0].qual, expected_qual[count], 0.01);
        }
        else
        {
        	BOOST_CHECK(std::isnan(v.calls[0].qual));
        }

        float gq = v.calls[0].formats["GQX"].asFloat();
		if(!gq) gq = v.calls[0].formats["GQ"].asFloat();
        BOOST_CHECK_CLOSE(gq, expected_gq[count], 0.01);
        BOOST_CHECK_EQUAL(v.calls[0].ad[0], expected_ad[2*count]);
        BOOST_CHECK_EQUAL(v.calls[0].ad[1], expected_ad[2*count+1]);
        BOOST_CHECK_EQUAL(v.calls[0].ad_ref, expected_ad_r[count]);
        BOOST_CHECK_EQUAL(v.calls[0].ad_other, expected_ad_o[count]);
        BOOST_CHECK_EQUAL(v.calls[0].dp, expected_dp[count]);

        ++count;
    }
    BOOST_CHECK(count > 0);
}

BOOST_AUTO_TEST_CASE(variantReadingFormats2)
{
    boost::filesystem::path p(__FILE__);
    boost::filesystem::path tp = p.parent_path()
                                   .parent_path()   // test
                                   .parent_path()   // c++
                                    / boost::filesystem::path("data")
                                    / boost::filesystem::path("test_formats.vcf.gz");

    std::cout << "Reading " << tp << std::endl;

    VariantReader r;
    r.addSample(tp.string().c_str(), "NA12877");
    r.addSample(tp.string().c_str(), "NA12878");

    size_t count = 0;

    const char * expected[] = {
        "chr1:826157-826160 826157-826160:C 0/0 PASS 1/0 PASS",
        "chr1:826160-826160 826160-826160:CTTTGAT 0/1 153 PASS 1/0 153 PASS",
        "chr1:826161-826167 826161-826167:T 0/0 PASS 1/0 PASS",
        "chr1:826168-826180 826168-826180:T 0/0 PASS 1/0 PASS",
        "chr1:826180-826182 826180-826182:C 826180-826182:T 0/1 191 PASS 1/0 191 PASS"
    };

    float expected_qual[] = {
    	-1,
    	153.0,
    	-1,
    	-1,
		191,
    	-1,
    	153.0,
    	-1,
    	-1,
		191,
    };

    float expected_gq[] = {
    	138.0,
		193.0,
		102.0,
		96.0,
		231.0,
    	139.0,
		194.0,
		103.0,
		97.0,
		232.0,
    };

    float expected_ad[] = {
    	0, 0,
		34, 5,
    	0, 0,
    	0, 0,
		38, 7,
    	0, 0,
		5, 34,		// GT is the other way around in this one
    	0, 0,
    	0, 0,
		8, 39,
    };

    float expected_ad_r[] = {
    	-1,
	 	34,
    	-1,
    	-1,
		 38,
    	-1,
	 	34,
    	-1,
    	-1,
		 39,
    };

    float expected_ad_o[] = {
    	0,
	 	0,
    	0,
    	0,
		0,
    	0,
	 	0,
    	0,
    	0,
		1,
    };

    int expected_dp[] = {
    	47,
		55,
		35,
		33,
		48,
    	47,
		55,
		36,
		37,
		49,
    };

    while(r.advance())
    {
        Variants & v = r.current();
        std::ostringstream ss;
        ss << v;
        // std::cerr << ss.str() << "\n";
        BOOST_CHECK(count < sizeof(expected)/sizeof(const char *));
        BOOST_CHECK_EQUAL(ss.str(), expected[count]);

        for (int i = 0; i < 2; ++i)
        {
        	// std::cerr << "\t" << i << "\n";
	        if(expected_qual[count] >= 0)
	        {
	        	BOOST_CHECK_CLOSE(v.calls[i].qual, expected_qual[count + 5*i], 0.01);
	        }
	        else
	        {
	        	BOOST_CHECK(std::isnan(v.calls[i].qual));
	        }
            float gq = v.calls[i].formats["GQX"].asFloat();
            if(gq <= 0) gq = v.calls[i].formats["GQ"].asFloat();
	        BOOST_CHECK_CLOSE(gq, expected_gq[count + 5*i], 0.01);
	        BOOST_CHECK_EQUAL(v.calls[i].ad[0], expected_ad[2*count + 10*i]);
	        BOOST_CHECK_EQUAL(v.calls[i].ad[1], expected_ad[2*count+1 + 10*i]);
	        BOOST_CHECK_EQUAL(v.calls[i].ad_ref, expected_ad_r[count + 5*i]);
	        BOOST_CHECK_EQUAL(v.calls[i].ad_other, expected_ad_o[count + 5*i]);
	        BOOST_CHECK_EQUAL(v.calls[i].dp, expected_dp[count + 5*i]);
        }
        ++count;
    }
    BOOST_CHECK(count > 0);
}

// Exercise HTSlib's actual packed widths, sample strides and vector terminators.
BOOST_AUTO_TEST_CASE(formatPublicApiBoundaries)
{
    auto hdr = bcfhelpers::ph(bcf_hdr_init("w"));
    bcf_hdr_append(hdr.get(), "##contig=<ID=custom,length=20>");
    bcf_hdr_append(hdr.get(), "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">");
    bcf_hdr_append(hdr.get(), "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">");
    bcf_hdr_append(hdr.get(), "##FORMAT=<ID=F,Number=.,Type=Float,Description=\"F\">");
    bcf_hdr_add_sample(hdr.get(), "truth");
    bcf_hdr_add_sample(hdr.get(), "query");
    BOOST_REQUIRE_EQUAL(bcf_hdr_sync(hdr.get()), 0);
    auto rec = bcfhelpers::pb(bcf_init());
    rec->rid = 0;
    BOOST_REQUIRE_EQUAL(bcf_update_alleles_str(hdr.get(), rec.get(), "A,C"), 0);
    for (int value : {10, 1000, 100000})
    {
        int32_t ad[] = {value, value + 1, value + 2, value + 3};
        BOOST_REQUIRE_EQUAL(bcf_update_format_int32(hdr.get(), rec.get(), "AD", ad, 4), 0);
        int guarded[] = {12345, -1, -1, 54321};
        bcfhelpers::getAD(hdr.get(), rec.get(), 1, guarded + 1, 2);
        BOOST_CHECK_EQUAL(guarded[0], 12345);
        BOOST_CHECK_EQUAL(guarded[1], value + 2);
        BOOST_CHECK_EQUAL(guarded[2], value + 3);
        BOOST_CHECK_EQUAL(guarded[3], 54321);
        BOOST_CHECK_THROW(bcfhelpers::getAD(hdr.get(), rec.get(), 1, guarded + 1, 1), bcfhelpers::importexception);
        BOOST_CHECK_EQUAL(guarded[1], value + 2);
        BOOST_CHECK_EQUAL(guarded[2], value + 3);
        BOOST_CHECK_THROW(bcfhelpers::getAD(hdr.get(), rec.get(), 1, nullptr, 0), bcfhelpers::importexception);
        BOOST_CHECK_EQUAL(bcfhelpers::getFormatString(hdr.get(), rec.get(), "AD", 1), std::to_string(value + 2));
    }
    int32_t padded[] = {4, 5, bcf_int32_missing, bcf_int32_vector_end};
    bcf_update_format_int32(hdr.get(), rec.get(), "AD", padded, 4);
    const auto missing = bcfhelpers::getFormatInts(hdr.get(), rec.get(), "AD", 1);
    BOOST_REQUIRE_EQUAL(missing.size(), 1);
    BOOST_CHECK_EQUAL(missing[0], 0); // Preserve the existing numeric missing-value convention.
    BOOST_CHECK(bcfhelpers::getFormatInts(hdr.get(), rec.get(), "AD", -1).empty());
    BOOST_CHECK(bcfhelpers::getFormatInts(hdr.get(), rec.get(), "AD", 2).empty());
    float floats[] = {1.25f, 2.5f, 3.75f, 0};
    bcf_float_set_vector_end(floats[3]);
    bcf_update_format_float(hdr.get(), rec.get(), "F", floats, 4);
    auto values = bcfhelpers::getFormatFloats(hdr.get(), rec.get(), "F", 1);
    BOOST_REQUIRE_EQUAL(values.size(), 1);
    BOOST_CHECK_CLOSE(values[0], 3.75f, 0.001);
    bcf_float_set_missing(floats[2]);
    bcf_update_format_float(hdr.get(), rec.get(), "F", floats, 4);
    values = bcfhelpers::getFormatFloats(hdr.get(), rec.get(), "F", 1);
    BOOST_REQUIRE_EQUAL(values.size(), 1);
    BOOST_CHECK(std::isnan(values[0]));
    bcf_update_format_int32(hdr.get(), rec.get(), "AD", nullptr, 0);
    int untouched = 99;
    bcfhelpers::getAD(hdr.get(), rec.get(), 1, &untouched, 1);
    BOOST_CHECK_EQUAL(untouched, 99);
}

BOOST_AUTO_TEST_CASE(genotypePublicApiBoundaries)
{
    auto hdr = bcfhelpers::ph(bcf_hdr_init("w"));
    bcf_hdr_append(hdr.get(), "##contig=<ID=custom,length=20>");
    bcf_hdr_append(hdr.get(), "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">");
    bcf_hdr_add_sample(hdr.get(), "truth");
    bcf_hdr_add_sample(hdr.get(), "query");
    BOOST_REQUIRE_EQUAL(bcf_hdr_sync(hdr.get()), 0);
    auto rec = bcfhelpers::pb(bcf_init());
    rec->rid = 0;
    bcf_update_alleles_str(hdr.get(), rec.get(), "A,C");
    int32_t encoded[] = {bcf_gt_unphased(0), bcf_gt_phased(1), bcf_gt_unphased(1), bcf_int32_vector_end};
    bcf_update_genotypes(hdr.get(), rec.get(), encoded, 4);
    int gt[MAX_GT], count = -1;
    bool phased = false;
    bcfhelpers::getGT(hdr.get(), rec.get(), 0, gt, count, phased);
    BOOST_CHECK_EQUAL(count, 2);
    BOOST_CHECK_EQUAL(gt[0], 0);
    BOOST_CHECK_EQUAL(gt[1], 1);
    BOOST_CHECK(phased);
    bcfhelpers::getGT(hdr.get(), rec.get(), 1, gt, count, phased);
    BOOST_CHECK_EQUAL(count, 1);
    BOOST_CHECK_EQUAL(gt[0], 1);
    BOOST_CHECK_EQUAL(gt[1], -1);
    BOOST_CHECK(!phased);
    encoded[2] = bcf_gt_missing;
    bcf_update_genotypes(hdr.get(), rec.get(), encoded, 4);
    bcfhelpers::getGT(hdr.get(), rec.get(), 1, gt, count, phased);
    BOOST_CHECK_EQUAL(count, 1);
    BOOST_CHECK_EQUAL(gt[0], -1);
    encoded[2] = bcf_gt_unphased(2);
    bcf_update_genotypes(hdr.get(), rec.get(), encoded, 4);
    BOOST_CHECK_THROW(bcfhelpers::getGT(hdr.get(), rec.get(), 1, gt, count, phased), bcfhelpers::importexception);
}

#include "BlockQuantify.hh"

BOOST_AUTO_TEST_CASE(ga4ghIndependentSampleDecisions)
{
    boost::filesystem::path source(__FILE__);
    const auto reference = source.parent_path().parent_path().parent_path() / "data" / "microhg19.fa";
    FastaFile fasta(reference.string().c_str());
    auto hdr = bcfhelpers::ph(bcf_hdr_init("w"));
    bcfhelpers::bcfHeaderFromReference(hdr.get(), reference.string().c_str());
    bcf_hdr_append(hdr.get(), "##INFO=<ID=Regions,Number=.,Type=String,Description=\"Regions\">");
    bcf_hdr_append(hdr.get(), "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Preserved BD\">");
    bcf_hdr_append(hdr.get(), "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Preserved BK\">");
    bcf_hdr_append(hdr.get(), "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Preserved BI\">");
    bcf_hdr_append(hdr.get(), "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Preserved QQ\">");
    bcf_hdr_add_sample(hdr.get(), "TRUTH");
    bcf_hdr_add_sample(hdr.get(), "QUERY");
    BOOST_REQUIRE_EQUAL(bcf_hdr_sync(hdr.get()), 0);
    auto quantify = makeQuantifier(hdr.get(), fasta, "ga4gh", "");
    quantify->updateHeader(hdr.get());
    auto bd = bcf_hdr_get_hrec(hdr.get(), BCF_HL_FMT, "ID", "BD", nullptr);
    BOOST_REQUIRE(bd);
    BOOST_CHECK_EQUAL(bd->vals[bcf_hrec_find_key(bd, "Description")], "\"Preserved BD\"");
    auto rec = bcfhelpers::pb(bcf_init());
    rec->rid = 0;
    rec->pos = 0;
    bcf_update_alleles_str(hdr.get(), rec.get(), "A,C");
    int32_t gt[] = {bcf_gt_unphased(0), bcf_gt_unphased(1), bcf_gt_unphased(0), bcf_gt_unphased(1)};
    bcf_update_genotypes(hdr.get(), rec.get(), gt, 4);
    bcfhelpers::setFormatStrings(hdr.get(), rec.get(), "BD", {"FN", "FP"});
    bcfhelpers::setFormatStrings(hdr.get(), rec.get(), "BK", {"gm", "gm"});
    bcfhelpers::setFormatFloats(hdr.get(), rec.get(), "QQ", {10, 20});
    bcf_update_info_string(hdr.get(), rec.get(), "Regions", "CONF");
    quantify->add(bcf_dup(rec.get()));
    quantify->count();
    const auto output = quantify->getVariants().front();
    BOOST_CHECK_EQUAL(bcfhelpers::getFormatString(hdr.get(), output, "BD", 0), "FN");
    BOOST_CHECK_EQUAL(bcfhelpers::getFormatString(hdr.get(), output, "BD", 1), "FP");
    BOOST_CHECK_EQUAL(bcfhelpers::getFormatString(hdr.get(), output, "BVT", 0), "SNP");
    BOOST_CHECK_EQUAL(bcfhelpers::getFormatString(hdr.get(), output, "BVT", 1), "SNP");
    const auto &rocs = quantify->getRocs();
    const auto all = rocs.find("a:SNP:ALL");
    BOOST_REQUIRE(all != rocs.end());
    const auto counts = all->second.getTotals();
    BOOST_CHECK_EQUAL(counts.fn(), 1);
    BOOST_CHECK_EQUAL(counts.fp(), 1);
    BOOST_CHECK_EQUAL(counts.tp(), 0);
    BOOST_CHECK_EQUAL(counts.tp2(), 0);
}

BOOST_AUTO_TEST_CASE(benchmarkSuperlocusMissingQuality)
{
    const auto reference = boost::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data" / "microhg19.fa";
    FastaFile fasta(reference.string().c_str());
    auto hdr = bcfhelpers::ph(bcf_hdr_init("w"));
    bcfhelpers::bcfHeaderFromReference(hdr.get(), reference.string().c_str());
    bcf_hdr_append(hdr.get(), "##INFO=<ID=BS,Number=1,Type=Integer,Description=\"Superlocus\">");
    bcf_hdr_append(hdr.get(), "##INFO=<ID=Regions,Number=.,Type=String,Description=\"Regions\">");
    bcf_hdr_add_sample(hdr.get(), "TRUTH");
    bcf_hdr_add_sample(hdr.get(), "QUERY");
    bcf_hdr_sync(hdr.get());
    auto quantify = makeQuantifier(hdr.get(), fasta, "ga4gh", "");
    quantify->updateHeader(hdr.get());
    for (int i = 0; i < 3; ++i)
    {
        auto rec = bcfhelpers::pb(bcf_init());
        rec->rid = 0;
        rec->pos = i;
        bcf_update_alleles_str(hdr.get(), rec.get(), "A,C");
        int32_t gt[] = {bcf_gt_unphased(0), bcf_gt_unphased(1), bcf_gt_unphased(0), bcf_gt_unphased(1)};
        if (i == 2) gt[2] = gt[3] = bcf_gt_missing;
        bcf_update_genotypes(hdr.get(), rec.get(), gt, 4);
        bcfhelpers::setFormatStrings(hdr.get(), rec.get(), "BD", {"TP", i == 2 ? "." : "TP"});
        bcfhelpers::setFormatStrings(hdr.get(), rec.get(), "BK", {"gm", "gm"});
        float qq[] = {bcfhelpers::missing_float(), 40};
        if (i != 0) bcf_float_set_missing(qq[1]);
        bcf_update_format_float(hdr.get(), rec.get(), "QQ", qq, 2);
        const int32_t bs = 1;
        bcf_update_info_int32(hdr.get(), rec.get(), "BS", &bs, 1);
        bcf_update_info_string(hdr.get(), rec.get(), "Regions", "CONF");
        quantify->add(bcf_dup(rec.get()));
    }
    quantify->count();
    for (const auto rec : quantify->getVariants())
        BOOST_CHECK_EQUAL(bcfhelpers::getFormatFloat(hdr.get(), rec, "QQ", 0), 40);
    BOOST_CHECK(std::isnan(bcfhelpers::getFormatFloat(hdr.get(), quantify->getVariants().back(), "QQ", 1)));
}

BOOST_AUTO_TEST_CASE(readerSeekPublicApiBoundaries)
{
    boost::filesystem::path source(__FILE__);
    const auto data = source.parent_path().parent_path().parent_path() / "data";
    VariantReader reader;
    reader.addSample((data / "test_formats.vcf.gz").string().c_str(), "NA12877");
    reader.rewind("chr1", 826157);
    BOOST_CHECK(reader.advance());
    reader.rewind("chr1", 2000000000LL);
    BOOST_CHECK(!reader.advance());
    reader.rewind("absent_custom_contig", 0);
    BOOST_CHECK(!reader.advance());
    reader.rewind("chr1", 826157);
    BOOST_CHECK(reader.advance());
    VariantReader unindexed;
    BOOST_CHECK_THROW(unindexed.addSample((data / "test_formats.vcf").string().c_str(), "NA12877"), std::exception);
}

BOOST_AUTO_TEST_CASE(readerMultipleFilesAndLongCoordinates)
{
    const char *tmp = std::getenv("TMPDIR");
    const std::string base = std::string(tmp ? tmp : ".") + "/multi-reader";
    const int64_t boundary = int64_t(1) << 31;
    const auto write_file = [&](const std::string &filename, const char *chrom, int64_t pos) {
        auto hdr = bcfhelpers::ph(bcf_hdr_init("w"));
        const std::string contig = std::string("##contig=<ID=") + chrom + ",length=4000000000>";
        bcf_hdr_append(hdr.get(), contig.c_str());
        bcf_hdr_append(hdr.get(), "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">");
        bcf_hdr_add_sample(hdr.get(), "sample");
        bcf_hdr_sync(hdr.get());
        auto rec = bcfhelpers::pb(bcf_init());
        rec->rid = 0;
        rec->pos = pos;
        bcf_update_alleles_str(hdr.get(), rec.get(), "A,C");
        int32_t gt[] = {bcf_gt_unphased(0), bcf_gt_unphased(1)};
        bcf_update_genotypes(hdr.get(), rec.get(), gt, 2);
        htsFile *out = bcf_open(filename.c_str(), "wz");
        BOOST_REQUIRE(out);
        BOOST_REQUIRE_EQUAL(bcf_hdr_write(out, hdr.get()), 0);
        BOOST_REQUIRE_EQUAL(bcf_write(out, hdr.get(), rec.get()), 0);
        BOOST_REQUIRE_EQUAL(bcf_close(out), 0);
        BOOST_REQUIRE_EQUAL(bcf_index_build(filename.c_str(), 14), 0);
    };
    const std::string first = base + "-first.vcf.gz", second = base + "-second.vcf.gz";
    write_file(first, "first", 10);
    write_file(second, "second", boundary);
    for (int filter = 0; filter < 3; ++filter)
    {
        VariantReader reader;
        if (filter == 1) reader.setRegions("second:2147483648-2147483650", false);
        if (filter == 2) reader.setTargets("second:2147483648-2147483650", false);
        reader.addSample(first.c_str(), "sample");
        reader.addSample(second.c_str(), "sample");
        reader.rewind("second", boundary - 1);
        BOOST_REQUIRE(reader.advance());
        BOOST_CHECK_EQUAL(reader.current().chr, "second");
        BOOST_CHECK_EQUAL(reader.current().pos, boundary);
        BOOST_CHECK(!reader.advance());
        reader.rewind("second", boundary + 1);
        BOOST_CHECK(!reader.advance());
        reader.rewind("absent", 0);
        BOOST_CHECK(!reader.advance());
        reader.rewind("second", boundary);
        BOOST_REQUIRE(reader.advance());
        BOOST_CHECK_EQUAL(reader.current().pos, boundary);
    }
}

BOOST_AUTO_TEST_CASE(readerDeclaredEmptyContig)
{
    const char *tmp = std::getenv("TMPDIR");
    const std::string filename = std::string(tmp ? tmp : ".") + "/declared-empty.vcf.gz";
    auto hdr = bcfhelpers::ph(bcf_hdr_init("w"));
    bcf_hdr_append(hdr.get(), "##contig=<ID=empty,length=100>");
    bcf_hdr_append(hdr.get(), "##contig=<ID=present,length=100>");
    bcf_hdr_append(hdr.get(), "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">");
    bcf_hdr_add_sample(hdr.get(), "sample");
    bcf_hdr_sync(hdr.get());
    auto rec = bcfhelpers::pb(bcf_init());
    rec->rid = 1;
    rec->pos = 10;
    bcf_update_alleles_str(hdr.get(), rec.get(), "A,C");
    int32_t gt[] = {bcf_gt_unphased(0), bcf_gt_unphased(1)};
    bcf_update_genotypes(hdr.get(), rec.get(), gt, 2);
    htsFile *out = bcf_open(filename.c_str(), "wz");
    BOOST_REQUIRE(out);
    BOOST_REQUIRE_EQUAL(bcf_hdr_write(out, hdr.get()), 0);
    BOOST_REQUIRE_EQUAL(bcf_write(out, hdr.get(), rec.get()), 0);
    bcf_close(out);
    BOOST_REQUIRE_EQUAL(bcf_index_build(filename.c_str(), 14), 0);
    VariantReader reader;
    reader.addSample(filename.c_str(), "sample");
    reader.rewind("empty", 0);
    BOOST_CHECK(!reader.advance());
    reader.rewind("present", 0);
    BOOST_CHECK(reader.advance());
    reader.rewind("absent", 0);
    BOOST_CHECK(!reader.advance());
    reader.rewind("present", 0);
    BOOST_CHECK(reader.advance());
}
