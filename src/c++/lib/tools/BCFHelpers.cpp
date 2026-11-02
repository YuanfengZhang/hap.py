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
 *  \brief Implementation of BCF file format helpers
 *
 * \file BCFHelpers.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "helpers/BCFHelpers.hh"

#include "Error.hh"

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include <memory>
#include <limits>
#include <set>
#include <boost/algorithm/string.hpp>

/**
 * @brief Helper to get out GT fields
 */

namespace bcfhelpers
{

    namespace _impl
    {
        /** C++-ified version of bcf_get_format_values */
        template <typename target_type_t>
        struct bcf_get_numeric_format
        {
            bcf_get_numeric_format() {}
            /* return true when successful
             */
            void operator()(const bcf_hdr_t *hdr, bcf1_t *line,
                                       const char *tag, int isample,
                                       std::vector<target_type_t> & dest) const
            {
                dest.clear();
                const int samples = bcf_hdr_nsamples(hdr);
                if (isample < 0 || isample >= samples) return;
                const int id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag);
                if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, id)) return;
                const int type = bcf_hdr_id2type(hdr, BCF_HL_FMT, id);
                int capacity = 0;
                if (type == BCF_HT_INT)
                {
                    int32_t *raw = nullptr;
                    const int count = bcf_get_format_int32(hdr, line, tag, &raw, &capacity);
                    std::unique_ptr<int32_t, decltype(&free)> values(raw, &free);
                    if (count <= 0) return;
                    const int stride = count / samples;
                    for (int i = 0; i < stride; ++i)
                    {
                        const int32_t value = raw[isample * stride + i];
                        if (value == bcf_int32_vector_end) break;
                        dest.push_back(value == bcf_int32_missing ? target_type_t(0) : target_type_t(value));
                    }
                }
                else if (type == BCF_HT_REAL)
                {
                    float *raw = nullptr;
                    const int count = bcf_get_format_float(hdr, line, tag, &raw, &capacity);
                    std::unique_ptr<float, decltype(&free)> values(raw, &free);
                    if (count <= 0) return;
                    const int stride = count / samples;
                    for (int i = 0; i < stride; ++i)
                    {
                        const float value = raw[isample * stride + i];
                        if (bcf_float_is_vector_end(value)) break;
                        if (bcf_float_is_missing(value))
                            dest.push_back(std::numeric_limits<target_type_t>::is_integer
                                           ? target_type_t(0)
                                           : std::numeric_limits<target_type_t>::quiet_NaN());
                        else if (std::numeric_limits<target_type_t>::is_integer &&
                                 (!std::isfinite(value) || value < std::numeric_limits<target_type_t>::lowest() ||
                                  static_cast<double>(value) > std::numeric_limits<target_type_t>::max()))
                            throw importexception("Non-finite or out-of-range numeric FORMAT value");
                        else dest.push_back(target_type_t(value));
                    }
                }
            }
        };

        template<typename T>
        std::vector<T> numericInfo(bcf_hdr_t *header, bcf1_t *line, const char *field)
        {
            std::vector<T> result;
            const int id = bcf_hdr_id2int(header, BCF_DT_ID, field);
            if (!bcf_hdr_idinfo_exists(header, BCF_HL_INFO, id)) return result;
            const int type = bcf_hdr_id2type(header, BCF_HL_INFO, id);
            int capacity = 0;
            if (type == BCF_HT_INT)
            {
                int32_t *raw = nullptr;
                const int count = bcf_get_info_int32(header, line, field, &raw, &capacity);
                std::unique_ptr<int32_t, decltype(&free)> values(raw, &free);
                for (int i = 0; i < count; ++i)
                {
                    if (raw[i] == bcf_int32_vector_end) break;
                    result.push_back(static_cast<T>(raw[i]));
                }
            }
            else if (type == BCF_HT_REAL)
            {
                const bcf_info_t *info = bcf_get_info(header, line, field);
                if (info && (info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 ||
                             info->type == BCF_BT_INT32))
                {
                    // A float declaration may use integer BCF storage; convert values through the public decoder, not by reinterpreting bits.
                    const int64_t missing = info->type == BCF_BT_INT8 ? bcf_int8_missing :
                                            info->type == BCF_BT_INT16 ? bcf_int16_missing : bcf_int32_missing;
                    const int64_t end = info->type == BCF_BT_INT8 ? bcf_int8_vector_end :
                                        info->type == BCF_BT_INT16 ? bcf_int16_vector_end : bcf_int32_vector_end;
                    uint8_t *cursor = info->vptr;
                    for (int i = 0; i < info->len; ++i)
                    {
                        const int64_t value = bcf_dec_int1(cursor, info->type, &cursor);
                        if (value == end) break;
                        if (value == missing)
                            result.push_back(std::numeric_limits<T>::is_integer ?
                                             static_cast<T>(bcf_int32_missing) : std::numeric_limits<T>::quiet_NaN());
                        else result.push_back(static_cast<T>(value));
                    }
                    return result;
                }
                float *raw = nullptr;
                const int count = bcf_get_info_float(header, line, field, &raw, &capacity);
                std::unique_ptr<float, decltype(&free)> values(raw, &free);
                for (int i = 0; i < count; ++i)
                {
                    if (bcf_float_is_vector_end(raw[i])) break;
                    if (std::numeric_limits<T>::is_integer &&
                        (!std::isfinite(raw[i]) || raw[i] < std::numeric_limits<T>::lowest() ||
                         static_cast<double>(raw[i]) > std::numeric_limits<T>::max()))
                        throw importexception(std::string("Non-finite or out-of-range numeric INFO value: ") + field + " at " + std::to_string(line->pos + 1) + " value " + std::to_string(raw[i]));
                    result.push_back(static_cast<T>(raw[i]));
                }
            }
            else if (type == BCF_HT_STR)
            {
                char *raw = nullptr;
                const int count = bcf_get_info_string(header, line, field, &raw, &capacity);
                std::unique_ptr<char, decltype(&free)> value(raw, &free);
                if (count > 0)
                {
                    const double number = std::strtod(raw, nullptr);
                    if (std::numeric_limits<T>::is_integer &&
                        (!std::isfinite(number) || number < std::numeric_limits<T>::lowest() ||
                         number > std::numeric_limits<T>::max()))
                        throw importexception("Out-of-range string INFO value");
                    result.push_back(static_cast<T>(number));
                }
            }
            return result;
        }
    }

    bool seek(bcf_srs_t *readers, const char *chromosome, int64_t position)
    {
        readers->errnum = static_cast<bcf_sr_error>(0);
        // In 1.24 _reader_seek returns zero on success, despite the header comment.
        // The indexed region dictionary, not just VCF header declarations, determines presence.
        if (chromosome && readers->regions &&
            bcf_sr_regions_seek(readers->regions, chromosome) < 0)
            return false;
        bcf_sr_seek(readers, chromosome, position < 0 ? 0 : position);
        if (readers->errnum)
            throw importexception("Cannot seek VCF/BCF: " + std::string(bcf_sr_strerror(readers->errnum)));
        return true;
    }

    void bcfHeaderFromReference(bcf_hdr_t * header, const char * reference)
    {
        if (!reference || !*reference)
            throw importexception("A readable FASTA reference is required for the output sequence dictionary");
        std::unique_ptr<faidx_t, decltype(&fai_destroy)> index(fai_load(reference), &fai_destroy);
        if (!index) throw importexception("Cannot load FASTA reference/index: " + std::string(reference));
        bcf_hdr_append(header, ("##reference=" + std::string(reference)).c_str());
        for (int i = 0; i < faidx_nseq(index.get()); ++i)
        {
            const char *name = faidx_iseq(index.get(), i);
            const auto length = faidx_seq_len64(index.get(), name);
            const std::string line = "##contig=<ID=" + std::string(name) + ",length=" + std::to_string(length) + ">";
            if (bcf_hdr_append(header, line.c_str()) != 0)
                throw importexception("Cannot append reference contig: " + std::string(name));
        }
        bcf_hdr_append(header, "##INFO=<ID=END,Number=.,Type=Integer,Description=\"SV end position\">");
        bcf_hdr_append(header, "##INFO=<ID=IMPORT_FAIL,Number=.,Type=Flag,Description=\"Flag to identify variants that could not be imported.\">");
        bcf_hdr_append(header, "##FORMAT=<ID=AGT,Number=1,Type=String,Description=\"Genotypes at ambiguous locations\">");
        bcf_hdr_append(header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
        bcf_hdr_append(header, "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">");
        bcf_hdr_append(header, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
        bcf_hdr_append(header, "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depths\">");
        bcf_hdr_append(header, "##FORMAT=<ID=ADO,Number=.,Type=Integer,Description=\"Summed depth of non-called alleles.\">");
    }


    /** extract chrom / pos / length */
    void getLocation(bcf_hdr_t * hdr, bcf1_t * rec, int64_t & refstart, int64_t & refend)
    {
        bcf_unpack(rec, BCF_UN_STR);
        refstart = rec->pos;
        refend = refstart;

        int endfield = getInfoInt(hdr, rec, "END", -1);

        if(endfield > 0)
        {
            // if there is an end field, don't validate the ref allele
            refend = endfield-1;
        }
        else if(rec->d.allele)
        {
            refend = refstart + strlen(rec->d.allele[0]) - 1;
            if (strchr(rec->d.allele[0], '.') || strchr(rec->d.allele[0], '-'))
            {
                // length might be inaccurate now
                refend = refstart;
                throw bcfhelpers::importexception(
                    std::string("[W] Unsupported REF allele with undefined length: ") +
                    std::string(rec->d.allele[0]));
            }
        }
    }

    /**
     * @brief Retrieve an info field as an integer
     *
     * @param result the default to return if the field is not present
     */
    std::string getInfoString(bcf_hdr_t *header, bcf1_t *line, const char *field, const char *def_result)
    {
        const int id = bcf_hdr_id2int(header, BCF_DT_ID, field);
        if (!bcf_hdr_idinfo_exists(header, BCF_HL_INFO, id)) return def_result;
        const int type = bcf_hdr_id2type(header, BCF_HL_INFO, id);
        if (type == BCF_HT_FLAG)
            return getInfoFlag(header, line, field) ? "NULL" : def_result;
        if (type == BCF_HT_INT)
        {
            const auto values = getInfoInts(header, line, field);
            return values.empty() ? def_result : std::to_string(values.front());
        }
        if (type == BCF_HT_REAL)
        {
            const auto values = getInfoFloats(header, line, field);
            if (values.empty()) return def_result;
            char text[256];
            snprintf(text, sizeof(text), "%g", values.front());
            return text;
        }
        char *raw = nullptr;
        int capacity = 0;
        const int count = bcf_get_info_string(header, line, field, &raw, &capacity);
        std::unique_ptr<char, decltype(&free)> value(raw, &free);
        return count >= 0 ? std::string(raw, count) : std::string(def_result);
    }

    int getInfoInt(bcf_hdr_t *header, bcf1_t *line, const char *field, int result)
    {
        const auto values = getInfoInts(header, line, field);
        return values.empty() ? result : values.front();
    }

    std::vector<int> getInfoInts(bcf_hdr_t *header, bcf1_t *line, const char *field)
    {
        return _impl::numericInfo<int>(header, line, field);
    }

    float getInfoFloat(bcf_hdr_t *header, bcf1_t *line, const char *field)
    {
        const auto values = getInfoFloats(header, line, field);
        return values.empty() ? std::numeric_limits<float>::quiet_NaN() : values.front();
    }

    std::vector<float> getInfoFloats(bcf_hdr_t *header, bcf1_t *line, const char *field)
    {
        return _impl::numericInfo<float>(header, line, field);
    }

    /**
     * @brief Retrieve an info flag
     *
     * @return true of the flag is set
     */
    bool getInfoFlag(bcf_hdr_t * hdr, bcf1_t * line, const char * field)
    {
        bcf_unpack(line, BCF_UN_INFO);
        return bcf_get_info_flag(hdr, line, field, nullptr, 0) == 1;
    }

    /**
     * @brief Read the GT field
     */
    void getGT(bcf_hdr_t * header, bcf1_t * line, int isample, int * gt, int & ngt, bool & phased)
    {
        ngt = 0;
        phased = false;
        std::fill(gt, gt + MAX_GT, -1);
        const int samples = bcf_hdr_nsamples(header);
        if (isample < 0 || isample >= samples) return;
        int32_t *raw = nullptr;
        int capacity = 0;
        const int count = bcf_get_genotypes(header, line, &raw, &capacity);
        std::unique_ptr<int32_t, decltype(&free)> values(raw, &free);
        if (count <= 0) return;
        const int stride = count / samples;
        const int32_t *sample = raw + isample * stride;
        while (ngt < stride && sample[ngt] != bcf_int32_vector_end) ++ngt;
        if (ngt > MAX_GT)
        {
            throw importexception("GT ploidy exceeds supported capacity at " + std::string(bcf_hdr_id2name(header, line->rid)) + ":" + std::to_string(line->pos + 1));
        }
        for (int i = 0; i < ngt; ++i)
        {
            if (bcf_gt_is_missing(sample[i])) continue;
            const int allele = bcf_gt_allele(sample[i]);
            if (allele < 0 || allele >= line->n_allele)
                throw importexception("Invalid GT allele index " + std::to_string(allele) + " outside REF/ALT alleles at " + std::string(bcf_hdr_id2name(header, line->rid)) + ":" + std::to_string(line->pos + 1) + " sample " + std::to_string(isample));
            gt[i] = allele;
            phased = phased || (i > 0 && bcf_gt_is_phased(sample[i]));
        }
    }

    /** read GQ(X) -- will use in this order: GQX, GQ, -1
     *
     * This is somewhat Illumina-specific, TODO: refactor to avoid using this function in a general setting
     * */
    void getGQ(const bcf_hdr_t * header, bcf1_t * line, int isample, float & gq)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<float> gf;

        std::vector<float> values;
        gf(header, line, "GQX", isample, values);
        if(values.empty())
        {
            gf(header, line, "GQ", isample, values);
        }
        if(values.size() > 1)
        {
            std::cerr << "[W] too many GQ fields at " << bcf_hdr_id2name(header, line->rid) << ":" << line->pos << "\n";
        }
        if(values.empty())
        {
            gq = -1;
        }
        else
        {
            gq = values[0];
        }
    }

    /** read AD */
    void getAD(const bcf_hdr_t * header, bcf1_t * line, int isample, int *ad, int max_ad)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<int> gf;

        std::vector<int> values;
        gf(header, line, "AD", isample, values);
        if (max_ad < 0 || values.size() > static_cast<size_t>(max_ad) ||
            (!ad && !values.empty()))
        {
            throw importexception("AD field exceeds destination capacity at " + std::string(bcf_hdr_id2name(header, line->rid)) + ":" + std::to_string(line->pos + 1) + " sample " + std::to_string(isample) + ": capacity=" + std::to_string(max_ad) + " values=" + std::to_string(values.size()) + "; destination was not written");
        }
        for(size_t q = 0; q < values.size(); ++q)
        {
            ad[q] = values[q];
        }
    }

    /** read DP(I) -- will use in this order: DP, DPI, -1 */
    void getDP(const bcf_hdr_t * header, bcf1_t * line, int isample, int & dp)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<int> gf;

        std::vector<int> values;
        gf(header, line, "DP", isample, values);
        if(values.empty())
        {
            gf(header, line, "DPI", isample, values);
        }
        if(values.size() > 1)
        {
            std::cerr << "[W] too many DP fields at " << bcf_hdr_id2name(header, line->rid) << ":" << line->pos << "\n";
        }
        if(values.empty())
        {
            dp = 0;
        }
        else
        {
            dp = values[0];
        }
    }

    /** read a format field as a single int. */
    int getFormatInt(const bcf_hdr_t * header, bcf1_t * line, const char * field, int isample, int defaultresult)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<int> gf;

        std::vector<int> values;
        gf(header, line, field, isample, values);
        if(values.size() > 1)
        {
            std::ostringstream os;
            os << "[W] too many " << field << " fields at " << bcf_hdr_id2name(header, line->rid) << ":" << line->pos;
            throw importexception(os.str());
        }
        if(values.size() == 1)
        {
            return values[0];
        }
        return defaultresult;
    }

    std::vector<int> getFormatInts(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<int> gf;
        std::vector<int> result;
        gf(header, line, field, isample, result);
        return result;
    }

    /** read a format field as a single float. default return value is NaN */
    float getFormatFloat(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample)
    {
        using namespace _impl;
        float result = std::numeric_limits<float>::quiet_NaN();
        static const bcf_get_numeric_format<float> gf;

        std::vector<float> values;
        gf(header, line, field, isample, values);
        if(values.size() > 1)
        {
            std::ostringstream os;
            os << "[W] too many " << field << " fields at " << bcf_hdr_id2name(header, line->rid) << ":" << line->pos;
            throw importexception(os.str());
        }
        if(values.size() == 1)
        {
            result = values[0];
        }
        return result;
    }

    std::vector<float> getFormatFloats(const bcf_hdr_t *header, bcf1_t *line, const char *field, int isample)
    {
        using namespace _impl;
        static const bcf_get_numeric_format<float> gf;
        std::vector<float> result;
        gf(header, line, field, isample, result);
        return result;
    }

    /** read a format field as a single double. result will not be overwritten on failure */
    std::string getFormatString(const bcf_hdr_t * hdr, bcf1_t * line, const char * field, int isample, const char * result)
    {
        if (isample < 0 || isample >= bcf_hdr_nsamples(hdr)) return result;
        const int id = bcf_hdr_id2int(hdr, BCF_DT_ID, field);
        if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, id)) return result;
        const int type = bcf_hdr_id2type(hdr, BCF_HL_FMT, id);
        if (type == BCF_HT_INT)
        {
            const auto values = getFormatInts(hdr, line, field, isample);
            return values.empty() ? result : std::to_string(values.front());
        }
        if (type == BCF_HT_REAL)
        {
            const auto values = getFormatFloats(hdr, line, field, isample);
            return values.empty() ? result : std::to_string(values.front());
        }
        char **raw = nullptr;
        int capacity = 0;
        const int count = bcf_get_format_string(hdr, line, field, &raw, &capacity);
        std::string value = result;
        if (count > 0) value = raw[isample];
        if (raw) { free(raw[0]); free(raw); }
        return value;
    }

    /** update format string for a single sample.  */
    void setFormatStrings(const bcf_hdr_t * hdr, bcf1_t * line, const char * field,
                         const std::vector<std::string> & formats)
    {
        // TODO this can probably be done faster / better
        std::unique_ptr<const char *[]> p_fmts = std::unique_ptr<const char *[]>(new const char *[line->n_sample]);
        bcf_unpack(line, BCF_UN_FMT);
        bool any_nonempty = false;
        if(formats.size() != line->n_sample)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field
               << " " << bcf_hdr_id2name(hdr, line->rid) << ":" << line->pos
               << " -- we have " << formats.size() << " for " << line->n_sample << " samples";
            throw importexception(os.str());
        }
        for (int si = 0; si < line->n_sample; ++si)
        {
            p_fmts.get()[si] = formats[si].c_str();
            if(!formats[si].empty())
            {
                any_nonempty = true;
            }
        }
        int res;
        if(any_nonempty)
        {
            res = bcf_update_format_string(hdr, line, field, p_fmts.get(), line->n_sample);
        }
        else
        {
            res = bcf_update_format_string(hdr, line, field, NULL, 0);
        }

        if(res != 0)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field << " " << bcf_hdr_id2name(hdr, line->rid) << ":" << line->pos
               << " -- we have " << formats.size() << " for " << line->n_sample << " samples";
            throw importexception(os.str());
        }
    }

    /** update format with single float values.  */
    void setFormatFloats(const bcf_hdr_t * header, bcf1_t * line, const char * field,
                         const std::vector<float> & value)
    {
        if(value.empty())
        {
            int res = bcf_update_format(header, line, field, NULL, 0, 0);
            if(res != 0)
            {
                std::ostringstream os;
                os << "[W] cannot update format " << field << " " << bcf_hdr_id2name(header, line->rid) << ":" << line->pos
                   << " -- we have " << value.size() << " for " << line->n_sample << " samples";

                throw importexception(os.str());
            }
            return;
        }
        std::unique_ptr<float[]> p_dbl = std::unique_ptr<float[]>(new float[line->n_sample]);

        for(size_t i = 0; i < line->n_sample; ++i)
        {
            if(i < value.size())
            {
                p_dbl.get()[i] = value[i];
            }
            else
            {
                p_dbl.get()[i] = std::numeric_limits<float>::quiet_NaN();
            }
        }

        int res = bcf_update_format_float(header, line, field, p_dbl.get(), line->n_sample);
        if(res != 0)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field << " " << bcf_hdr_id2name(header, line->rid) << ":" << line->pos
               << " -- we have " << value.size() << " for " << line->n_sample << " samples";
            throw importexception(os.str());
        }
    }

    void setFormatInts(const bcf_hdr_t * header, bcf1_t * line, const char * field,
                       const std::vector<int> & value)
    {
        if(value.empty())
        {
            int res = bcf_update_format(header, line, field, NULL, 0, 0);
            if(res != 0)
            {
                std::ostringstream os;
                os << "[W] cannot update format " << field << " " << bcf_hdr_id2name(header, line->rid) << ":" << line->pos
                   << " -- we have " << value.size() << " for " << line->n_sample << " samples";

                throw importexception(os.str());
            }
            return;
        }
        std::unique_ptr<int[]> p_dbl = std::unique_ptr<int[]>(new int[line->n_sample]);

        for(size_t i = 0; i < line->n_sample; ++i)
        {
            if(i < value.size())
            {
                p_dbl.get()[i] = value[i];
            }
            else
            {
                p_dbl.get()[i] = bcf_int32_missing;
            }
        }

        int res = bcf_update_format_int32(header, line, field, p_dbl.get(), line->n_sample);
        if(res != 0)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field << " " << bcf_hdr_id2name(header, line->rid) << ":" << line->pos
               << " -- we have " << value.size() << " for " << line->n_sample << " samples";
            throw importexception(os.str());
        }
    }

    /** return sample names from header */
    std::list<std::string> getSampleNames(const bcf_hdr_t * hdr)
    {
        std::list<std::string> l;
        for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
        {
            std::string samplename = hdr->samples[i];
            if(samplename == "*")
            {
                std::cerr << "Skipping sample named '*'" << "\n";
                continue;
            }
            l.push_back(samplename);
        }
        return l;
    }

    /** return number of reference padding bases */
    int isRefPadded(bcf1_t * line)
    {
        bcf_unpack(line, BCF_UN_SHR);
        if(line->n_allele == 1)
        {
            return 0;
        }

        if(!line->d.allele)
        {
            throw bcfhelpers::importexception("Empty REF allele -- this is probably because parsing this record failed.");
        }

        const char * ref = line->d.allele[0];
        const int reflen = (int) strlen(ref);

        int max_match = reflen;
        for(int al = 1; al < line->n_allele; ++al)
        {
            const char * alt = line->d.allele[al];
            // symbolic or missing ALT
            if(strcmp(alt, ".") == 0 || *alt == '<')
            {
                max_match = 0;
                break;
            }
            int rpos = 0;
            for(rpos = 0; rpos < reflen; ++rpos, ++alt)
            {
                const char rb = *(ref + rpos);
                if(*alt == 0 || *alt != rb)
                {
                    break;
                }
            }
            max_match = std::min(rpos, max_match);
        }
        return max_match;
    }

    std::pair<AlleleType, std::string> classifyAlleleString(std::string allele_string)
    {
        boost::to_upper(allele_string);

        if (allele_string == "<DEL>" || allele_string == "*")
        {
            return std::make_pair(AlleleType::SYMBOLIC_DEL, std::string());
        }
        else if (allele_string == "" || allele_string == ".")
        {
            return std::make_pair(AlleleType::MISSING, std::string("."));
        }
        else if(allele_string.find('<') != std::string::npos || allele_string.find('>') != std::string::npos)
        {
            return std::make_pair(AlleleType::SYMBOLIC_OTHER, allele_string);
        }
        else
        {
            // check alleles
            for (size_t qq = 0; qq < allele_string.size(); ++qq)
            {
                /* A	A	Adenine */
                /* C	C	Cytosine */
                /* G	G	Guanine */
                /* T	T	Thymine */
                /* U	U	Uracil */
                /* R	A or G	puRine */
                /* Y	C, T or U	pYrimidines */
                /* K	G, T or U	bases which are Ketones */
                /* M	A or C	bases with aMino groups */
                /* S	C or G	Strong interaction */
                /* W	A, T or U	Weak interaction */
                /* B	not A (i.e. C, G, T or U)	B comes after A */
                /* D	not C (i.e. A, G, T or U)	D comes after C */
                /* H	not G (i.e., A, C, T or U)	H comes after G */
                /* V	neither T nor U (i.e. A, C or G)	V comes after U */
                /* N	A C G T U	Nucleic acid */
                /* X	masked */
                /* -/.	gap of indeterminate length - not supported */
                if (allele_string[qq] != 'A' && \
                   allele_string[qq] != 'C' && \
                   allele_string[qq] != 'G' && \
                   allele_string[qq] != 'T' && \
                   allele_string[qq] != 'U' && \
                   allele_string[qq] != 'R' && \
                   allele_string[qq] != 'Y' && \
                   allele_string[qq] != 'K' && \
                   allele_string[qq] != 'M' && \
                   allele_string[qq] != 'S' && \
                   allele_string[qq] != 'W' && \
                   allele_string[qq] != 'B' && \
                   allele_string[qq] != 'D' && \
                   allele_string[qq] != 'H' && \
                   allele_string[qq] != 'V' && \
                   allele_string[qq] != 'N' && \
                   allele_string[qq] != 'X')
                {
                    return std::make_pair(AlleleType::UNKNOWN, allele_string);
                }
            }
            return std::make_pair(AlleleType::NUC, allele_string);
        }
    }

    std::ostream & operator<<(std::ostream & o, const AlleleType at)
    {
        switch(at)
        {
            case AlleleType::NUC:
                o << "NUC";
                break;
            case AlleleType::MISSING:
                o << "MISSING";
                break;
            case AlleleType::SYMBOLIC_DEL:
                o << "SYMBOLIC_DEL";
                break;
            case AlleleType::SYMBOLIC_OTHER:
                o << "SYMBOLIC_OTHER";
                break;
            case AlleleType::UNKNOWN:
                o << "UNKNOWN";
                break;
        }
        return o;
    }
} // namespace bcfhelpers
