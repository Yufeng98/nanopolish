//---------------------------------------------------------
// Copyright 2015 Ontario Institute for Cancer Research
// Written by Jared Simpson (jared.simpson@oicr.on.ca)
//---------------------------------------------------------
//
// nanopolish_emissions -- Emission distributions and
// related functions for the HMMs
//
#ifndef NANOPOLISH_EMISSIONS_H
#define NANOPOLISH_EMISSIONS_H

#include <math.h>
#include "nanopolish_common.h"
#include "nanopolish_squiggle_read.h"

#include "fixed.h"
#include "FixPointCS/Cpp/Fixed64.h"
#include "FloatX/src/floatx.hpp"
#include "flexfloat.hpp"
typedef flexfloat<6, 21> floatc;

using namespace flx;
// typedef floatx<6, 9> floatc;

using namespace Fixed64;
using namespace numeric;
typedef Fixed<18, 14> fixed;
typedef Fixed<20, 12> fixed_long;

//#define DEBUG_HMM_EMISSION 1

// From SO: http://stackoverflow.com/questions/10847007/using-the-gaussian-probability-density-function-in-c
static const float inv_sqrt_2pi = 0.3989422804014327;
inline float normal_pdf(float x, const GaussianParameters& g)
{
    float a = (x - g.mean) / g.stdv;
    return inv_sqrt_2pi / g.stdv * exp(-0.5f * a * a);
}

inline float normal_pdf(float x, const PoreModelStateParams& s)
{
    float a = (x - s.level_mean) / s.level_stdv;
    return inv_sqrt_2pi / s.level_stdv * exp(-0.5f * a * a);
}

inline float z_score(const SquiggleRead& read,
                     const PoreModel& pore_model,
                     uint32_t kmer_rank,
                     uint32_t event_idx,
                     uint8_t strand)
{
    float level = read.get_drift_scaled_level(event_idx, strand);
    GaussianParameters gp = read.get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    return (level - gp.mean) / gp.stdv;
}

static const float log_inv_sqrt_2pi = log(0.3989422804014327);

inline float log_normal_pdf(float x, const PoreModelStateParams& s)
{
    float a = (x - s.level_mean) / s.level_stdv;
    return log_inv_sqrt_2pi - s.level_log_stdv + (-0.5f * a * a);
}

inline float log_normal_pdf(float x, const GaussianParameters& g)
{
    float a = (x - g.mean) / g.stdv;
    return log_inv_sqrt_2pi - g.log_stdv + (-0.5f * a * a);
}

inline float log_normal_pdf_print(float x, const GaussianParameters& g)
{
    float a = (x - g.mean) / g.stdv;
    fprintf(stderr, "fp a %f a_square %f\n", a, a * a);
    return log_inv_sqrt_2pi - g.log_stdv + (-0.5f * a * a);
}

inline float log_probability_match_r9(const SquiggleRead& read,
                                      const PoreModel& pore_model,
                                      uint32_t kmer_rank,
                                      uint32_t event_idx,
                                      uint8_t strand)
{
    // event level mean, scaled with the drift value
    float level = read.get_drift_scaled_level(event_idx, strand);
    GaussianParameters gp = read.get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    float lp = log_normal_pdf(level, gp);
    return lp;
}

inline float log_probability_match_r9_print(const SquiggleRead& read,
                                      const PoreModel& pore_model,
                                      uint32_t kmer_rank,
                                      uint32_t event_idx,
                                      uint8_t strand)
{
    // event level mean, scaled with the drift value
    float level = read.get_drift_scaled_level(event_idx, strand);
    GaussianParameters gp = read.get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    float lp = log_normal_pdf_print(level, gp);
    return lp;
}

inline fixed f_32_log_normal_pdf(fixed x, const GaussianParameters& g) {
    fixed gp_mean = g.f_32_mean;
    fixed gp_stdv = g.f_32_stdv;
    fixed gp_log_stdv = g.f_32_log_stdv;
    fixed f_log_inv_sqrt_2pi = log(0.3989422804014327);
    fixed a = (x - gp_mean) / gp_stdv;
    return f_log_inv_sqrt_2pi - gp_log_stdv + (-0.5f * a * a);
}

inline fixed f_32_log_normal_pdf_a(fixed x, const GaussianParameters& g) {
    fixed gp_mean = g.f_32_mean;
    fixed gp_stdv = g.f_32_stdv;
    fixed gp_log_stdv = g.f_32_log_stdv;
    fixed f_log_inv_sqrt_2pi = log(0.3989422804014327);
    fixed a = (x - gp_mean) / gp_stdv;
    return a;
}

inline floatc fp_log_normal_pdf(floatc x, const GaussianParameters& g) {
    floatc gp_mean = g.fp_mean;
    floatc gp_stdv = g.fp_stdv;
    floatc gp_log_stdv = g.fp_log_stdv;
    floatc f_log_inv_sqrt_2pi = log(0.3989422804014327);
    floatc a = (x - gp_mean) / gp_stdv;
    return f_log_inv_sqrt_2pi - gp_log_stdv + (-0.5f * a * a);
}

inline FP_LONG f_log_normal_pdf(FP_LONG x, const GaussianParameters& g) {
    FP_LONG gp_mean = g.f_mean;
    FP_LONG gp_stdv = g.f_stdv;
    FP_LONG gp_log_stdv = g.f_log_stdv;
    FP_LONG log_inv_sqrt_2pi = FromFloat(log(0.3989422804014327));
    FP_LONG a = DivPrecise(Sub(x, gp_mean), gp_stdv);
    FP_LONG a_square = Mul(a, a);
    FP_LONG a_square_tmp = Mul(a_square, FromFloat(-0.5));
    FP_LONG result_tmp = Sub(log_inv_sqrt_2pi, gp_log_stdv);
    // fprintf(stderr, "gp_stdv %f a %f\n", ToFloat(gp_stdv), ToFloat(a));
    return Add(result_tmp, a_square_tmp);
}

inline FP_LONG f_log_normal_pdf_print(FP_LONG x, const GaussianParameters& g) {
    FP_LONG gp_mean = g.f_mean;
    FP_LONG gp_stdv = g.f_stdv;
    FP_LONG gp_log_stdv = g.f_log_stdv;
    FP_LONG log_inv_sqrt_2pi = FromFloat(log(0.3989422804014327));
    FP_LONG a = DivPrecise(Sub(x, gp_mean), gp_stdv);
    FP_LONG a_square = Mul(a, a);
    FP_LONG a_square_tmp = Mul(a_square, FromFloat(-0.5));
    FP_LONG result_tmp = Sub(log_inv_sqrt_2pi, gp_log_stdv);
    fprintf(stderr, "fix a %f a_square %f\n", ToFloat(a), ToFloat(a_square));
    return Add(result_tmp, a_square_tmp);
}

inline fixed f_32_log_probability_match_r9(const SquiggleRead& read,
                                      const PoreModel& pore_model,
                                      uint32_t kmer_rank,
                                      uint32_t event_idx,
                                      uint8_t strand)
{
    // event level mean, scaled with the drift value
    fixed level = read.f_32_get_drift_scaled_level(event_idx, strand);
    GaussianParameters gp = read.f_32_get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    fixed lp = f_32_log_normal_pdf(level, gp);
    return lp;
}

inline fixed f_32_log_probability_match_r9_a(const SquiggleRead& read,
                                      const PoreModel& pore_model,
                                      uint32_t kmer_rank,
                                      uint32_t event_idx,
                                      uint8_t strand)
{
    // event level mean, scaled with the drift value
    fixed level = read.f_32_get_drift_scaled_level(event_idx, strand);
    GaussianParameters gp = read.f_32_get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    fixed a = f_32_log_normal_pdf_a(level, gp);
    return a;
}

inline floatc fp_log_probability_match_r9(const SquiggleRead& read,
                                      const PoreModel& pore_model,
                                      uint32_t kmer_rank,
                                      uint32_t event_idx,
                                      uint8_t strand)
{
    // event level mean, scaled with the drift value
    floatc level = read.fp_get_drift_scaled_level(event_idx, strand);
    GaussianParameters gp = read.fp_get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    floatc lp = fp_log_normal_pdf(level, gp);
    return lp;
}

inline FP_LONG f_log_probability_match_r9(const SquiggleRead& read,
                                      const PoreModel& pore_model,
                                      uint32_t kmer_rank,
                                      uint32_t event_idx,
                                      uint8_t strand)
{
    // event level mean, scaled with the drift value
    FP_LONG level = read.f_get_drift_scaled_level(event_idx, strand);
    // fprintf(stderr, "level %f\n", ToFloat(level));
    GaussianParameters gp = read.f_get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    FP_LONG lp = f_log_normal_pdf(level, gp);
    return lp;
}

inline FP_LONG f_log_probability_match_r9_print(const SquiggleRead& read,
                                      const PoreModel& pore_model,
                                      uint32_t kmer_rank,
                                      uint32_t event_idx,
                                      uint8_t strand)
{
    // event level mean, scaled with the drift value
    FP_LONG level = read.f_get_drift_scaled_level(event_idx, strand);
    GaussianParameters gp = read.f_get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    FP_LONG lp = f_log_normal_pdf_print(level, gp);
    return lp;
}

inline float log_probability_match_r7(const SquiggleRead& read,
                                      const PoreModel& pore_model,
                                      uint32_t kmer_rank,
                                      uint32_t event_idx,
                                      uint8_t strand,
                                      float state_scale = 1.0f,
                                      float log_state_scale = 0.0f)
{
    float level = read.get_drift_scaled_level(event_idx, strand);
    GaussianParameters gp = read.get_scaled_gaussian_from_pore_model_state(pore_model, strand, kmer_rank);
    gp.stdv *= state_scale;
    gp.log_stdv += log_state_scale;
    float lp = log_normal_pdf(level, gp);
    return lp;
}

inline float log_probability_event_insert_r7(const SquiggleRead& read,
                                             const PoreModel& pore_model,
                                             uint32_t kmer_rank,
                                             uint32_t event_idx,
                                             uint8_t strand)
{
    static const float scale = 1.75f;
    static const float log_scale = log(scale);

    return log_probability_match_r7(read, pore_model, kmer_rank, event_idx, strand, scale, log_scale);
}

inline float log_probability_background(const SquiggleRead&,
                                        uint32_t,
                                        uint8_t)
{
    return -3.0f;
}

#endif
