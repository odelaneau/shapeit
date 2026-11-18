/*******************************************************************************
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef _HAPLOTYPE_SEGMENT_DOUBLE_H
#define _HAPLOTYPE_SEGMENT_DOUBLE_H

#include <objects/compute_job.h>
#include <objects/hmm_parameters.h>
#include <utils/otools.h>

#ifdef __AVX2__
#include <immintrin.h>
#elif defined(__ARM_NEON)
#include <arm_neon.h>
#endif
#include <boost/align/aligned_allocator.hpp>

template <typename T>
using aligned_vector32 =
    std::vector<T, boost::alignment::aligned_allocator<T, 32>>;

class haplotype_segment_double {
private:
  // EXTERNAL DATA
  hmm_parameters &M;
  genotype *G;
  bitmatrix Hhap, Hvar;

  // COORDINATES & CONSTANTS
  int segment_first;
  int segment_last;
  int locus_first;
  int locus_last;
  int ambiguous_first;
  int ambiguous_last;
  int missing_first;
  int missing_last;
  int transition_first;
  int transition_last;
  unsigned int n_cond_haps;
  unsigned int n_missing;

  // CURSORS
  int curr_segment_index;
  int curr_segment_locus;
  int curr_abs_locus;
  int prev_abs_locus;
  int curr_rel_locus;
  int curr_rel_locus_offset;
  int curr_abs_ambiguous;
  int curr_abs_transition;
  int curr_abs_missing;
  int curr_rel_missing;

  // DYNAMIC ARRAYS
  double probSumT;
  aligned_vector32<double> prob;
  aligned_vector32<double> probSumK;
  aligned_vector32<double> probSumH;
  std::vector<aligned_vector32<double>> Alpha;
  std::vector<aligned_vector32<double>> AlphaSum;
  std::vector<int> AlphaLocus;
  aligned_vector32<double> AlphaSumSum;
  std::vector<aligned_vector32<double>> AlphaMissing;
  std::vector<aligned_vector32<double>> AlphaSumMissing;
  double HProbs[HAP_NUMBER * HAP_NUMBER] __attribute__((aligned(32)));
  double DProbs[HAP_NUMBER * HAP_NUMBER * HAP_NUMBER * HAP_NUMBER]
      __attribute__((aligned(32)));

  // STATIC ARRAYS
  double sumHProbs;
  double sumDProbs;
  double g0[HAP_NUMBER], g1[HAP_NUMBER];
  double nt, yt;

  // INLINED AND UNROLLED ROUTINES
  void INIT_HOM();
  void INIT_AMB();
  void INIT_MIS();
  bool RUN_HOM(char);
  void RUN_AMB();
  void RUN_MIS();
  void COLLAPSE_HOM();
  void COLLAPSE_AMB();
  void COLLAPSE_MIS();
  void SUMK();
  void IMPUTE(std::vector<float> &);
  bool TRANS_HAP();
  bool TRANS_DIP_MULT();
  bool TRANS_DIP_ADD();
  void SET_FIRST_TRANS(std::vector<double> &);
  int SET_OTHER_TRANS(std::vector<double> &);

public:
  // CONSTRUCTOR/DESTRUCTOR
  haplotype_segment_double(genotype *, bitmatrix &, std::vector<unsigned int> &,
                           window &, hmm_parameters &);
  ~haplotype_segment_double();

  // void fetch();
  void forward();
  int backward(std::vector<double> &, std::vector<float> &);
};

/*******************************************************************************/
/*****************			HOMOZYGOUS GENOTYPE
 * ************************/
/*******************************************************************************/

inline void haplotype_segment_double::INIT_HOM() {
  bool ag =
      VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
#if defined(__AVX2__)
  __m256d _sum0 = _mm256_set1_pd(0.0f);
  __m256d _sum1 = _mm256_set1_pd(0.0f);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    __m256d _prob0 = _mm256_set1_pd((ag == ah) ? 1.0f : M.ed / M.ee);
    __m256d _prob1 = _mm256_set1_pd((ag == ah) ? 1.0f : M.ed / M.ee);
    _sum0 = _mm256_add_pd(_sum0, _prob0);
    _sum1 = _mm256_add_pd(_sum1, _prob1);
    _mm256_store_pd(&prob[i + 0], _prob0);
    _mm256_store_pd(&prob[i + 4], _prob1);
  }
  _mm256_store_pd(&probSumH[0], _sum0);
  _mm256_store_pd(&probSumH[4], _sum1);
#elif defined(__ARM_NEON)
  float64x2_t _sum = vdupq_n_f64(0.0);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float64x2_t _prob = vdupq_n_f64((ag == ah) ? 1.0 : M.ed / M.ee);
    _sum = vaddq_f64(_sum, _prob);
    vst1q_f64(&prob[i + 0], _prob);
    vst1q_f64(&prob[i + 2], _prob);
    vst1q_f64(&prob[i + 4], _prob);
    vst1q_f64(&prob[i + 6], _prob);
  }
  vst1q_f64(&probSumH[0], _sum);
  vst1q_f64(&probSumH[2], _sum);
  vst1q_f64(&probSumH[4], _sum);
  vst1q_f64(&probSumH[6], _sum);
#else
  double sum0 = 0.0, sum1 = 0.0;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    double p = (ag == ah) ? 1.0 : M.ed / M.ee;
    for (int h = 0; h < HAP_NUMBER; h++)
      prob[i + h] = p;
    sum0 += 4.0 * p;
    sum1 += 4.0 * p;
  }
  for (int h = 0; h < 4; h++)
    probSumH[h] = sum0;
  for (int h = 4; h < HAP_NUMBER; h++)
    probSumH[h] = sum1;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline bool haplotype_segment_double::RUN_HOM(char rare_allele) {
  bool ag =
      VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
  if (rare_allele < 0 || ag == rare_allele) {
#if defined(__AVX2__)
    __m256d _sum0 = _mm256_set1_pd(0.0f);
    __m256d _sum1 = _mm256_set1_pd(0.0f);
    __m256d _factor = _mm256_set1_pd(yt / (n_cond_haps * probSumT));
    __m256d _tFreq0 = _mm256_load_pd(&probSumH[0]);
    __m256d _tFreq1 = _mm256_load_pd(&probSumH[4]);
    _tFreq0 = _mm256_mul_pd(_tFreq0, _factor);
    _tFreq1 = _mm256_mul_pd(_tFreq1, _factor);
    __m256d _nt = _mm256_set1_pd(nt / probSumT);
    __m256d _mismatch = _mm256_set1_pd(M.ed / M.ee);
    for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
      bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
      __m256d _prob0 = _mm256_load_pd(&prob[i]);
      __m256d _prob1 = _mm256_load_pd(&prob[i + 4]);
      _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
      _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
      if (ag != ah) {
        _prob0 = _mm256_mul_pd(_prob0, _mismatch);
        _prob1 = _mm256_mul_pd(_prob1, _mismatch);
      }
      _sum0 = _mm256_add_pd(_sum0, _prob0);
      _sum1 = _mm256_add_pd(_sum1, _prob1);
      _mm256_store_pd(&prob[i], _prob0);
      _mm256_store_pd(&prob[i + 4], _prob1);
    }
    _mm256_store_pd(&probSumH[0], _sum0);
    _mm256_store_pd(&probSumH[4], _sum1);
#elif defined(__ARM_NEON)
    float64x2_t _sum0 = vdupq_n_f64(0.0);
    float64x2_t _sum1 = vdupq_n_f64(0.0);
    float64x2_t _sum2 = vdupq_n_f64(0.0);
    float64x2_t _sum3 = vdupq_n_f64(0.0);
    float64x2_t _factor = vdupq_n_f64(yt / (n_cond_haps * probSumT));
    float64x2_t _tf0 = vld1q_f64(&probSumH[0]);
    float64x2_t _tf1 = vld1q_f64(&probSumH[2]);
    float64x2_t _tf2 = vld1q_f64(&probSumH[4]);
    float64x2_t _tf3 = vld1q_f64(&probSumH[6]);
    _tf0 = vmulq_f64(_tf0, _factor);
    _tf1 = vmulq_f64(_tf1, _factor);
    _tf2 = vmulq_f64(_tf2, _factor);
    _tf3 = vmulq_f64(_tf3, _factor);
    float64x2_t _nt = vdupq_n_f64(nt / probSumT);
    float64x2_t _mismatch = vdupq_n_f64(M.ed / M.ee);
    for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
      bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
      float64x2_t _p0 = vld1q_f64(&prob[i]);
      float64x2_t _p1 = vld1q_f64(&prob[i + 2]);
      float64x2_t _p2 = vld1q_f64(&prob[i + 4]);
      float64x2_t _p3 = vld1q_f64(&prob[i + 6]);
      _p0 = vmlaq_f64(_tf0, _p0, _nt);
      _p1 = vmlaq_f64(_tf1, _p1, _nt);
      _p2 = vmlaq_f64(_tf2, _p2, _nt);
      _p3 = vmlaq_f64(_tf3, _p3, _nt);
      if (ag != ah) {
        _p0 = vmulq_f64(_p0, _mismatch);
        _p1 = vmulq_f64(_p1, _mismatch);
        _p2 = vmulq_f64(_p2, _mismatch);
        _p3 = vmulq_f64(_p3, _mismatch);
      }
      _sum0 = vaddq_f64(_sum0, _p0);
      _sum1 = vaddq_f64(_sum1, _p1);
      _sum2 = vaddq_f64(_sum2, _p2);
      _sum3 = vaddq_f64(_sum3, _p3);
      vst1q_f64(&prob[i], _p0);
      vst1q_f64(&prob[i + 2], _p1);
      vst1q_f64(&prob[i + 4], _p2);
      vst1q_f64(&prob[i + 6], _p3);
    }
    vst1q_f64(&probSumH[0], _sum0);
    vst1q_f64(&probSumH[2], _sum1);
    vst1q_f64(&probSumH[4], _sum2);
    vst1q_f64(&probSumH[6], _sum3);
#else
    double sum0 = 0.0, sum1 = 0.0;
    double factor = yt / (n_cond_haps * probSumT);
    double tFreq0[4], tFreq1[4];
    for (int h = 0; h < 4; h++) {
      tFreq0[h] = probSumH[h] * factor;
      tFreq1[h] = probSumH[h + 4] * factor;
    }
    double nt_ = nt / probSumT;
    double mismatch = M.ed / M.ee;
    for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
      bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
      for (int h = 0; h < 4; h++) {
        prob[i + h] = prob[i + h] * nt_ + tFreq0[h];
        prob[i + 4 + h] = prob[i + 4 + h] * nt_ + tFreq1[h];
        if (ag != ah) {
          prob[i + h] *= mismatch;
          prob[i + 4 + h] *= mismatch;
        }
      }
      for (int h = 0; h < 4; h++) {
        sum0 += prob[i + h];
        sum1 += prob[i + 4 + h];
      }
    }
    for (int h = 0; h < 4; h++)
      probSumH[h] = sum0;
    for (int h = 4; h < HAP_NUMBER; h++)
      probSumH[h] = sum1;
#endif
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
               probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
    return true;
  }
  return false;
}

inline void haplotype_segment_double::COLLAPSE_HOM() {
  bool ag =
      VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
#if defined(__AVX2__)
  __m256d _sum0 = _mm256_set1_pd(0.0f);
  __m256d _sum1 = _mm256_set1_pd(0.0f);
  __m256d _tFreq =
      _mm256_set1_pd(yt / n_cond_haps); // Check divide by probSumT here!
  __m256d _nt = _mm256_set1_pd(nt / probSumT);
  __m256d _mismatch = _mm256_set1_pd(M.ed / M.ee);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    __m256d _prob0 = _mm256_set1_pd(probSumK[k]);
    __m256d _prob1 = _mm256_set1_pd(probSumK[k]);
    _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq);
    _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq);
    if (ag != ah) {
      _prob0 = _mm256_mul_pd(_prob0, _mismatch);
      _prob1 = _mm256_mul_pd(_prob1, _mismatch);
    }
    _sum0 = _mm256_add_pd(_sum0, _prob0);
    _sum1 = _mm256_add_pd(_sum1, _prob1);
    _mm256_store_pd(&prob[i], _prob0);
    _mm256_store_pd(&prob[i + 4], _prob1);
  }
  _mm256_store_pd(&probSumH[0], _sum0);
  _mm256_store_pd(&probSumH[4], _sum1);
#elif defined(__ARM_NEON)
  float64x2_t _sum = vdupq_n_f64(0.0);
  float64x2_t _tFreq =
      vdupq_n_f64(yt / n_cond_haps); // Check divide by probSumT here!
  float64x2_t _nt = vdupq_n_f64(nt / probSumT);
  float64x2_t _mismatch = vdupq_n_f64(M.ed / M.ee);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float64x2_t _prob = vmlaq_f64(_tFreq, vdupq_n_f64(probSumK[k]), _nt);
    if (ag != ah)
      _prob = vmulq_f64(_prob, _mismatch);
    _sum = vaddq_f64(_sum, _prob);
    vst1q_f64(&prob[i], _prob);
    vst1q_f64(&prob[i + 2], _prob);
    vst1q_f64(&prob[i + 4], _prob);
    vst1q_f64(&prob[i + 6], _prob);
  }
  vst1q_f64(&probSumH[0], _sum);
  vst1q_f64(&probSumH[2], _sum);
  vst1q_f64(&probSumH[4], _sum);
  vst1q_f64(&probSumH[6], _sum);
#else
  double sum0 = 0.0, sum1 = 0.0;
  double tFreq = yt / n_cond_haps; // Check divide by probSumT here!
  double nt_ = nt / probSumT;
  double mismatch = M.ed / M.ee;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    double p0 = probSumK[k] * nt_ + tFreq;
    double p1 = probSumK[k] * nt_ + tFreq;
    if (ag != ah) {
      p0 *= mismatch;
      p1 *= mismatch;
    }
    for (int h = 0; h < 4; h++)
      prob[i + h] = p0;
    for (int h = 4; h < HAP_NUMBER; h++)
      prob[i + h] = p1;
    sum0 += 4.0 * p0;
    sum1 += 4.0 * p1;
  }
  for (int h = 0; h < 4; h++)
    probSumH[h] = sum0;
  for (int h = 4; h < HAP_NUMBER; h++)
    probSumH[h] = sum1;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************			HETEROZYGOUS GENOTYPE
 * ********************/
/*******************************************************************************/

inline void haplotype_segment_double::INIT_AMB() {
  unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
  for (int h = 0; h < HAP_NUMBER; h++) {
    g0[h] = HAP_GET(amb_code, h) ? M.ed / M.ee : 1.0f;
    g1[h] = HAP_GET(amb_code, h) ? 1.0f : M.ed / M.ee;
  }
#ifdef __AVX2__
  __m256d _sum0 = _mm256_set1_pd(0.0f);
  __m256d _sum1 = _mm256_set1_pd(0.0f);
  __m256d _emit0[2], _emit1[2];
  _emit0[0] = _mm256_loadu_pd(&g0[0]);
  _emit0[1] = _mm256_loadu_pd(&g1[0]);
  _emit1[0] = _mm256_loadu_pd(&g0[4]);
  _emit1[1] = _mm256_loadu_pd(&g1[4]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    __m256d _prob0 = _emit0[ah];
    __m256d _prob1 = _emit1[ah];
    _sum0 = _mm256_add_pd(_sum0, _prob0);
    _sum1 = _mm256_add_pd(_sum1, _prob1);
    _mm256_store_pd(&prob[i], _prob0);
    _mm256_store_pd(&prob[i + 4], _prob1);
  }
  _mm256_store_pd(&probSumH[0], _sum0);
  _mm256_store_pd(&probSumH[4], _sum1);
#elif defined(__ARM_NEON)
  float64x2_t _sum0 = vdupq_n_f64(0.0);
  float64x2_t _sum1 = vdupq_n_f64(0.0);
  float64x2_t _sum2 = vdupq_n_f64(0.0);
  float64x2_t _sum3 = vdupq_n_f64(0.0);
  float64x2_t _e00 = vld1q_f64(&g0[0]), _e01 = vld1q_f64(&g0[2]);
  float64x2_t _e02 = vld1q_f64(&g1[0]), _e03 = vld1q_f64(&g1[2]);
  float64x2_t _e10 = vld1q_f64(&g0[4]), _e11 = vld1q_f64(&g0[6]);
  float64x2_t _e12 = vld1q_f64(&g1[4]), _e13 = vld1q_f64(&g1[6]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float64x2_t _p0 = ah ? _e02 : _e00;
    float64x2_t _p1 = ah ? _e03 : _e01;
    float64x2_t _p2 = ah ? _e12 : _e10;
    float64x2_t _p3 = ah ? _e13 : _e11;
    _sum0 = vaddq_f64(_sum0, _p0);
    _sum1 = vaddq_f64(_sum1, _p1);
    _sum2 = vaddq_f64(_sum2, _p2);
    _sum3 = vaddq_f64(_sum3, _p3);
    vst1q_f64(&prob[i], _p0);
    vst1q_f64(&prob[i + 2], _p1);
    vst1q_f64(&prob[i + 4], _p2);
    vst1q_f64(&prob[i + 6], _p3);
  }
  vst1q_f64(&probSumH[0], _sum0);
  vst1q_f64(&probSumH[2], _sum1);
  vst1q_f64(&probSumH[4], _sum2);
  vst1q_f64(&probSumH[6], _sum3);
#else
  double sum0 = 0.0, sum1 = 0.0;
  for (int h = 0; h < 4; h++) {
    sum0 += g0[h];
    sum1 += g1[h];
  }
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    for (int h = 0; h < 4; h++)
      prob[i + h] = (ah ? g1[h] : g0[h]);
    for (int h = 4; h < HAP_NUMBER; h++)
      prob[i + h] = (ah ? g1[h - 4] : g0[h - 4]);
  }
  for (int h = 0; h < 4; h++)
    probSumH[h] = sum0;
  for (int h = 4; h < HAP_NUMBER; h++)
    probSumH[h] = sum1;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline void haplotype_segment_double::RUN_AMB() {
  unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
  for (int h = 0; h < HAP_NUMBER; h++) {
    g0[h] = HAP_GET(amb_code, h) ? M.ed / M.ee : 1.0f;
    g1[h] = HAP_GET(amb_code, h) ? 1.0f : M.ed / M.ee;
  }
#ifdef __AVX2__
  __m256d _sum0 = _mm256_set1_pd(0.0f);
  __m256d _sum1 = _mm256_set1_pd(0.0f);
  __m256d _factor = _mm256_set1_pd(yt / (n_cond_haps * probSumT));
  __m256d _tFreq0 = _mm256_load_pd(&probSumH[0]);
  __m256d _tFreq1 = _mm256_load_pd(&probSumH[4]);
  _tFreq0 = _mm256_mul_pd(_tFreq0, _factor);
  _tFreq1 = _mm256_mul_pd(_tFreq1, _factor);
  __m256d _nt = _mm256_set1_pd(nt / probSumT);
  __m256d _emit0[2], _emit1[2];
  _emit0[0] = _mm256_loadu_pd(&g0[0]);
  _emit0[1] = _mm256_loadu_pd(&g1[0]);
  _emit1[0] = _mm256_loadu_pd(&g0[4]);
  _emit1[1] = _mm256_loadu_pd(&g1[4]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    __m256d _prob0 = _mm256_load_pd(&prob[i + 0]);
    __m256d _prob1 = _mm256_load_pd(&prob[i + 4]);
    _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
    _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
    _prob0 = _mm256_mul_pd(_prob0, _emit0[ah]);
    _prob1 = _mm256_mul_pd(_prob1, _emit1[ah]);
    _sum0 = _mm256_add_pd(_sum0, _prob0);
    _sum1 = _mm256_add_pd(_sum1, _prob1);
    _mm256_store_pd(&prob[i + 0], _prob0);
    _mm256_store_pd(&prob[i + 4], _prob1);
  }
  _mm256_store_pd(&probSumH[0], _sum0);
  _mm256_store_pd(&probSumH[4], _sum1);
#elif defined(__ARM_NEON)
  float64x2_t _sum0 = vdupq_n_f64(0.0);
  float64x2_t _sum1 = vdupq_n_f64(0.0);
  float64x2_t _sum2 = vdupq_n_f64(0.0);
  float64x2_t _sum3 = vdupq_n_f64(0.0);
  float64x2_t _factor = vdupq_n_f64(yt / (n_cond_haps * probSumT));
  float64x2_t _tf0 = vld1q_f64(&probSumH[0]);
  float64x2_t _tf1 = vld1q_f64(&probSumH[2]);
  float64x2_t _tf2 = vld1q_f64(&probSumH[4]);
  float64x2_t _tf3 = vld1q_f64(&probSumH[6]);
  _tf0 = vmulq_f64(_tf0, _factor);
  _tf1 = vmulq_f64(_tf1, _factor);
  _tf2 = vmulq_f64(_tf2, _factor);
  _tf3 = vmulq_f64(_tf3, _factor);
  float64x2_t _nt = vdupq_n_f64(nt / probSumT);
  float64x2_t _e00 = vld1q_f64(&g0[0]), _e01 = vld1q_f64(&g0[2]);
  float64x2_t _e02 = vld1q_f64(&g1[0]), _e03 = vld1q_f64(&g1[2]);
  float64x2_t _e10 = vld1q_f64(&g0[4]), _e11 = vld1q_f64(&g0[6]);
  float64x2_t _e12 = vld1q_f64(&g1[4]), _e13 = vld1q_f64(&g1[6]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float64x2_t _p0 = vld1q_f64(&prob[i]);
    float64x2_t _p1 = vld1q_f64(&prob[i + 2]);
    float64x2_t _p2 = vld1q_f64(&prob[i + 4]);
    float64x2_t _p3 = vld1q_f64(&prob[i + 6]);
    _p0 = vmulq_f64(vmlaq_f64(_tf0, _p0, _nt), ah ? _e02 : _e00);
    _p1 = vmulq_f64(vmlaq_f64(_tf1, _p1, _nt), ah ? _e03 : _e01);
    _p2 = vmulq_f64(vmlaq_f64(_tf2, _p2, _nt), ah ? _e12 : _e10);
    _p3 = vmulq_f64(vmlaq_f64(_tf3, _p3, _nt), ah ? _e13 : _e11);
    _sum0 = vaddq_f64(_sum0, _p0);
    _sum1 = vaddq_f64(_sum1, _p1);
    _sum2 = vaddq_f64(_sum2, _p2);
    _sum3 = vaddq_f64(_sum3, _p3);
    vst1q_f64(&prob[i], _p0);
    vst1q_f64(&prob[i + 2], _p1);
    vst1q_f64(&prob[i + 4], _p2);
    vst1q_f64(&prob[i + 6], _p3);
  }
  vst1q_f64(&probSumH[0], _sum0);
  vst1q_f64(&probSumH[2], _sum1);
  vst1q_f64(&probSumH[4], _sum2);
  vst1q_f64(&probSumH[6], _sum3);
#else
  double sum0 = 0.0, sum1 = 0.0;
  double factor = yt / (n_cond_haps * probSumT);
  double tFreq0[4], tFreq1[4];
  for (int h = 0; h < 4; h++) {
    tFreq0[h] = probSumH[h] * factor;
    tFreq1[h] = probSumH[h + 4] * factor;
  }
  double nt_ = nt / probSumT;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    for (int h = 0; h < 4; h++) {
      prob[i + h] = (prob[i + h] * nt_ + tFreq0[h]) * (ah ? g1[h] : g0[h]);
      prob[i + 4 + h] =
          (prob[i + 4 + h] * nt_ + tFreq1[h]) * (ah ? g1[h - 4] : g0[h - 4]);
      sum0 += prob[i + h];
      sum1 += prob[i + 4 + h];
    }
  }
  for (int h = 0; h < 4; h++)
    probSumH[h] = sum0;
  for (int h = 4; h < HAP_NUMBER; h++)
    probSumH[h] = sum1;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline void haplotype_segment_double::COLLAPSE_AMB() {
  unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
  for (int h = 0; h < HAP_NUMBER; h++) {
    g0[h] = HAP_GET(amb_code, h) ? M.ed / M.ee : 1.0f;
    g1[h] = HAP_GET(amb_code, h) ? 1.0f : M.ed / M.ee;
  }
#ifdef __AVX2__
  __m256d _sum0 = _mm256_set1_pd(0.0f);
  __m256d _sum1 = _mm256_set1_pd(0.0f);
  __m256d _tFreq = _mm256_set1_pd(yt / n_cond_haps);
  __m256d _nt = _mm256_set1_pd(nt / probSumT);
  __m256d _emit0[2], _emit1[2];
  _emit0[0] = _mm256_loadu_pd(&g0[0]);
  _emit0[1] = _mm256_loadu_pd(&g1[0]);
  _emit1[0] = _mm256_loadu_pd(&g0[4]);
  _emit1[1] = _mm256_loadu_pd(&g1[4]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    __m256d _prob0 = _mm256_set1_pd(probSumK[k]);
    __m256d _prob1 = _mm256_set1_pd(probSumK[k]);
    _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq);
    _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq);
    _prob0 = _mm256_mul_pd(_prob0, _emit0[ah]);
    _prob1 = _mm256_mul_pd(_prob1, _emit1[ah]);
    _sum0 = _mm256_add_pd(_sum0, _prob0);
    _sum1 = _mm256_add_pd(_sum1, _prob1);
    _mm256_store_pd(&prob[i + 0], _prob0);
    _mm256_store_pd(&prob[i + 4], _prob1);
  }
  _mm256_store_pd(&probSumH[0], _sum0);
  _mm256_store_pd(&probSumH[4], _sum1);
#elif defined(__ARM_NEON)
  float64x2_t _sum0 = vdupq_n_f64(0.0);
  float64x2_t _sum1 = vdupq_n_f64(0.0);
  float64x2_t _sum2 = vdupq_n_f64(0.0);
  float64x2_t _sum3 = vdupq_n_f64(0.0);
  float64x2_t _tFreq = vdupq_n_f64(yt / n_cond_haps);
  float64x2_t _nt = vdupq_n_f64(nt / probSumT);
  float64x2_t _e00 = vld1q_f64(&g0[0]), _e01 = vld1q_f64(&g0[2]);
  float64x2_t _e02 = vld1q_f64(&g1[0]), _e03 = vld1q_f64(&g1[2]);
  float64x2_t _e10 = vld1q_f64(&g0[4]), _e11 = vld1q_f64(&g0[6]);
  float64x2_t _e12 = vld1q_f64(&g1[4]), _e13 = vld1q_f64(&g1[6]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float64x2_t _base = vmlaq_f64(_tFreq, vdupq_n_f64(probSumK[k]), _nt);
    float64x2_t _p0 = vmulq_f64(_base, ah ? _e02 : _e00);
    float64x2_t _p1 = vmulq_f64(_base, ah ? _e03 : _e01);
    float64x2_t _p2 = vmulq_f64(_base, ah ? _e12 : _e10);
    float64x2_t _p3 = vmulq_f64(_base, ah ? _e13 : _e11);
    _sum0 = vaddq_f64(_sum0, _p0);
    _sum1 = vaddq_f64(_sum1, _p1);
    _sum2 = vaddq_f64(_sum2, _p2);
    _sum3 = vaddq_f64(_sum3, _p3);
    vst1q_f64(&prob[i], _p0);
    vst1q_f64(&prob[i + 2], _p1);
    vst1q_f64(&prob[i + 4], _p2);
    vst1q_f64(&prob[i + 6], _p3);
  }
  vst1q_f64(&probSumH[0], _sum0);
  vst1q_f64(&probSumH[2], _sum1);
  vst1q_f64(&probSumH[4], _sum2);
  vst1q_f64(&probSumH[6], _sum3);
#else
  double sum0 = 0.0, sum1 = 0.0;
  double tFreq = yt / n_cond_haps;
  double nt_ = nt / probSumT;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    double base = probSumK[k] * nt_ + tFreq;
    for (int h = 0; h < 4; h++) {
      prob[i + h] = base * (ah ? g1[h] : g0[h]);
      sum0 += prob[i + h];
    }
    for (int h = 4; h < HAP_NUMBER; h++) {
      prob[i + h] = base * (ah ? g1[h - 4] : g0[h - 4]);
      sum1 += prob[i + h];
    }
  }
  for (int h = 0; h < 4; h++)
    probSumH[h] = sum0;
  for (int h = 4; h < HAP_NUMBER; h++)
    probSumH[h] = sum1;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************			MISSING GENOTYPE
 * ************************/
/*******************************************************************************/

inline void haplotype_segment_double::INIT_MIS() {
  fill(prob.begin(), prob.end(), 1.0f / (HAP_NUMBER * n_cond_haps));
  fill(probSumH.begin(), probSumH.end(), 1.0f / HAP_NUMBER);
  probSumT = 1.0f;
}

inline void haplotype_segment_double::RUN_MIS() {
#if defined(__AVX2__)
  __m256d _sum0 = _mm256_set1_pd(0.0f);
  __m256d _sum1 = _mm256_set1_pd(0.0f);
  __m256d _factor = _mm256_set1_pd(yt / (n_cond_haps * probSumT));
  __m256d _tFreq0 = _mm256_load_pd(&probSumH[0]);
  __m256d _tFreq1 = _mm256_load_pd(&probSumH[4]);
  _tFreq0 = _mm256_mul_pd(_tFreq0, _factor);
  _tFreq1 = _mm256_mul_pd(_tFreq1, _factor);
  __m256d _nt = _mm256_set1_pd(nt / probSumT);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    __m256d _prob0 = _mm256_load_pd(&prob[i]);
    __m256d _prob1 = _mm256_load_pd(&prob[i + 4]);
    _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
    _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
    _sum0 = _mm256_add_pd(_sum0, _prob0);
    _sum1 = _mm256_add_pd(_sum1, _prob1);
    _mm256_store_pd(&prob[i], _prob0);
    _mm256_store_pd(&prob[i + 4], _prob1);
  }
  _mm256_store_pd(&probSumH[0], _sum0);
  _mm256_store_pd(&probSumH[4], _sum1);
#elif defined(__ARM_NEON)
  float64x2_t _sum0 = vdupq_n_f64(0.0);
  float64x2_t _sum1 = vdupq_n_f64(0.0);
  float64x2_t _sum2 = vdupq_n_f64(0.0);
  float64x2_t _sum3 = vdupq_n_f64(0.0);
  float64x2_t _factor = vdupq_n_f64(yt / (n_cond_haps * probSumT));
  float64x2_t _tf0 = vld1q_f64(&probSumH[0]);
  float64x2_t _tf1 = vld1q_f64(&probSumH[2]);
  float64x2_t _tf2 = vld1q_f64(&probSumH[4]);
  float64x2_t _tf3 = vld1q_f64(&probSumH[6]);
  _tf0 = vmulq_f64(_tf0, _factor);
  _tf1 = vmulq_f64(_tf1, _factor);
  _tf2 = vmulq_f64(_tf2, _factor);
  _tf3 = vmulq_f64(_tf3, _factor);
  float64x2_t _nt = vdupq_n_f64(nt / probSumT);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    float64x2_t _p0 = vld1q_f64(&prob[i]);
    float64x2_t _p1 = vld1q_f64(&prob[i + 2]);
    float64x2_t _p2 = vld1q_f64(&prob[i + 4]);
    float64x2_t _p3 = vld1q_f64(&prob[i + 6]);
    _p0 = vmlaq_f64(_tf0, _p0, _nt);
    _p1 = vmlaq_f64(_tf1, _p1, _nt);
    _p2 = vmlaq_f64(_tf2, _p2, _nt);
    _p3 = vmlaq_f64(_tf3, _p3, _nt);
    _sum0 = vaddq_f64(_sum0, _p0);
    _sum1 = vaddq_f64(_sum1, _p1);
    _sum2 = vaddq_f64(_sum2, _p2);
    _sum3 = vaddq_f64(_sum3, _p3);
    vst1q_f64(&prob[i], _p0);
    vst1q_f64(&prob[i + 2], _p1);
    vst1q_f64(&prob[i + 4], _p2);
    vst1q_f64(&prob[i + 6], _p3);
  }
  vst1q_f64(&probSumH[0], _sum0);
  vst1q_f64(&probSumH[2], _sum1);
  vst1q_f64(&probSumH[4], _sum2);
  vst1q_f64(&probSumH[6], _sum3);
#else
  double sum0 = 0.0, sum1 = 0.0;
  double factor = yt / (n_cond_haps * probSumT);
  double tFreq0[4], tFreq1[4];
  for (int h = 0; h < 4; h++) {
    tFreq0[h] = probSumH[h] * factor;
    tFreq1[h] = probSumH[h + 4] * factor;
  }
  double nt_ = nt / probSumT;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    for (int h = 0; h < 4; h++) {
      prob[i + h] = prob[i + h] * nt_ + tFreq0[h];
      prob[i + 4 + h] = prob[i + 4 + h] * nt_ + tFreq1[h];
      sum0 += prob[i + h];
      sum1 += prob[i + 4 + h];
    }
  }
  for (int h = 0; h < 4; h++)
    probSumH[h] = sum0;
  for (int h = 4; h < HAP_NUMBER; h++)
    probSumH[h] = sum1;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline void haplotype_segment_double::COLLAPSE_MIS() {
#if defined(__AVX2__)
  __m256d _sum0 = _mm256_set1_pd(0.0f);
  __m256d _sum1 = _mm256_set1_pd(0.0f);
  __m256d _tFreq = _mm256_set1_pd(yt / n_cond_haps);
  __m256d _nt = _mm256_set1_pd(nt / probSumT);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    __m256d _prob0 = _mm256_set1_pd(probSumK[k]);
    __m256d _prob1 = _mm256_set1_pd(probSumK[k]);
    _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq);
    _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq);
    _sum0 = _mm256_add_pd(_sum0, _prob0);
    _sum1 = _mm256_add_pd(_sum1, _prob1);
    _mm256_store_pd(&prob[i], _prob0);
    _mm256_store_pd(&prob[i + 4], _prob1);
  }
  _mm256_store_pd(&probSumH[0], _sum0);
  _mm256_store_pd(&probSumH[4], _sum1);
#elif defined(__ARM_NEON)
  float64x2_t _sum = vdupq_n_f64(0.0);
  float64x2_t _tFreq = vdupq_n_f64(yt / n_cond_haps);
  float64x2_t _nt = vdupq_n_f64(nt / probSumT);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    float64x2_t _prob = vmlaq_f64(_tFreq, vdupq_n_f64(probSumK[k]), _nt);
    _sum = vaddq_f64(_sum, _prob);
    vst1q_f64(&prob[i], _prob);
    vst1q_f64(&prob[i + 2], _prob);
    vst1q_f64(&prob[i + 4], _prob);
    vst1q_f64(&prob[i + 6], _prob);
  }
  vst1q_f64(&probSumH[0], _sum);
  vst1q_f64(&probSumH[2], _sum);
  vst1q_f64(&probSumH[4], _sum);
  vst1q_f64(&probSumH[6], _sum);
#else
  double sum0 = 0.0, sum1 = 0.0;
  double tFreq = yt / n_cond_haps;
  double nt_ = nt / probSumT;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    double p = probSumK[k] * nt_ + tFreq;
    for (int h = 0; h < 4; h++)
      prob[i + h] = p;
    for (int h = 4; h < HAP_NUMBER; h++)
      prob[i + h] = p;
    sum0 += 4.0 * p;
    sum1 += 4.0 * p;
  }
  for (int h = 0; h < 4; h++)
    probSumH[h] = sum0;
  for (int h = 4; h < HAP_NUMBER; h++)
    probSumH[h] = sum1;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************					SUM Ks
 * ************************/
/*******************************************************************************/

inline void haplotype_segment_double::SUMK() {
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    probSumK[k] = prob[i + 0] + prob[i + 1] + prob[i + 2] + prob[i + 3] +
                  prob[i + 4] + prob[i + 5] + prob[i + 6] + prob[i + 7];
  }
}

/*******************************************************************************/
/*****************		TRANSITION COMPUTATIONS
 * ************************/
/*******************************************************************************/

inline bool haplotype_segment_double::TRANS_HAP() {
  sumHProbs = 0.0f;
  unsigned int curr_rel_segment_index = curr_segment_index - segment_first;
  yt = M.getForwardTransProb(AlphaLocus[curr_rel_segment_index - 1],
                             prev_abs_locus);
  nt = 1.0f - yt;
  double fact1 = nt / AlphaSumSum[curr_rel_segment_index - 1];
  for (int h1 = 0; h1 < HAP_NUMBER; h1++) {
    double fact2 = (AlphaSum[curr_rel_segment_index - 1][h1] /
                    AlphaSumSum[curr_rel_segment_index - 1]) *
                   yt / n_cond_haps;
#ifdef __AVX2__
    __m256d _sum0 = _mm256_set1_pd(0.0f);
    __m256d _sum1 = _mm256_set1_pd(0.0f);
    for (int k = 0; k < n_cond_haps; k++) {
      __m256d _alpha = _mm256_set1_pd(
          Alpha[curr_rel_segment_index - 1][k * HAP_NUMBER + h1] * fact1 +
          fact2);
      __m256d _beta0 = _mm256_load_pd(&prob[k * HAP_NUMBER + 0]);
      __m256d _beta1 = _mm256_load_pd(&prob[k * HAP_NUMBER + 4]);
      _sum0 = _mm256_add_pd(_sum0, _mm256_mul_pd(_alpha, _beta0));
      _sum1 = _mm256_add_pd(_sum1, _mm256_mul_pd(_alpha, _beta1));
    }
    _mm256_store_pd(&HProbs[h1 * HAP_NUMBER + 0], _sum0);
    _mm256_store_pd(&HProbs[h1 * HAP_NUMBER + 4], _sum1);
#elif defined(__ARM_NEON)
    float64x2_t _sum0 = vdupq_n_f64(0.0);
    float64x2_t _sum1 = vdupq_n_f64(0.0);
    float64x2_t _sum2 = vdupq_n_f64(0.0);
    float64x2_t _sum3 = vdupq_n_f64(0.0);
    for (int k = 0; k < n_cond_haps; k++) {
      float64x2_t _alpha = vdupq_n_f64(
          Alpha[curr_rel_segment_index - 1][k * HAP_NUMBER + h1] * fact1 +
          fact2);
      float64x2_t _beta0 = vld1q_f64(&prob[k * HAP_NUMBER + 0]);
      float64x2_t _beta1 = vld1q_f64(&prob[k * HAP_NUMBER + 2]);
      float64x2_t _beta2 = vld1q_f64(&prob[k * HAP_NUMBER + 4]);
      float64x2_t _beta3 = vld1q_f64(&prob[k * HAP_NUMBER + 6]);
      _sum0 = vmlaq_f64(_sum0, _alpha, _beta0);
      _sum1 = vmlaq_f64(_sum1, _alpha, _beta1);
      _sum2 = vmlaq_f64(_sum2, _alpha, _beta2);
      _sum3 = vmlaq_f64(_sum3, _alpha, _beta3);
    }
    vst1q_f64(&HProbs[h1 * HAP_NUMBER + 0], _sum0);
    vst1q_f64(&HProbs[h1 * HAP_NUMBER + 2], _sum1);
    vst1q_f64(&HProbs[h1 * HAP_NUMBER + 4], _sum2);
    vst1q_f64(&HProbs[h1 * HAP_NUMBER + 6], _sum3);
#else
    double sum0 = 0.0, sum1 = 0.0;
    for (int k = 0; k < n_cond_haps; k++) {
      double alpha =
          Alpha[curr_rel_segment_index - 1][k * HAP_NUMBER + h1] * fact1 +
          fact2;
      for (int h = 0; h < 4; h++)
        sum0 += alpha * prob[k * HAP_NUMBER + h];
      for (int h = 4; h < HAP_NUMBER; h++)
        sum1 += alpha * prob[k * HAP_NUMBER + h];
    }
    for (int h = 0; h < 4; h++)
      HProbs[h1 * HAP_NUMBER + h] = sum0;
    for (int h = 4; h < HAP_NUMBER; h++)
      HProbs[h1 * HAP_NUMBER + h] = sum1;
#endif
    sumHProbs += HProbs[h1 * HAP_NUMBER + 0] + HProbs[h1 * HAP_NUMBER + 1] +
                 HProbs[h1 * HAP_NUMBER + 2] + HProbs[h1 * HAP_NUMBER + 3] +
                 HProbs[h1 * HAP_NUMBER + 4] + HProbs[h1 * HAP_NUMBER + 5] +
                 HProbs[h1 * HAP_NUMBER + 6] + HProbs[h1 * HAP_NUMBER + 7];
  }
  return (std::isnan(sumHProbs) || std::isinf(sumHProbs) ||
          sumHProbs < std::numeric_limits<double>::min());
}

inline bool haplotype_segment_double::TRANS_DIP_MULT() {
  sumDProbs = 0.0f;
  double scaling = 1.0 / sumHProbs;
  for (int pd = 0, t = 0; pd < 64; ++pd) {
    if (DIP_GET(G->Diplotypes[curr_segment_index - 1], pd)) {
      for (int nd = 0; nd < 64; ++nd) {
        if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
          DProbs[t] =
              (((double)HProbs[DIP_HAP0(pd) * HAP_NUMBER + DIP_HAP0(nd)]) *
               scaling) *
              ((double)(HProbs[DIP_HAP1(pd) * HAP_NUMBER + DIP_HAP1(nd)]) *
               scaling);
          sumDProbs += DProbs[t];
          t++;
        }
      }
    }
  }
  return (std::isnan(sumDProbs) || std::isinf(sumDProbs) ||
          sumDProbs < std::numeric_limits<double>::min());
}

inline bool haplotype_segment_double::TRANS_DIP_ADD() {
  sumDProbs = 0.0f;
  double scaling = 1.0 / sumHProbs;
  for (int pd = 0, t = 0; pd < 64; ++pd) {
    if (DIP_GET(G->Diplotypes[curr_segment_index - 1], pd)) {
      for (int nd = 0; nd < 64; ++nd) {
        if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
          DProbs[t] = DProbs[t] =
              (((double)HProbs[DIP_HAP0(pd) * HAP_NUMBER + DIP_HAP0(nd)]) *
               scaling) +
              ((double)(HProbs[DIP_HAP1(pd) * HAP_NUMBER + DIP_HAP1(nd)]) *
               scaling);
          sumDProbs += DProbs[t];
          t++;
        }
      }
    }
  }
  return (std::isnan(sumDProbs) || std::isinf(sumDProbs) ||
          sumDProbs < std::numeric_limits<double>::min());
}

inline void
haplotype_segment_double::IMPUTE(std::vector<float> &missing_probabilities) {
#ifdef __AVX2__
  __m256d _sum0 = _mm256_set1_pd(0.0f);
  __m256d _sum1 = _mm256_set1_pd(0.0f);

  __m256d _sumA0[2], _sumA1[2];
  _sumA0[0] = _mm256_set1_pd(0.0f);
  _sumA0[1] = _mm256_set1_pd(0.0f);
  _sumA1[0] = _mm256_set1_pd(0.0f);
  _sumA1[1] = _mm256_set1_pd(0.0f);

  __m256d _alphaSum0 = _mm256_load_pd(&AlphaSumMissing[curr_rel_missing][0]);
  __m256d _alphaSum1 = _mm256_load_pd(&AlphaSumMissing[curr_rel_missing][1]);

  __m256d _ones = _mm256_set1_pd(1.0f);

  _alphaSum0 = _mm256_div_pd(_ones, _alphaSum0);
  _alphaSum1 = _mm256_div_pd(_ones, _alphaSum1);

  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    __m256d _prob0 = _mm256_load_pd(&prob[i]);
    __m256d _prob1 = _mm256_load_pd(&prob[i + 4]);

    __m256d _alpha0 = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i + 0]);
    __m256d _alpha1 = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i + 4]);

    _sum0 = _mm256_mul_pd(_mm256_mul_pd(_alpha0, _alphaSum0), _prob0);
    _sum1 = _mm256_mul_pd(_mm256_mul_pd(_alpha1, _alphaSum1), _prob1);

    _sumA0[ah] = _mm256_add_pd(_sumA0[ah], _sum0);
    _sumA1[ah] = _mm256_add_pd(_sumA1[ah], _sum1);
  }
  double *prob0 = (double *)&_sumA0[0];
  double *prob1 = (double *)&_sumA0[1];
  for (int h = 0; h < 4; h++)
    missing_probabilities[curr_abs_missing * HAP_NUMBER + h] =
        prob1[h] / (prob0[h] + prob1[h]);
  prob0 = (double *)&_sumA1[0];
  prob1 = (double *)&_sumA1[1];
  for (int h = 4; h < HAP_NUMBER; h++)
    missing_probabilities[curr_abs_missing * HAP_NUMBER + h] =
        prob1[h] / (prob0[h] + prob1[h]);
#else
  double sumA0[2][HAP_NUMBER], sumA1[2][HAP_NUMBER];
  for (int a = 0; a < 2; a++)
    for (int h = 0; h < HAP_NUMBER; h++) {
      sumA0[a][h] = 0.0;
      sumA1[a][h] = 0.0;
    }
  double alphaSum0[HAP_NUMBER], alphaSum1[HAP_NUMBER];
  for (int h = 0; h < HAP_NUMBER; h++) {
    alphaSum0[h] = 1.0 / AlphaSumMissing[curr_rel_missing][h];
    alphaSum1[h] = 1.0 / AlphaSumMissing[curr_rel_missing][h + HAP_NUMBER];
  }
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    for (int h = 0; h < HAP_NUMBER; h++) {
      double s0 =
          AlphaMissing[curr_rel_missing][i + h] * alphaSum0[h] * prob[i + h];
      double s1 =
          AlphaMissing[curr_rel_missing][i + h] * alphaSum1[h] * prob[i + h];
      sumA0[ah][h] += s0;
      sumA1[ah][h] += s1;
    }
  }
  for (int h = 0; h < HAP_NUMBER; h++) {
    missing_probabilities[curr_abs_missing * HAP_NUMBER + h] =
        (float)(sumA1[1][h] / (sumA0[0][h] + sumA1[0][h]));
  }
#endif
}

#endif
