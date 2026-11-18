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

#ifndef _HAPLOTYPE_SEGMENT_SINGLE_H
#define _HAPLOTYPE_SEGMENT_SINGLE_H

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

class haplotype_segment_single {
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
  float probSumT;
  aligned_vector32<float> prob;
  aligned_vector32<float> probSumK;
  aligned_vector32<float> probSumH;
  std::vector<aligned_vector32<float>> Alpha;
  std::vector<aligned_vector32<float>> AlphaSum;
  std::vector<int> AlphaLocus;
  aligned_vector32<float> AlphaSumSum;
  std::vector<aligned_vector32<float>> AlphaMissing;
  std::vector<aligned_vector32<float>> AlphaSumMissing;
  float HProbs[HAP_NUMBER * HAP_NUMBER] __attribute__((aligned(32)));
  double DProbs[HAP_NUMBER * HAP_NUMBER * HAP_NUMBER * HAP_NUMBER]
      __attribute__((aligned(32)));

  // STATIC ARRAYS
  float sumHProbs;
  double sumDProbs;
  float g0[HAP_NUMBER], g1[HAP_NUMBER];
  float nt, yt;

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
  haplotype_segment_single(genotype *, bitmatrix &, std::vector<unsigned int> &,
                           window &, hmm_parameters &);
  ~haplotype_segment_single();

  // void fetch();
  void forward();
  int backward(std::vector<double> &, std::vector<float> &);
};

/*******************************************************************************/
/*****************			HOMOZYGOUS GENOTYPE
 * ************************/
/*******************************************************************************/

inline void haplotype_segment_single::INIT_HOM() {
  bool ag =
      VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
#if defined(__AVX2__)
  __m256 _sum = _mm256_set1_ps(0.0f);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);

    // std::cout << curr_rel_locus << " " << curr_rel_locus_offset << " " << k
    // << " " << n_cond_haps << std::endl;

    __m256 _prob = _mm256_set1_ps((ag == ah) ? 1.0f : M.ed / M.ee);
    _sum = _mm256_add_ps(_sum, _prob);
    _mm256_store_ps(&prob[i], _prob);
  }
  _mm256_store_ps(&probSumH[0], _sum);
#elif defined(__ARM_NEON)
  float32x4_t _sum0 = vdupq_n_f32(0.0f);
  float32x4_t _sum1 = vdupq_n_f32(0.0f);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float32x4_t _prob = vdupq_n_f32((ag == ah) ? 1.0f : M.ed / M.ee);
    _sum0 = vaddq_f32(_sum0, _prob);
    _sum1 = vaddq_f32(_sum1, _prob);
    vst1q_f32(&prob[i], _prob);
    vst1q_f32(&prob[i + 4], _prob);
  }
  vst1q_f32(&probSumH[0], _sum0);
  vst1q_f32(&probSumH[4], _sum1);
#else
  float sum = 0.0f;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float p = (ag == ah) ? 1.0f : M.ed / M.ee;
    for (int h = 0; h < HAP_NUMBER; h++)
      prob[i + h] = p;
    sum += p;
  }
  for (int h = 0; h < HAP_NUMBER; h++)
    probSumH[h] = sum;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline bool haplotype_segment_single::RUN_HOM(char rare_allele) {
  bool ag =
      VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
  if (rare_allele < 0 || ag == rare_allele) {
#if defined(__AVX2__)
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
    __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
    _tFreq = _mm256_mul_ps(_tFreq, _factor);
    __m256 _nt = _mm256_set1_ps(nt / probSumT);
    __m256 _mismatch = _mm256_set1_ps(M.ed / M.ee);
    for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
      bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
      __m256 _prob = _mm256_load_ps(&prob[i]);
      _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
      if (ag != ah)
        _prob = _mm256_mul_ps(_prob, _mismatch);
      _sum = _mm256_add_ps(_sum, _prob);
      _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
#elif defined(__ARM_NEON)
    float32x4_t _sum0 = vdupq_n_f32(0.0f);
    float32x4_t _sum1 = vdupq_n_f32(0.0f);
    float32x4_t _factor = vdupq_n_f32(yt / (n_cond_haps * probSumT));
    float32x4_t _tFreq0 = vld1q_f32(&probSumH[0]);
    float32x4_t _tFreq1 = vld1q_f32(&probSumH[4]);
    _tFreq0 = vmulq_f32(_tFreq0, _factor);
    _tFreq1 = vmulq_f32(_tFreq1, _factor);
    float32x4_t _nt = vdupq_n_f32(nt / probSumT);
    float32x4_t _mismatch = vdupq_n_f32(M.ed / M.ee);
    for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
      bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
      float32x4_t _prob0 = vld1q_f32(&prob[i]);
      float32x4_t _prob1 = vld1q_f32(&prob[i + 4]);
      _prob0 = vmlaq_f32(_tFreq0, _prob0, _nt);
      _prob1 = vmlaq_f32(_tFreq1, _prob1, _nt);
      if (ag != ah) {
        _prob0 = vmulq_f32(_prob0, _mismatch);
        _prob1 = vmulq_f32(_prob1, _mismatch);
      }
      _sum0 = vaddq_f32(_sum0, _prob0);
      _sum1 = vaddq_f32(_sum1, _prob1);
      vst1q_f32(&prob[i], _prob0);
      vst1q_f32(&prob[i + 4], _prob1);
    }
    vst1q_f32(&probSumH[0], _sum0);
    vst1q_f32(&probSumH[4], _sum1);
#else
    float sum = 0.0f;
    float factor = yt / (n_cond_haps * probSumT);
    float tFreq[HAP_NUMBER];
    for (int h = 0; h < HAP_NUMBER; h++)
      tFreq[h] = probSumH[h] * factor;
    float nt_ = nt / probSumT;
    float mismatch = M.ed / M.ee;
    for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
      bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
      for (int h = 0; h < HAP_NUMBER; h++) {
        prob[i + h] = prob[i + h] * nt_ + tFreq[h];
        if (ag != ah)
          prob[i + h] *= mismatch;
      }
      for (int h = 0; h < HAP_NUMBER; h++)
        sum += prob[i + h];
    }
    for (int h = 0; h < HAP_NUMBER; h++)
      probSumH[h] = sum;
#endif
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
               probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
    return true;
  }
  return false;
}

inline void haplotype_segment_single::COLLAPSE_HOM() {
  bool ag =
      VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
#if defined(__AVX2__)
  __m256 _sum = _mm256_set1_ps(0.0f);
  __m256 _tFreq =
      _mm256_set1_ps(yt / n_cond_haps); // Check divide by probSumT here!
  __m256 _nt = _mm256_set1_ps(nt / probSumT);
  __m256 _mismatch = _mm256_set1_ps(M.ed / M.ee);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    __m256 _prob = _mm256_set1_ps(probSumK[k]);
    _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
    if (ag != ah)
      _prob = _mm256_mul_ps(_prob, _mismatch);
    _sum = _mm256_add_ps(_sum, _prob);
    _mm256_store_ps(&prob[i], _prob);
  }
  _mm256_store_ps(&probSumH[0], _sum);
#elif defined(__ARM_NEON)
  float32x4_t _sum0 = vdupq_n_f32(0.0f);
  float32x4_t _sum1 = vdupq_n_f32(0.0f);
  float32x4_t _tFreq =
      vdupq_n_f32(yt / n_cond_haps); // Check divide by probSumT here!
  float32x4_t _nt = vdupq_n_f32(nt / probSumT);
  float32x4_t _mismatch = vdupq_n_f32(M.ed / M.ee);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float32x4_t _prob = vdupq_n_f32(probSumK[k]);
    _prob = vmlaq_f32(_tFreq, _prob, _nt);
    if (ag != ah)
      _prob = vmulq_f32(_prob, _mismatch);
    _sum0 = vaddq_f32(_sum0, _prob);
    _sum1 = vaddq_f32(_sum1, _prob);
    vst1q_f32(&prob[i], _prob);
    vst1q_f32(&prob[i + 4], _prob);
  }
  vst1q_f32(&probSumH[0], _sum0);
  vst1q_f32(&probSumH[4], _sum1);
#else
  float sum = 0.0f;
  float tFreq = yt / n_cond_haps; // Check divide by probSumT here!
  float nt_ = nt / probSumT;
  float mismatch = M.ed / M.ee;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float p = probSumK[k] * nt_ + tFreq;
    if (ag != ah)
      p *= mismatch;
    for (int h = 0; h < HAP_NUMBER; h++)
      prob[i + h] = p;
    sum += p;
  }
  for (int h = 0; h < HAP_NUMBER; h++)
    probSumH[h] = sum;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************			HETEROZYGOUS GENOTYPE
 * ********************/
/*******************************************************************************/

inline void haplotype_segment_single::INIT_AMB() {
  unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
  for (int h = 0; h < HAP_NUMBER; h++) {
    g0[h] = HAP_GET(amb_code, h) ? M.ed / M.ee : 1.0f;
    g1[h] = HAP_GET(amb_code, h) ? 1.0f : M.ed / M.ee;
  }
#if defined(__AVX2__)
  __m256 _sum = _mm256_set1_ps(0.0f);
  __m256 _emit[2];
  _emit[0] = _mm256_loadu_ps(&g0[0]);
  _emit[1] = _mm256_loadu_ps(&g1[0]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    __m256 _prob = _emit[ah];
    _sum = _mm256_add_ps(_sum, _prob);
    _mm256_store_ps(&prob[i], _prob);
  }
  _mm256_store_ps(&probSumH[0], _sum);
#elif defined(__ARM_NEON)
  float32x4_t _sum0 = vdupq_n_f32(0.0f);
  float32x4_t _sum1 = vdupq_n_f32(0.0f);
  float32x4_t _e00 = vld1q_f32(&g0[0]), _e01 = vld1q_f32(&g1[0]);
  float32x4_t _e10 = vld1q_f32(&g0[4]), _e11 = vld1q_f32(&g1[4]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float32x4_t _p0 = ah ? _e01 : _e00;
    float32x4_t _p1 = ah ? _e11 : _e10;
    _sum0 = vaddq_f32(_sum0, _p0);
    _sum1 = vaddq_f32(_sum1, _p1);
    vst1q_f32(&prob[i], _p0);
    vst1q_f32(&prob[i + 4], _p1);
  }
  vst1q_f32(&probSumH[0], _sum0);
  vst1q_f32(&probSumH[4], _sum1);
#else
  float sum = 0.0f;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    for (int h = 0; h < HAP_NUMBER; h++)
      prob[i + h] = ah ? g1[h] : g0[h];
    for (int h = 0; h < HAP_NUMBER; h++)
      sum += ah ? g1[h] : g0[h];
  }
  for (int h = 0; h < HAP_NUMBER; h++)
    probSumH[h] = sum;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline void haplotype_segment_single::RUN_AMB() {
  unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
  for (int h = 0; h < HAP_NUMBER; h++) {
    g0[h] = HAP_GET(amb_code, h) ? M.ed / M.ee : 1.0f;
    g1[h] = HAP_GET(amb_code, h) ? 1.0f : M.ed / M.ee;
  }
#if defined(__AVX2__)
  __m256 _sum = _mm256_set1_ps(0.0f);
  __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
  __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
  _tFreq = _mm256_mul_ps(_tFreq, _factor);
  __m256 _nt = _mm256_set1_ps(nt / probSumT);
  __m256 _emit[2];
  _emit[0] = _mm256_loadu_ps(&g0[0]);
  _emit[1] = _mm256_loadu_ps(&g1[0]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    __m256 _prob = _mm256_load_ps(&prob[i]);
    _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
    _prob = _mm256_mul_ps(_prob, _emit[ah]);
    _sum = _mm256_add_ps(_sum, _prob);
    _mm256_store_ps(&prob[i], _prob);
  }
  _mm256_store_ps(&probSumH[0], _sum);
#elif defined(__ARM_NEON)
  float32x4_t _sum0 = vdupq_n_f32(0.0f);
  float32x4_t _sum1 = vdupq_n_f32(0.0f);
  float32x4_t _factor = vdupq_n_f32(yt / (n_cond_haps * probSumT));
  float32x4_t _tFreq0 = vld1q_f32(&probSumH[0]);
  float32x4_t _tFreq1 = vld1q_f32(&probSumH[4]);
  _tFreq0 = vmulq_f32(_tFreq0, _factor);
  _tFreq1 = vmulq_f32(_tFreq1, _factor);
  float32x4_t _nt = vdupq_n_f32(nt / probSumT);
  float32x4_t _e00 = vld1q_f32(&g0[0]), _e01 = vld1q_f32(&g1[0]);
  float32x4_t _e10 = vld1q_f32(&g0[4]), _e11 = vld1q_f32(&g1[4]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float32x4_t _p0 = vld1q_f32(&prob[i]);
    float32x4_t _p1 = vld1q_f32(&prob[i + 4]);
    _p0 = vmulq_f32(vmlaq_f32(_tFreq0, _p0, _nt), ah ? _e01 : _e00);
    _p1 = vmulq_f32(vmlaq_f32(_tFreq1, _p1, _nt), ah ? _e11 : _e10);
    _sum0 = vaddq_f32(_sum0, _p0);
    _sum1 = vaddq_f32(_sum1, _p1);
    vst1q_f32(&prob[i], _p0);
    vst1q_f32(&prob[i + 4], _p1);
  }
  vst1q_f32(&probSumH[0], _sum0);
  vst1q_f32(&probSumH[4], _sum1);
#else
  float sum = 0.0f;
  float factor = yt / (n_cond_haps * probSumT);
  float tFreq[HAP_NUMBER];
  for (int h = 0; h < HAP_NUMBER; h++)
    tFreq[h] = probSumH[h] * factor;
  float nt_ = nt / probSumT;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    for (int h = 0; h < HAP_NUMBER; h++) {
      prob[i + h] = (prob[i + h] * nt_ + tFreq[h]) * (ah ? g1[h] : g0[h]);
    }
    for (int h = 0; h < HAP_NUMBER; h++)
      sum += prob[i + h];
  }
  for (int h = 0; h < HAP_NUMBER; h++)
    probSumH[h] = sum;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*
inline
void haplotype_segment_single::RUN_AMB() {
        unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
        for (int h = 0 ; h < HAP_NUMBER ; h ++) {
                g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
                g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
        }
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
        for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
                bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
                __m256d _prob0 = _mm256_load_pd(&prob[i+0]);
                __m256d _prob1 = _mm256_load_pd(&prob[i+4]);
                _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
                _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
                _prob0 = _mm256_mul_pd(_prob0, _emit0[ah]);
                _prob1 = _mm256_mul_pd(_prob1, _emit1[ah]);
                _sum0 = _mm256_add_pd(_sum0, _prob0);
                _sum1 = _mm256_add_pd(_sum1, _prob1);
                _mm256_store_pd(&prob[i+0], _prob0);
                _mm256_store_pd(&prob[i+4], _prob1);
        }
        _mm256_store_pd(&probSumH[0], _sum0);
        _mm256_store_pd(&probSumH[4], _sum1);
        probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
*/

inline void haplotype_segment_single::COLLAPSE_AMB() {
  unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
  for (int h = 0; h < HAP_NUMBER; h++) {
    g0[h] = HAP_GET(amb_code, h) ? M.ed / M.ee : 1.0f;
    g1[h] = HAP_GET(amb_code, h) ? 1.0f : M.ed / M.ee;
  }
#if defined(__AVX2__)
  __m256 _sum = _mm256_set1_ps(0.0f);
  __m256 _tFreq = _mm256_set1_ps(yt / n_cond_haps);
  __m256 _nt = _mm256_set1_ps(nt / probSumT);
  __m256 _emit[2];
  _emit[0] = _mm256_loadu_ps(&g0[0]);
  _emit[1] = _mm256_loadu_ps(&g1[0]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    __m256 _prob = _mm256_set1_ps(probSumK[k]);
    _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
    _prob = _mm256_mul_ps(_prob, _emit[ah]);
    _sum = _mm256_add_ps(_sum, _prob);
    _mm256_store_ps(&prob[i], _prob);
  }
  _mm256_store_ps(&probSumH[0], _sum);
#elif defined(__ARM_NEON)
  float32x4_t _sum0 = vdupq_n_f32(0.0f);
  float32x4_t _sum1 = vdupq_n_f32(0.0f);
  float32x4_t _tFreq = vdupq_n_f32(yt / n_cond_haps);
  float32x4_t _nt = vdupq_n_f32(nt / probSumT);
  float32x4_t _e00 = vld1q_f32(&g0[0]), _e01 = vld1q_f32(&g1[0]);
  float32x4_t _e10 = vld1q_f32(&g0[4]), _e11 = vld1q_f32(&g1[4]);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float32x4_t _p0 = vmulq_f32(
        vmlaq_f32(_tFreq, vdupq_n_f32(probSumK[k]), _nt), ah ? _e01 : _e00);
    float32x4_t _p1 = vmulq_f32(
        vmlaq_f32(_tFreq, vdupq_n_f32(probSumK[k]), _nt), ah ? _e11 : _e10);
    _sum0 = vaddq_f32(_sum0, _p0);
    _sum1 = vaddq_f32(_sum1, _p1);
    vst1q_f32(&prob[i], _p0);
    vst1q_f32(&prob[i + 4], _p1);
  }
  vst1q_f32(&probSumH[0], _sum0);
  vst1q_f32(&probSumH[4], _sum1);
#else
  float sum = 0.0f;
  float tFreq = yt / n_cond_haps;
  float nt_ = nt / probSumT;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float p = probSumK[k] * nt_ + tFreq;
    for (int h = 0; h < HAP_NUMBER; h++)
      prob[i + h] = p * (ah ? g1[h] : g0[h]);
    for (int h = 0; h < HAP_NUMBER; h++)
      sum += prob[i + h];
  }
  for (int h = 0; h < HAP_NUMBER; h++)
    probSumH[h] = sum;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************			MISSING GENOTYPE
 * ************************/
/*******************************************************************************/

inline void haplotype_segment_single::INIT_MIS() {
  fill(prob.begin(), prob.end(), 1.0f / (HAP_NUMBER * n_cond_haps));
  fill(probSumH.begin(), probSumH.end(), 1.0f / HAP_NUMBER);
  probSumT = 1.0f;
}

inline void haplotype_segment_single::RUN_MIS() {
#if defined(__AVX2__)
  __m256 _sum = _mm256_set1_ps(0.0f);
  __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
  __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
  _tFreq = _mm256_mul_ps(_tFreq, _factor);
  __m256 _nt = _mm256_set1_ps(nt / probSumT);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    __m256 _prob = _mm256_load_ps(&prob[i]);
    _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
    _sum = _mm256_add_ps(_sum, _prob);
    _mm256_store_ps(&prob[i], _prob);
  }
  _mm256_store_ps(&probSumH[0], _sum);
#elif defined(__ARM_NEON)
  float32x4_t _sum0 = vdupq_n_f32(0.0f);
  float32x4_t _sum1 = vdupq_n_f32(0.0f);
  float32x4_t _factor = vdupq_n_f32(yt / (n_cond_haps * probSumT));
  float32x4_t _tFreq0 = vld1q_f32(&probSumH[0]);
  float32x4_t _tFreq1 = vld1q_f32(&probSumH[4]);
  _tFreq0 = vmulq_f32(_tFreq0, _factor);
  _tFreq1 = vmulq_f32(_tFreq1, _factor);
  float32x4_t _nt = vdupq_n_f32(nt / probSumT);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    float32x4_t _p0 = vld1q_f32(&prob[i]);
    float32x4_t _p1 = vld1q_f32(&prob[i + 4]);
    _p0 = vmlaq_f32(_tFreq0, _p0, _nt);
    _p1 = vmlaq_f32(_tFreq1, _p1, _nt);
    _sum0 = vaddq_f32(_sum0, _p0);
    _sum1 = vaddq_f32(_sum1, _p1);
    vst1q_f32(&prob[i], _p0);
    vst1q_f32(&prob[i + 4], _p1);
  }
  vst1q_f32(&probSumH[0], _sum0);
  vst1q_f32(&probSumH[4], _sum1);
#else
  float sum = 0.0f;
  float factor = yt / (n_cond_haps * probSumT);
  float tFreq[HAP_NUMBER];
  for (int h = 0; h < HAP_NUMBER; h++)
    tFreq[h] = probSumH[h] * factor;
  float nt_ = nt / probSumT;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    for (int h = 0; h < HAP_NUMBER; h++)
      prob[i + h] = prob[i + h] * nt_ + tFreq[h];
    for (int h = 0; h < HAP_NUMBER; h++)
      sum += prob[i + h];
  }
  for (int h = 0; h < HAP_NUMBER; h++)
    probSumH[h] = sum;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline void haplotype_segment_single::COLLAPSE_MIS() {
#if defined(__AVX2__)
  __m256 _sum = _mm256_set1_ps(0.0f);
  __m256 _tFreq = _mm256_set1_ps(yt / n_cond_haps);
  __m256 _nt = _mm256_set1_ps(nt / probSumT);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    __m256 _prob = _mm256_set1_ps(probSumK[k]);
    _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
    _sum = _mm256_add_ps(_sum, _prob);
    _mm256_store_ps(&prob[i], _prob);
  }
  _mm256_store_ps(&probSumH[0], _sum);
#elif defined(__ARM_NEON)
  float32x4_t _sum0 = vdupq_n_f32(0.0f);
  float32x4_t _sum1 = vdupq_n_f32(0.0f);
  float32x4_t _tFreq = vdupq_n_f32(yt / n_cond_haps);
  float32x4_t _nt = vdupq_n_f32(nt / probSumT);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    float32x4_t _prob = vmlaq_f32(_tFreq, vdupq_n_f32(probSumK[k]), _nt);
    _sum0 = vaddq_f32(_sum0, _prob);
    _sum1 = vaddq_f32(_sum1, _prob);
    vst1q_f32(&prob[i], _prob);
    vst1q_f32(&prob[i + 4], _prob);
  }
  vst1q_f32(&probSumH[0], _sum0);
  vst1q_f32(&probSumH[4], _sum1);
#else
  float sum = 0.0f;
  float tFreq = yt / n_cond_haps;
  float nt_ = nt / probSumT;
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    float p = probSumK[k] * nt_ + tFreq;
    for (int h = 0; h < HAP_NUMBER; h++)
      prob[i + h] = p;
    sum += p;
  }
  for (int h = 0; h < HAP_NUMBER; h++)
    probSumH[h] = sum;
#endif
  probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] +
             probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************					SUM Ks
 * ************************/
/*******************************************************************************/

inline void haplotype_segment_single::SUMK() {
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    probSumK[k] = prob[i + 0] + prob[i + 1] + prob[i + 2] + prob[i + 3] +
                  prob[i + 4] + prob[i + 5] + prob[i + 6] + prob[i + 7];
  }
}

/*******************************************************************************/
/*****************		TRANSITION COMPUTATIONS
 * ************************/
/*******************************************************************************/

inline bool haplotype_segment_single::TRANS_HAP() {
  sumHProbs = 0.0f;
  unsigned int curr_rel_segment_index = curr_segment_index - segment_first;
  yt = M.getForwardTransProb(AlphaLocus[curr_rel_segment_index - 1],
                             prev_abs_locus);
  nt = 1.0f - yt;
  float fact1 = nt / AlphaSumSum[curr_rel_segment_index - 1];
  for (int h1 = 0; h1 < HAP_NUMBER; h1++) {
    float fact2 = (AlphaSum[curr_rel_segment_index - 1][h1] /
                   AlphaSumSum[curr_rel_segment_index - 1]) *
                  yt / n_cond_haps;
#if defined(__AVX2__)
    __m256 _sum = _mm256_set1_ps(0.0f);
    for (int k = 0; k < n_cond_haps; k++) {
      __m256 _alpha = _mm256_set1_ps(
          Alpha[curr_rel_segment_index - 1][k * HAP_NUMBER + h1] * fact1 +
          fact2);
      __m256 _beta = _mm256_load_ps(&prob[k * HAP_NUMBER]);
      _sum = _mm256_add_ps(_sum, _mm256_mul_ps(_alpha, _beta));
    }
    _mm256_store_ps(&HProbs[h1 * HAP_NUMBER], _sum);
#elif defined(__ARM_NEON)
    float32x4_t _sum0 = vdupq_n_f32(0.0f);
    float32x4_t _sum1 = vdupq_n_f32(0.0f);
    for (int k = 0; k < n_cond_haps; k++) {
      float32x4_t _alpha = vdupq_n_f32(
          Alpha[curr_rel_segment_index - 1][k * HAP_NUMBER + h1] * fact1 +
          fact2);
      float32x4_t _b0 = vld1q_f32(&prob[k * HAP_NUMBER]);
      float32x4_t _b1 = vld1q_f32(&prob[k * HAP_NUMBER + 4]);
      _sum0 = vmlaq_f32(_sum0, _alpha, _b0);
      _sum1 = vmlaq_f32(_sum1, _alpha, _b1);
    }
    vst1q_f32(&HProbs[h1 * HAP_NUMBER], _sum0);
    vst1q_f32(&HProbs[h1 * HAP_NUMBER + 4], _sum1);
#else
    float sum = 0.0f;
    for (int k = 0; k < n_cond_haps; k++) {
      float alpha =
          Alpha[curr_rel_segment_index - 1][k * HAP_NUMBER + h1] * fact1 +
          fact2;
      for (int h = 0; h < HAP_NUMBER; h++)
        sum += alpha * prob[k * HAP_NUMBER + h];
    }
    for (int h = 0; h < HAP_NUMBER; h++)
      HProbs[h1 * HAP_NUMBER + h] = sum;
#endif
    sumHProbs += HProbs[h1 * HAP_NUMBER + 0] + HProbs[h1 * HAP_NUMBER + 1] +
                 HProbs[h1 * HAP_NUMBER + 2] + HProbs[h1 * HAP_NUMBER + 3] +
                 HProbs[h1 * HAP_NUMBER + 4] + HProbs[h1 * HAP_NUMBER + 5] +
                 HProbs[h1 * HAP_NUMBER + 6] + HProbs[h1 * HAP_NUMBER + 7];
  }
  return (std::isnan(sumHProbs) || std::isinf(sumHProbs) ||
          sumHProbs < std::numeric_limits<float>::min());
}

inline bool haplotype_segment_single::TRANS_DIP_MULT() {
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

inline bool haplotype_segment_single::TRANS_DIP_ADD() {
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
haplotype_segment_single::IMPUTE(std::vector<float> &missing_probabilities) {
#if defined(__AVX2__)
  __m256 _sum = _mm256_set1_ps(0.0f);
  __m256 _sumA[2];
  _sumA[0] = _mm256_set1_ps(0.0f);
  _sumA[1] = _mm256_set1_ps(0.0f);
  __m256 _alphaSum = _mm256_load_ps(&AlphaSumMissing[curr_rel_missing][0]);
  __m256 _ones = _mm256_set1_ps(1.0f);
  _alphaSum = _mm256_div_ps(_ones, _alphaSum);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    __m256 _prob = _mm256_load_ps(&prob[i]);
    __m256 _alpha = _mm256_load_ps(&AlphaMissing[curr_rel_missing][i]);
    _sum = _mm256_mul_ps(_mm256_mul_ps(_alpha, _alphaSum), _prob);
    _sumA[ah] = _mm256_add_ps(_sumA[ah], _sum);
  }
  float *prob0 = (float *)&_sumA[0];
  float *prob1 = (float *)&_sumA[1];
  for (int h = 0; h < HAP_NUMBER; h++) {
    missing_probabilities[curr_abs_missing * HAP_NUMBER + h] =
        prob1[h] / (prob0[h] + prob1[h]);
  }
#elif defined(__ARM_NEON)
  float32x4_t _sumA0[2], _sumA1[2];
  _sumA0[0] = vdupq_n_f32(0.0f);
  _sumA0[1] = vdupq_n_f32(0.0f);
  _sumA1[0] = vdupq_n_f32(0.0f);
  _sumA1[1] = vdupq_n_f32(0.0f);
  float32x4_t _as0 = vld1q_f32(&AlphaSumMissing[curr_rel_missing][0]);
  float32x4_t _as1 = vld1q_f32(&AlphaSumMissing[curr_rel_missing][4]);
  _as0 = vdivq_f32(vdupq_n_f32(1.0f), _as0);
  _as1 = vdivq_f32(vdupq_n_f32(1.0f), _as1);
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    float32x4_t _prob0 = vld1q_f32(&prob[i]);
    float32x4_t _prob1 = vld1q_f32(&prob[i + 4]);
    float32x4_t _a0 = vld1q_f32(&AlphaMissing[curr_rel_missing][i]);
    float32x4_t _a1 = vld1q_f32(&AlphaMissing[curr_rel_missing][i + 4]);
    float32x4_t _s0 = vmulq_f32(vmulq_f32(_a0, _as0), _prob0);
    float32x4_t _s1 = vmulq_f32(vmulq_f32(_a1, _as1), _prob1);
    _sumA0[ah] = vaddq_f32(_sumA0[ah], _s0);
    _sumA1[ah] = vaddq_f32(_sumA1[ah], _s1);
  }
  float *prob0 = (float *)&_sumA0[0];
  float *prob1 = (float *)&_sumA0[1];
  for (int h = 0; h < 4; h++)
    missing_probabilities[curr_abs_missing * HAP_NUMBER + h] =
        prob1[h] / (prob0[h] + prob1[h]);
  prob0 = (float *)&_sumA1[0];
  prob1 = (float *)&_sumA1[1];
  for (int h = 4; h < HAP_NUMBER; h++)
    missing_probabilities[curr_abs_missing * HAP_NUMBER + h] =
        prob1[h] / (prob0[h] + prob1[h]);
#else
  float sumA[2][HAP_NUMBER];
  for (int a = 0; a < 2; a++)
    for (int h = 0; h < HAP_NUMBER; h++)
      sumA[a][h] = 0.0f;
  float alphaSum[HAP_NUMBER];
  for (int h = 0; h < HAP_NUMBER; h++)
    alphaSum[h] = 1.0f / AlphaSumMissing[curr_rel_missing][h];
  for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
    bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
    for (int h = 0; h < HAP_NUMBER; h++) {
      float s =
          AlphaMissing[curr_rel_missing][i + h] * alphaSum[h] * prob[i + h];
      sumA[ah][h] += s;
    }
  }
  for (int h = 0; h < HAP_NUMBER; h++) {
    missing_probabilities[curr_abs_missing * HAP_NUMBER + h] =
        sumA[1][h] / (sumA[0][h] + sumA[1][h]);
  }
#endif
}

#endif
