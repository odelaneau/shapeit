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

#ifndef _HAPLOTYPE_SEGMENT_H
#define _HAPLOTYPE_SEGMENT_H

#include <utils/otools.h>
#include <objects/compute_job.h>
#include <objects/hmm_parameters.h>

#include <immintrin.h>
#include <type_traits>
#include <boost/align/aligned_allocator.hpp>

template <typename T>
using aligned_vector32 = std::vector<T, boost::alignment::aligned_allocator < T, 32 > >;

/*******************************************************************************/
/***  simd8<T>: HAP_NUMBER (8) lanes of T.                                    */
/***  float  -> one __m256 ; double -> a pair of __m256d (lanes 0-3, 4-7).    */
/***  Thin by-value wrapper; at -O3 it lowers to the same intrinsics as the   */
/***  hand-written code, so the float (default) path is unchanged.            */
/*******************************************************************************/

template <typename T> struct simd8;

template <> struct simd8<float> {
	__m256 v;
	static inline simd8 splat(float x) { return { _mm256_set1_ps(x) }; }
	static inline simd8 load (const float * p) { return { _mm256_load_ps(p) }; }
	static inline simd8 loadu(const float * p) { return { _mm256_loadu_ps(p) }; }
	inline void store(float * p) const { _mm256_store_ps(p, v); }
	friend inline simd8 fmadd(simd8 a, simd8 b, simd8 c) { return { _mm256_fmadd_ps(a.v, b.v, c.v) }; }
	friend inline simd8 operator*(simd8 a, simd8 b) { return { _mm256_mul_ps(a.v, b.v) }; }
	friend inline simd8 operator+(simd8 a, simd8 b) { return { _mm256_add_ps(a.v, b.v) }; }
	friend inline simd8 operator/(simd8 a, simd8 b) { return { _mm256_div_ps(a.v, b.v) }; }
};

template <> struct simd8<double> {
	__m256d lo, hi;
	static inline simd8 splat(double x) { return { _mm256_set1_pd(x), _mm256_set1_pd(x) }; }
	static inline simd8 load (const double * p) { return { _mm256_load_pd(p), _mm256_load_pd(p + 4) }; }
	static inline simd8 loadu(const double * p) { return { _mm256_loadu_pd(p), _mm256_loadu_pd(p + 4) }; }
	inline void store(double * p) const { _mm256_store_pd(p, lo); _mm256_store_pd(p + 4, hi); }
	friend inline simd8 fmadd(simd8 a, simd8 b, simd8 c) { return { _mm256_fmadd_pd(a.lo, b.lo, c.lo), _mm256_fmadd_pd(a.hi, b.hi, c.hi) }; }
	friend inline simd8 operator*(simd8 a, simd8 b) { return { _mm256_mul_pd(a.lo, b.lo), _mm256_mul_pd(a.hi, b.hi) }; }
	friend inline simd8 operator+(simd8 a, simd8 b) { return { _mm256_add_pd(a.lo, b.lo), _mm256_add_pd(a.hi, b.hi) }; }
	friend inline simd8 operator/(simd8 a, simd8 b) { return { _mm256_div_pd(a.lo, b.lo), _mm256_div_pd(a.hi, b.hi) }; }
};

//Left-to-right scalar reduction of the 8-lane probSumH; keeps the exact
//summation order the hand-written kernels used (float path stays bit-exact).
template <typename T> inline
T hsum8(const aligned_vector32<T> & a) {
	return a[0] + a[1] + a[2] + a[3] + a[4] + a[5] + a[6] + a[7];
}

//B1: bit lookup on a pre-hoisted matrix row pointer. Identical result to
//bitmatrix::get(row, col) with col=k, but the row base is computed once by the
//caller. The `__restrict` on the caller's row pointer lets the compiler keep the
//shared byte across the unrolled k's despite the float store in the same loop.
static inline bool hget(const unsigned char * __restrict r, unsigned int k) {
	return (r[k>>3] >> (7 - (k&7))) & 1;
}

template <typename T>
class haplotype_segment {
private:
	//EXTERNAL DATA
	hmm_parameters & M;
	genotype * G;
	bitmatrix Hhap, Hvar;

	//COORDINATES & CONSTANTS
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

	//CURSORS
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


	//DYNAMIC ARRAYS
	T probSumT;
	aligned_vector32 < T > prob;
	aligned_vector32 < T > probSumK;
	aligned_vector32 < T > probSumH;
	std::vector < aligned_vector32 < T > > Alpha;
	std::vector < aligned_vector32 < T > > AlphaSum;
	std::vector < int > AlphaLocus;
	aligned_vector32 < T > AlphaSumSum;
	std::vector < aligned_vector32 < T > > AlphaMissing;
	std::vector < aligned_vector32 < T > > AlphaSumMissing;
	T HProbs [HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));
	double DProbs [HAP_NUMBER * HAP_NUMBER * HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));

	//STATIC ARRAYS
	T sumHProbs;
	double sumDProbs;
	T g0[HAP_NUMBER], g1[HAP_NUMBER];
	T nt, yt;

	//INLINED AND UNROLLED ROUTINES
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
	void IMPUTE(std::vector < float > & );
	bool TRANS_HAP();
	bool TRANS_DIP_MULT();
	bool TRANS_DIP_ADD();
	void SET_FIRST_TRANS(std::vector < double > & );
	int SET_OTHER_TRANS(std::vector < double > & );

public:
	//CONSTRUCTOR/DESTRUCTOR
	haplotype_segment(genotype *, bitmatrix &, std::vector < unsigned int > &, window &, hmm_parameters &);
	~haplotype_segment();

	//void fetch();
	void forward();
	int backward(std::vector < double > &, std::vector < float > &);
};

/*******************************************************************************/
/*****************			HOMOZYGOUS GENOTYPE			************************/
/*******************************************************************************/

template <typename T> inline
void haplotype_segment<T>::INIT_HOM() {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	const unsigned char * __restrict hrow = Hvar.row_ptr(curr_rel_locus+curr_rel_locus_offset);
	simd8<T> _sum = simd8<T>::splat(0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = hget(hrow, k);
		simd8<T> _prob = simd8<T>::splat((ag==ah)?1.0f:M.ed/M.ee);
		_sum = _sum + _prob;
		_prob.store(&prob[i]);
	}
	_sum.store(&probSumH[0]);
	probSumT = hsum8(probSumH);
}

template <typename T> inline
bool haplotype_segment<T>::RUN_HOM(char rare_allele) {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	if (rare_allele < 0 || ag == rare_allele) {
		const unsigned char * __restrict hrow = Hvar.row_ptr(curr_rel_locus+curr_rel_locus_offset);
		simd8<T> _sum0 = simd8<T>::splat(0.0f), _sum1 = simd8<T>::splat(0.0f);
		simd8<T> _sum2 = simd8<T>::splat(0.0f), _sum3 = simd8<T>::splat(0.0f);
		simd8<T> _tFreq = simd8<T>::load(&probSumH[0]) * simd8<T>::splat(yt / (n_cond_haps * probSumT));
		simd8<T> _nt = simd8<T>::splat(nt / probSumT);
		simd8<T> _emit[2] = { simd8<T>::splat(1.0f), simd8<T>::splat(M.ed/M.ee) };
		unsigned int n4 = n_cond_haps & ~3u, k = 0; int i = 0;
		for( ; k < n4 ; k += 4, i += 4*HAP_NUMBER) {
			bool ah0 = hget(hrow, k+0);
			bool ah1 = hget(hrow, k+1);
			bool ah2 = hget(hrow, k+2);
			bool ah3 = hget(hrow, k+3);
			simd8<T> _p0 = fmadd(simd8<T>::load(&prob[i+0*HAP_NUMBER]), _nt, _tFreq) * _emit[ag!=ah0];
			simd8<T> _p1 = fmadd(simd8<T>::load(&prob[i+1*HAP_NUMBER]), _nt, _tFreq) * _emit[ag!=ah1];
			simd8<T> _p2 = fmadd(simd8<T>::load(&prob[i+2*HAP_NUMBER]), _nt, _tFreq) * _emit[ag!=ah2];
			simd8<T> _p3 = fmadd(simd8<T>::load(&prob[i+3*HAP_NUMBER]), _nt, _tFreq) * _emit[ag!=ah3];
			_sum0 = _sum0 + _p0; _sum1 = _sum1 + _p1; _sum2 = _sum2 + _p2; _sum3 = _sum3 + _p3;
			_p0.store(&prob[i+0*HAP_NUMBER]); _p1.store(&prob[i+1*HAP_NUMBER]);
			_p2.store(&prob[i+2*HAP_NUMBER]); _p3.store(&prob[i+3*HAP_NUMBER]);
		}
		for( ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			bool ah = hget(hrow, k);
			simd8<T> _p = fmadd(simd8<T>::load(&prob[i]), _nt, _tFreq) * _emit[ag!=ah];
			_sum0 = _sum0 + _p;
			_p.store(&prob[i]);
		}
		simd8<T> _sum = (_sum0 + _sum1) + (_sum2 + _sum3);
		_sum.store(&probSumH[0]);
		probSumT = hsum8(probSumH);
		return true;
	}
	return false;
}

template <typename T> inline
void haplotype_segment<T>::COLLAPSE_HOM() {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	const unsigned char * __restrict hrow = Hvar.row_ptr(curr_rel_locus+curr_rel_locus_offset);
	simd8<T> _sum0 = simd8<T>::splat(0.0f), _sum1 = simd8<T>::splat(0.0f);
	simd8<T> _sum2 = simd8<T>::splat(0.0f), _sum3 = simd8<T>::splat(0.0f);
	simd8<T> _tFreq = simd8<T>::splat(yt / n_cond_haps);						//Check divide by probSumT here!
	simd8<T> _nt = simd8<T>::splat(nt / probSumT);
	simd8<T> _emit[2] = { simd8<T>::splat(1.0f), simd8<T>::splat(M.ed/M.ee) };
	unsigned int n4 = n_cond_haps & ~3u, k = 0; int i = 0;
	for( ; k < n4 ; k += 4, i += 4*HAP_NUMBER) {
		bool ah0 = hget(hrow, k+0);
		bool ah1 = hget(hrow, k+1);
		bool ah2 = hget(hrow, k+2);
		bool ah3 = hget(hrow, k+3);
		simd8<T> _p0 = fmadd(simd8<T>::splat(probSumK[k+0]), _nt, _tFreq) * _emit[ag!=ah0];
		simd8<T> _p1 = fmadd(simd8<T>::splat(probSumK[k+1]), _nt, _tFreq) * _emit[ag!=ah1];
		simd8<T> _p2 = fmadd(simd8<T>::splat(probSumK[k+2]), _nt, _tFreq) * _emit[ag!=ah2];
		simd8<T> _p3 = fmadd(simd8<T>::splat(probSumK[k+3]), _nt, _tFreq) * _emit[ag!=ah3];
		_sum0 = _sum0 + _p0; _sum1 = _sum1 + _p1; _sum2 = _sum2 + _p2; _sum3 = _sum3 + _p3;
		_p0.store(&prob[i+0*HAP_NUMBER]); _p1.store(&prob[i+1*HAP_NUMBER]);
		_p2.store(&prob[i+2*HAP_NUMBER]); _p3.store(&prob[i+3*HAP_NUMBER]);
	}
	for( ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = hget(hrow, k);
		simd8<T> _p = fmadd(simd8<T>::splat(probSumK[k]), _nt, _tFreq) * _emit[ag!=ah];
		_sum0 = _sum0 + _p;
		_p.store(&prob[i]);
	}
	simd8<T> _sum = (_sum0 + _sum1) + (_sum2 + _sum3);
	_sum.store(&probSumH[0]);
	probSumT = hsum8(probSumH);
}

/*******************************************************************************/
/*****************			HETEROZYGOUS GENOTYPE			********************/
/*******************************************************************************/

template <typename T> inline
void haplotype_segment<T>::INIT_AMB() {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	const unsigned char * __restrict hrow = Hvar.row_ptr(curr_rel_locus+curr_rel_locus_offset);
	simd8<T> _emit[2] = { simd8<T>::loadu(&g0[0]), simd8<T>::loadu(&g1[0]) };
	simd8<T> _sum = simd8<T>::splat(0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = hget(hrow, k);
		simd8<T> _prob = _emit[ah];
		_sum = _sum + _prob;
		_prob.store(&prob[i]);
	}
	_sum.store(&probSumH[0]);
	probSumT = hsum8(probSumH);
}

template <typename T> inline
void haplotype_segment<T>::RUN_AMB() {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	simd8<T> _sum0 = simd8<T>::splat(0.0f), _sum1 = simd8<T>::splat(0.0f);
	simd8<T> _sum2 = simd8<T>::splat(0.0f), _sum3 = simd8<T>::splat(0.0f);
	simd8<T> _tFreq = simd8<T>::load(&probSumH[0]) * simd8<T>::splat(yt / (n_cond_haps * probSumT));
	simd8<T> _nt = simd8<T>::splat(nt / probSumT);
	simd8<T> _emit[2] = { simd8<T>::loadu(&g0[0]), simd8<T>::loadu(&g1[0]) };
	const unsigned char * __restrict hrow = Hvar.row_ptr(curr_rel_locus+curr_rel_locus_offset);
	unsigned int n4 = n_cond_haps & ~3u, k = 0; int i = 0;
	for( ; k < n4 ; k += 4, i += 4*HAP_NUMBER) {
		bool ah0 = hget(hrow, k+0);
		bool ah1 = hget(hrow, k+1);
		bool ah2 = hget(hrow, k+2);
		bool ah3 = hget(hrow, k+3);
		simd8<T> _p0 = fmadd(simd8<T>::load(&prob[i+0*HAP_NUMBER]), _nt, _tFreq) * _emit[ah0];
		simd8<T> _p1 = fmadd(simd8<T>::load(&prob[i+1*HAP_NUMBER]), _nt, _tFreq) * _emit[ah1];
		simd8<T> _p2 = fmadd(simd8<T>::load(&prob[i+2*HAP_NUMBER]), _nt, _tFreq) * _emit[ah2];
		simd8<T> _p3 = fmadd(simd8<T>::load(&prob[i+3*HAP_NUMBER]), _nt, _tFreq) * _emit[ah3];
		_sum0 = _sum0 + _p0; _sum1 = _sum1 + _p1; _sum2 = _sum2 + _p2; _sum3 = _sum3 + _p3;
		_p0.store(&prob[i+0*HAP_NUMBER]); _p1.store(&prob[i+1*HAP_NUMBER]);
		_p2.store(&prob[i+2*HAP_NUMBER]); _p3.store(&prob[i+3*HAP_NUMBER]);
	}
	for( ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = hget(hrow, k);
		simd8<T> _p = fmadd(simd8<T>::load(&prob[i]), _nt, _tFreq) * _emit[ah];
		_sum0 = _sum0 + _p;
		_p.store(&prob[i]);
	}
	simd8<T> _sum = (_sum0 + _sum1) + (_sum2 + _sum3);
	_sum.store(&probSumH[0]);
	probSumT = hsum8(probSumH);
}

template <typename T> inline
void haplotype_segment<T>::COLLAPSE_AMB() {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	simd8<T> _sum0 = simd8<T>::splat(0.0f), _sum1 = simd8<T>::splat(0.0f);
	simd8<T> _sum2 = simd8<T>::splat(0.0f), _sum3 = simd8<T>::splat(0.0f);
	simd8<T> _tFreq = simd8<T>::splat(yt / n_cond_haps);
	simd8<T> _nt = simd8<T>::splat(nt / probSumT);
	simd8<T> _emit[2] = { simd8<T>::loadu(&g0[0]), simd8<T>::loadu(&g1[0]) };
	const unsigned char * __restrict hrow = Hvar.row_ptr(curr_rel_locus+curr_rel_locus_offset);
	unsigned int n4 = n_cond_haps & ~3u, k = 0; int i = 0;
	for( ; k < n4 ; k += 4, i += 4*HAP_NUMBER) {
		bool ah0 = hget(hrow, k+0);
		bool ah1 = hget(hrow, k+1);
		bool ah2 = hget(hrow, k+2);
		bool ah3 = hget(hrow, k+3);
		simd8<T> _p0 = fmadd(simd8<T>::splat(probSumK[k+0]), _nt, _tFreq) * _emit[ah0];
		simd8<T> _p1 = fmadd(simd8<T>::splat(probSumK[k+1]), _nt, _tFreq) * _emit[ah1];
		simd8<T> _p2 = fmadd(simd8<T>::splat(probSumK[k+2]), _nt, _tFreq) * _emit[ah2];
		simd8<T> _p3 = fmadd(simd8<T>::splat(probSumK[k+3]), _nt, _tFreq) * _emit[ah3];
		_sum0 = _sum0 + _p0; _sum1 = _sum1 + _p1; _sum2 = _sum2 + _p2; _sum3 = _sum3 + _p3;
		_p0.store(&prob[i+0*HAP_NUMBER]); _p1.store(&prob[i+1*HAP_NUMBER]);
		_p2.store(&prob[i+2*HAP_NUMBER]); _p3.store(&prob[i+3*HAP_NUMBER]);
	}
	for( ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = hget(hrow, k);
		simd8<T> _p = fmadd(simd8<T>::splat(probSumK[k]), _nt, _tFreq) * _emit[ah];
		_sum0 = _sum0 + _p;
		_p.store(&prob[i]);
	}
	simd8<T> _sum = (_sum0 + _sum1) + (_sum2 + _sum3);
	_sum.store(&probSumH[0]);
	probSumT = hsum8(probSumH);
}

/*******************************************************************************/
/*****************			MISSING GENOTYPE			************************/
/*******************************************************************************/

template <typename T> inline
void haplotype_segment<T>::INIT_MIS() {
	fill(prob.begin(), prob.end(), 1.0f/(HAP_NUMBER * n_cond_haps));
	fill(probSumH.begin(), probSumH.end(), 1.0f/HAP_NUMBER);
	probSumT = 1.0f;
}

template <typename T> inline
void haplotype_segment<T>::RUN_MIS() {
	simd8<T> _sum0 = simd8<T>::splat(0.0f), _sum1 = simd8<T>::splat(0.0f);
	simd8<T> _sum2 = simd8<T>::splat(0.0f), _sum3 = simd8<T>::splat(0.0f);
	simd8<T> _tFreq = simd8<T>::load(&probSumH[0]) * simd8<T>::splat(yt / (n_cond_haps * probSumT));
	simd8<T> _nt = simd8<T>::splat(nt / probSumT);
	unsigned int n4 = n_cond_haps & ~3u, k = 0; int i = 0;
	for( ; k < n4 ; k += 4, i += 4*HAP_NUMBER) {
		simd8<T> _p0 = fmadd(simd8<T>::load(&prob[i+0*HAP_NUMBER]), _nt, _tFreq);
		simd8<T> _p1 = fmadd(simd8<T>::load(&prob[i+1*HAP_NUMBER]), _nt, _tFreq);
		simd8<T> _p2 = fmadd(simd8<T>::load(&prob[i+2*HAP_NUMBER]), _nt, _tFreq);
		simd8<T> _p3 = fmadd(simd8<T>::load(&prob[i+3*HAP_NUMBER]), _nt, _tFreq);
		_sum0 = _sum0 + _p0; _sum1 = _sum1 + _p1; _sum2 = _sum2 + _p2; _sum3 = _sum3 + _p3;
		_p0.store(&prob[i+0*HAP_NUMBER]); _p1.store(&prob[i+1*HAP_NUMBER]);
		_p2.store(&prob[i+2*HAP_NUMBER]); _p3.store(&prob[i+3*HAP_NUMBER]);
	}
	for( ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		simd8<T> _p = fmadd(simd8<T>::load(&prob[i]), _nt, _tFreq);
		_sum0 = _sum0 + _p;
		_p.store(&prob[i]);
	}
	simd8<T> _sum = (_sum0 + _sum1) + (_sum2 + _sum3);
	_sum.store(&probSumH[0]);
	probSumT = hsum8(probSumH);
}

template <typename T> inline
void haplotype_segment<T>::COLLAPSE_MIS() {
	simd8<T> _sum0 = simd8<T>::splat(0.0f), _sum1 = simd8<T>::splat(0.0f);
	simd8<T> _sum2 = simd8<T>::splat(0.0f), _sum3 = simd8<T>::splat(0.0f);
	simd8<T> _tFreq = simd8<T>::splat(yt / n_cond_haps);
	simd8<T> _nt = simd8<T>::splat(nt / probSumT);
	unsigned int n4 = n_cond_haps & ~3u, k = 0; int i = 0;
	for( ; k < n4 ; k += 4, i += 4*HAP_NUMBER) {
		simd8<T> _p0 = fmadd(simd8<T>::splat(probSumK[k+0]), _nt, _tFreq);
		simd8<T> _p1 = fmadd(simd8<T>::splat(probSumK[k+1]), _nt, _tFreq);
		simd8<T> _p2 = fmadd(simd8<T>::splat(probSumK[k+2]), _nt, _tFreq);
		simd8<T> _p3 = fmadd(simd8<T>::splat(probSumK[k+3]), _nt, _tFreq);
		_sum0 = _sum0 + _p0; _sum1 = _sum1 + _p1; _sum2 = _sum2 + _p2; _sum3 = _sum3 + _p3;
		_p0.store(&prob[i+0*HAP_NUMBER]); _p1.store(&prob[i+1*HAP_NUMBER]);
		_p2.store(&prob[i+2*HAP_NUMBER]); _p3.store(&prob[i+3*HAP_NUMBER]);
	}
	for( ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		simd8<T> _p = fmadd(simd8<T>::splat(probSumK[k]), _nt, _tFreq);
		_sum0 = _sum0 + _p;
		_p.store(&prob[i]);
	}
	simd8<T> _sum = (_sum0 + _sum1) + (_sum2 + _sum3);
	_sum.store(&probSumH[0]);
	probSumT = hsum8(probSumH);
}

/*******************************************************************************/
/*****************					SUM Ks				************************/
/*******************************************************************************/

template <typename T> inline
void haplotype_segment<T>::SUMK() {
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		probSumK[k] = prob[i+0] + prob[i+1] + prob[i+2] + prob[i+3] + prob[i+4] + prob[i+5] + prob[i+6] + prob[i+7];
	}
}

/*******************************************************************************/
/*****************		TRANSITION COMPUTATIONS			************************/
/*******************************************************************************/

template <typename T> inline
bool haplotype_segment<T>::TRANS_HAP() {
	sumHProbs = 0.0f;
	unsigned int  curr_rel_segment_index = curr_segment_index-segment_first;
	yt = M.getForwardTransProb(AlphaLocus[curr_rel_segment_index - 1], prev_abs_locus);
	nt = 1.0f - yt;
	T fact1 = nt / AlphaSumSum[curr_rel_segment_index - 1];
	for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
		T fact2 = (AlphaSum[curr_rel_segment_index-1][h1]/AlphaSumSum[curr_rel_segment_index-1]) * yt / n_cond_haps;
		simd8<T> _sum = simd8<T>::splat(0.0f);
		for (int k = 0 ; k < n_cond_haps ; k ++) {
			simd8<T> _alpha = simd8<T>::splat(Alpha[curr_rel_segment_index-1][k*HAP_NUMBER + h1] * fact1 + fact2);
			simd8<T> _beta = simd8<T>::load(&prob[k*HAP_NUMBER]);
			_sum = _sum + _alpha * _beta;
		}
		_sum.store(&HProbs[h1*HAP_NUMBER]);
		sumHProbs += HProbs[h1*HAP_NUMBER+0]+HProbs[h1*HAP_NUMBER+1]+HProbs[h1*HAP_NUMBER+2]+HProbs[h1*HAP_NUMBER+3]+HProbs[h1*HAP_NUMBER+4]+HProbs[h1*HAP_NUMBER+5]+HProbs[h1*HAP_NUMBER+6]+HProbs[h1*HAP_NUMBER+7];
	}
	return (std::isnan(sumHProbs) || std::isinf(sumHProbs) || sumHProbs < std::numeric_limits<T>::min());
}

template <typename T> inline
bool haplotype_segment<T>::TRANS_DIP_MULT() {
	sumDProbs= 0.0f;
	double scaling = 1.0 / sumHProbs;
	for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
		if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
			for (int nd = 0 ; nd < 64 ; ++nd) {
				if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
					DProbs[t] = (((double)HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)]) * scaling) * ((double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)]) * scaling);
					sumDProbs += DProbs[t];
					t++;
				}
			}
		}
	}
	return (std::isnan(sumDProbs) || std::isinf(sumDProbs) || sumDProbs < std::numeric_limits<double>::min());
}

template <typename T> inline
bool haplotype_segment<T>::TRANS_DIP_ADD() {
	sumDProbs = 0.0f;
	double scaling = 1.0 / sumHProbs;
	for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
		if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
			for (int nd = 0 ; nd < 64 ; ++nd) {
				if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
					DProbs[t] = DProbs[t] = (((double)HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)]) * scaling) + ((double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)]) * scaling);
					sumDProbs += DProbs[t];
					t++;
				}
			}
		}
	}
	return (std::isnan(sumDProbs) || std::isinf(sumDProbs) || sumDProbs < std::numeric_limits<double>::min());
}

template <typename T> inline
void haplotype_segment<T>::IMPUTE(std::vector < float > & missing_probabilities) {
	const unsigned char * __restrict hrow = Hvar.row_ptr(curr_rel_locus+curr_rel_locus_offset);
	simd8<T> _sumA[2] = { simd8<T>::splat(0.0f), simd8<T>::splat(0.0f) };
	simd8<T> _alphaSum = simd8<T>::splat(1.0f) / simd8<T>::load(&AlphaSumMissing[curr_rel_missing][0]);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = hget(hrow, k);
		simd8<T> _prob = simd8<T>::load(&prob[i]);
		simd8<T> _alpha = simd8<T>::load(&AlphaMissing[curr_rel_missing][i]);
		_sumA[ah] = _sumA[ah] + (_alpha * _alphaSum) * _prob;
	}
	T * prob0 = reinterpret_cast<T*>(&_sumA[0]);
	T * prob1 = reinterpret_cast<T*>(&_sumA[1]);
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = prob1[h] / (prob0[h]+prob1[h]);
	}
}

#endif
