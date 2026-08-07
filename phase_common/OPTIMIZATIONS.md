# phase_common — Optimization Opportunities

Candidate optimizations for consolidating code, reducing runtime, and cutting memory
usage. Grouped by theme, each item notes the **location**, the **problem**, the
**proposed change**, and a rough **impact / effort / risk** estimate.

The three hot paths are:
1. **HMM forward/backward** — `models/haplotype_segment.{h,cpp}`
2. **PBWT selection & phasing sweep** — `containers/conditioning_set/*.cpp`
3. **Conditioning-state assembly** — `objects/compute_job.cpp`

---

## Status — implemented so far (2026-08-07)

- **A1 — DONE.** `haplotype_segment_single`/`_double` merged into one
  `template <typename T> class haplotype_segment` in `models/haplotype_segment.{h,cpp}`.
  The four old files were deleted. The type-dependent SIMD kernels are written once
  against a thin `simd8<T>` wrapper (see the follow-up refactor below).
- **B2(a) — DONE.** The six hot per-locus kernels (`RUN_HOM`, `COLLAPSE_HOM`,
  `RUN_AMB`, `COLLAPSE_AMB`, `RUN_MIS`, `COLLAPSE_MIS`) unroll the `k`-loop by 4
  into independent accumulators. **B2(b)** (AVX-512 / 2-haps-per-`__m512`) is **not**
  done.
- **B8 — DONE (bit-exact).** The homozygous mismatch `if (ag!=ah)` is now branchless
  via `_emit[2] = {splat(1.0f), splat(mismatch)}` indexed by `ag!=ah`.
- **B1(a) — DONE (bit-exact).** Added `bitmatrix::row_ptr()`; each kernel hoists a
  `__restrict` row pointer out of the `k`-loop and reads bits via `hget(hrow, k)`,
  so the invariant row base is computed once instead of re-derived every iteration.
- **B4 (Alpha) — DONE (bit-exact).** `forward()` replaces the per-segment
  `Alpha[seg] = prob` deep copy with an O(1) `prob.swap(Alpha[seg])`; the
  `AlphaMissing` copy is reordered before the swap. (The `AlphaMissing` copy itself
  is untouched — that's C2.)
- **`simd8<T>` refactor — DONE.** After A1/B2/B8 the header had grown to ~780 lines
  of duplicated `if constexpr (float)/else` kernels. Introduced a header-local
  `simd8<T>` value type ("8 lanes of `T`": `__m256` for float, a `__m256d` pair for
  double) plus a `hsum8()` reduction helper, and rewrote every kernel single-path.
  Result: **~497 lines**, one algorithm per kernel, float codegen/perf unchanged.

**Numerics note:** B2's unroll reorders the accumulation, so the **float** path is no
longer bit-identical to the pre-B2 build (Tier 2 — validate switch-error). B8 on its
own is bit-exact. The `simd8<T>` unification also gave the **double** (underflow-
fallback) path the same unroll, so it too shifted slightly; and it corrected an
apparent `AlphaSumMissing[1]`→`[4]` indexing typo in the old double `IMPUTE`.

> Line/file references in the **not-yet-done** items below predate these changes:
> `haplotype_segment_single.{h,cpp}` / `_double.{h,cpp}` are now the single merged
> `haplotype_segment.{h,cpp}`, and line numbers have shifted.

---

## A. Code consolidation

### A1. Merge `haplotype_segment_single` and `haplotype_segment_double` into one template — ✅ DONE
- **As implemented:** one `template <typename T> class haplotype_segment` in
  `models/haplotype_segment.{h,cpp}`; the four old files were removed. Type-dependent
  SIMD kernels are written once against the `simd8<T>` wrapper (not `if constexpr` —
  that was the intermediate step, since unified away). `forward`/`backward`/`SET_*`
  and all cursor logic are shared and `<float>`/`<double>` are explicitly instantiated.
- **Original problem/proposal below (for reference).**
- **Where:** `models/haplotype_segment_single.{h,cpp}` (116 + 231 lines) and
  `models/haplotype_segment_double.{h,cpp}` (116 + 230 lines).
- **Problem:** The two classes are ~95% identical. Every routine (`INIT_HOM`,
  `RUN_HOM`, `COLLAPSE_*`, `RUN_AMB`, `TRANS_*`, `IMPUTE`, `forward`, `backward`,
  ctor/dtor) exists twice, differing only in the scalar type (`float`/`double`) and
  the SIMD intrinsic width (`__m256` handling 8 floats vs. two `__m256d` handling
  4+4 doubles). This is the single largest maintenance liability in the module — a
  bug fix or algorithm change must be applied twice and kept in sync by hand
  (e.g. the double `DProbs[t] = DProbs[t] = ...` self-assign typo at
  `haplotype_segment_double.h:432` and `haplotype_segment_single.h:401`).
- **Proposal:** Convert to a single `template <typename T> class haplotype_segment`.
  Isolate the ~6 genuinely type-dependent SIMD kernels behind small `if constexpr`
  branches or a thin `simd_traits<T>` helper (broadcast, fmadd, horizontal-sum,
  store). All coordinate/cursor logic, `forward()`, `backward()`, `SET_*_TRANS`,
  and the diplotype loops become shared. Instantiate `<float>` and `<double>`.
- **Impact:** ~450 lines removed, one source of truth. No runtime change.
- **Effort:** Medium. **Risk:** Medium (numerics must be bit-verified against the
  current build on a regression set).

### A2. Factor the duplicated PBWT prefix-sort out of `select()` and `solve()`
- **Where:** `conditioning_set_selection.cpp:84-102` and
  `conditioning_set_solve.cpp:159-179`.
- **Problem:** The core PBWT update (the `u/v/p/q` split of arrays `A`/`C` into
  `A`/`B` + `C`/`D`, followed by the `std::copy` merge) is copy-pasted. `solve()`
  additionally maintains the rank array `R`.
- **Proposal:** Extract a single `pbwt_step(l, A, B, C, D)` inline helper (and an
  optional `R` update). Both callers invoke it.
- **Impact:** ~20 duplicated lines removed, guarantees the two sweeps stay
  algorithmically identical. **Effort:** Low. **Risk:** Low.

### A3. Remove dead / commented-out code
- **Where:** large commented `getMatchHetCount*` / `getMatchHets` block
  (`bitmatrix.cpp:79-189`), debug `std::cout`/`exit(1)` blocks
  (`bitmatrix.cpp:64-76`, `conditioning_set_solve.cpp:129`).
  *(The commented double-precision `RUN_AMB` variant that used to live in
  `haplotype_segment_single.h` is already gone — deleted by the single/double merge.)*
- **Proposal:** Delete. Git history preserves it.
- **Impact:** Readability only. **Effort:** Trivial. **Risk:** None.

---

## B. Runtime / speed

### B1. Kill the aliasing-blocked reloads in `Hvar.get()` in the innermost HMM loop — ✅ (a) DONE
- **As implemented (variant a, bit-exact):** added `bitmatrix::row_ptr(row)` and, in
  each of the 7 kernels that read `Hvar`, hoisted `const unsigned char * __restrict
  hrow = Hvar.row_ptr(curr_rel_locus+curr_rel_locus_offset);` out of the `k`-loop.
  Per-`k` access is now `hget(hrow, k)` = `(hrow[k>>3] >> (7-(k&7))) & 1`, identical
  bits to the old `get()`. The `__restrict` promises the `float` store doesn't touch
  the matrix, so the invariant `row*(n_cols>>3)` is computed once and the shared byte
  survives across the unrolled `k`'s. Variant (b) (byte expansion) not done.
- **Original problem/analysis below (for reference).**
- **Where:** every kernel in `haplotype_segment_single.h` /
  `_double.h` calls `Hvar.get(curr_rel_locus+curr_rel_locus_offset, k)` inside the
  `k`-loop (e.g. lines 128, 152, 174, 199, 285, 419).
- **What the compiler already handles (verified in -O3 asm):** the address
  arithmetic is *not* a bottleneck — `col/8` and `n_cols/8` compile to `shr`, `col%8`
  to `and $7`, and the bit test to a single `btl`. There is **no div/mod** here; the
  original rationale for this item was wrong.
- **Actual problem (verified):** because `bitmatrix::bytes` is `unsigned char*` and
  `unsigned char` may alias any type, the `_mm256_store_ps(&prob[i], …)` in the same
  loop forces the compiler to assume the `float` store could modify the matrix.
  Strict-aliasing LICM/CSE is therefore blocked, and the general-case loop reloads
  `Hvar.n_cols`, reloads the `bytes` pointer, **re-multiplies** the loop-invariant
  `row*(n_cols>>3)`, and re-fetches the *same* data byte for all 8 consecutive `k`
  that share it — every iteration:
  ```
  movq   (%rdi), %rax        ; reload n_cols
  shrq   $3, %rax
  imulq  %r11, %rax          ; recompute row*(n_cols>>3)  (invariant, not hoisted)
  addq   8(%rdi), %rax       ; reload bytes
  movzbl (%rax,%rcx), %ecx   ; reload the shared byte
  ```
- **Proposal:** Hoist `bytes` and `row*(n_cols>>3)` into loop-local variables (or
  mark the pointers `__restrict`) so the invariant base survives across iterations;
  and/or load each packed byte once and expand its 8 bits, building the SIMD
  emission vector via a blend of `_emit[0]`/`_emit[1]` from a byte-derived mask
  instead of 8 separate `get()` calls.
- **Impact:** Low–Medium — these are L1-cache hits plus one `imul` per iteration,
  running alongside the dependent `fmadd→mul→add` chain that is the true limiter
  (see B2). Real but not the headline win. **Effort:** Low (locals/`__restrict`) to
  Medium (byte expansion). **Risk:** Low–Medium.

### B2. Widen SIMD for the single-precision path (AVX-512, or process 2 haps/iter) — ✅ (a) DONE / (b) not done
- **As implemented:** proposal (a) only — the six hot per-locus kernels unroll the
  `k`-loop by 4 into independent accumulators `_sum0.._sum3`, combined as
  `(s0+s1)+(s2+s3)`. Runs on the existing `-mavx2 -mfma` build. **Not bit-exact**
  (reordered accumulation). Proposal (b) — AVX-512 2-haps-per-`__m512` — is **not**
  implemented; with the `simd8<T>` structure the *double* path could adopt `__m512d`
  as a near drop-in (wrapper-only), but the *float* 2-haps path needs a hap-pair
  abstraction and edits to every hot kernel (see D2).
- **Original problem/proposal below (for reference).**
- **Where:** single-precision kernels use one `__m256` for the 8 `HAP_NUMBER`
  lanes; the `k`-loop does one conditioning hap per iteration.
- **Problem:** With `HAP_NUMBER == 8`, a single `__m256` already fills exactly one
  hap's 8 probabilities, so there is no lane-level parallelism *across* haps. The
  loop is latency-bound on the dependent `fmadd`→`add` chain and the scalar
  `Hvar.get()`.
- **Proposal:** (a) Unroll the `k`-loop by 2–4 and accumulate into independent
  `_sum` registers to hide FMA latency (breaks the single-accumulator dependency
  chain). (b) On AVX-512 hardware, pack two haps into one `__m512` and halve the
  iteration count. Gate with `#ifdef __AVX512F__`.
- **Impact:** Medium–High on modern CPUs. **Effort:** Medium.
  **Risk:** Medium (numerics; needs `-mavx512f` build variant).
- **In plain words:** The loop works on 8 probabilities at a time, which exactly
  fill one SIMD register — so there is no spare room to work on a second
  haplotype at once, and each pass has to wait for the previous pass's
  `fmadd`→`add` to finish before it can start. The CPU sits idle in that queue.
  Two ways to keep it busy: (a) unroll the loop and keep 2–4 independent running
  totals so several `fmadd`s are in flight at once (works on any current CPU);
  (b) on AVX-512 chips, use a double-width register that holds 16 values, pack
  two haplotypes into it, and halve the number of passes.

### B3. Avoid the transcendental `expm1f` in the transition hot path
- **Where:** `hmm_parameters.cpp:58-66` (`getForwardTransProb`) and `68-74`.
  Called from `forward()`/`backward()` at every locus
  (`haplotype_segment_single.cpp:102,159`).
- **Problem:** The fast path (adjacent loci) returns the precomputed `t[prev_idx]`,
  but when loci are non-adjacent — which happens whenever rare alleles are skipped
  via `M.rare_allele` — it recomputes `expm1f(-0.04 * Neff * dist_cm / Nhap)`. A
  transcendental in the per-locus loop is expensive.
- **Proposal:** Precompute a cumulative genetic-distance / cumulative `-log(1-t)`
  array so any interval probability is a subtraction + one `exp`/`expm1`, or cache
  results keyed by `(prev_idx, curr_idx)` when the skip pattern is stable across the
  many MCMC iterations. At minimum, hoist the invariant `-0.04*Neff/Nhap` factor.
- **Impact:** Medium (scales with rare-allele density and iteration count).
  **Effort:** Low–Medium. **Risk:** Low.

### B4. Eliminate full-vector copies when storing Alpha / AlphaMissing — ✅ (Alpha) DONE
- **As implemented (bit-exact):** the segment-boundary `Alpha[seg] = prob` deep copy
  is replaced by an O(1) `prob.swap(Alpha[seg])` in `forward()`. This is safe because
  the forward recursion advances across segments through the reduced `probSumK`
  (`SUMK`→`COLLAPSE`), never re-reading a stored `Alpha`; the swapped-in stale `prob`
  is fully overwritten by the next segment's `COLLAPSE`. The per-missing-site
  `AlphaMissing[m] = prob` copy is now done **before** the swap (else it would read the
  swapped-out buffer). The `AlphaMissing` copy itself is **not** eliminated — that's
  the C2 reduction. Values stored are identical → bit-exact.
- **Original problem/proposal below (for reference).**
- **Where:** `haplotype_segment_single.cpp:122-128`:
  `Alpha[...] = prob;`, `AlphaSum[...] = probSumH;`, and
  `AlphaMissing[curr_rel_missing] = prob;`.
- **Problem:** `Alpha[seg] = prob` deep-copies `HAP_NUMBER * n_cond_haps` floats
  **once per segment**, and `AlphaMissing[m] = prob` copies the same size **once per
  missing site**. For large K and many missing sites this is a major fraction of
  `forward()` time and drives the allocation in B-memory below.
- **Proposal:** Compute forward `prob` directly into the `Alpha[seg]` slot
  (`std::swap` or write-in-place) instead of copy-after-the-fact. For
  `AlphaMissing`, store only what `IMPUTE` actually consumes (see M3) rather than
  the whole prob vector.
- **Impact:** Medium–High. **Effort:** Medium. **Risk:** Medium.

### B5. Reuse per-thread scratch in PBWT `select()` / `solve()`
- **Where:** `conditioning_set_selection.cpp:70-73` and
  `conditioning_set_solve.cpp:57-65` allocate `A,B,C,D` (and `R,G,Het,Mis,Amb`)
  of size `n_hap`/`n_ind` **on every `select(chunk)` / `solve(chunk)` call**.
- **Problem:** Repeated heap alloc/free of large vectors, once per chunk per
  iteration.
- **Proposal:** Hoist these into per-worker scratch buffers (indexed by
  `id_worker`, like `threadData` in the phaser) and `iota`/`fill`-reset them.
- **Impact:** Low–Medium (alloc churn + first-touch page faults).
  **Effort:** Low. **Risk:** Low.

### B6. Replace the coarse global mutex around stats with atomics
- **Where:** `phaser_algorithm.cpp:48-51,85-87` — `statH.push`, `statS.push`,
  and `Kbanned.pushIBD2` are serialized under `mutex_workers` inside the per-window
  loop for every job.
- **Problem:** With many threads this mutex is taken twice per window and can
  serialize otherwise-parallel HMM work; `statH/statS` pushes are trivial and could
  use thread-local accumulators merged at the end.
- **Proposal:** Give each worker thread-local `stats`/`Kbanned` buffers, merge once
  after `pthread_join`. Keep the mutex only for the genuinely shared job counter.
- **Impact:** Medium at high thread counts. **Effort:** Medium. **Risk:** Medium
  (must preserve `pushIBD2` semantics).

### B7. Skip the redundant duplicated conditional in `forward()`
- **Where:** `haplotype_segment_single.cpp:120-121`:
  ```
  if (curr_segment_locus == (G->Lengths[curr_segment_index] - 1)) SUMK();
  if (curr_segment_locus == G->Lengths[curr_segment_index] - 1) { ... }
  ```
- **Problem:** The identical condition is evaluated twice back-to-back.
- **Proposal:** Merge into a single `if` block.
- **Impact:** **None at runtime** — verified that GCC 11.4 `-O3` already CSEs the load
  and comparison (nothing writes `Lengths` between them). Keep only for readability.
  **Effort:** Trivial. **Risk:** None.

### B8. Make the homozygous mismatch branchless (bit-exact) — ✅ DONE
- **As implemented:** `RUN_HOM`/`COLLAPSE_HOM` use `simd8<T> _emit[2] = {splat(1.0f),
  splat(M.ed/M.ee)}` and `_p = _p * _emit[ag!=ah]` — no data-dependent branch, and
  multiplying by exactly `1.0f` is IEEE-exact, so output is byte-identical. (`INIT_HOM`
  builds `prob` directly and has no such multiply, so it needed no change.)
- **Original problem/proposal below (for reference).**
- **Where:** `haplotype_segment_single.h` — `INIT_HOM` (`_mm256_set1_ps((ag==ah)?…)`)
  and `RUN_HOM`/`COLLAPSE_HOM` (`if (ag!=ah) _prob = _mm256_mul_ps(_prob,_mismatch)`).
- **Problem (verified in `-O3` asm, GCC 11.4):** these compile to a *data-dependent
  branch* — `INIT_HOM` splits into a match tail and a `.L327` mismatch tail
  (`cmpb %al,%sil ; jne`), and `COLLAPSE_HOM` guards its `vmulps` with `je .L224`.
  `ah` (allele carriers) is effectively random, so the branch mispredicts ~50%
  (~15-cycle penalty). The AMB kernels are already branchless in source and GCC
  compiles them to a stack `_emit[]` table — i.e. the branchless form is exactly what
  the compiler prefers, it just can't legally introduce it here from a ternary/`if`.
- **Proposal:** Use the AMB pattern: `__m256 _emit[2]={set1(1.0f), set1(mismatch)};`
  then `_prob = _mm256_mul_ps(_prob, _emit[ag!=ah]);` (and `_prob = _emit[ag!=ah]` in
  `INIT_HOM`). Multiplying by exactly `1.0f` is IEEE-exact → **byte-identical output**.
- **Impact:** Medium (removes a ~50%-mispredicting branch from the per-locus loop).
  **Effort:** Low. **Risk:** Low (bit-exact).

### B9. Cache `M.ed/M.ee` as a member (bit-exact)
- **Where:** every kernel recomputes `M.ed/M.ee`. **Verified:** `INIT_HOM` emits a
  `vdivsd` **inside** the mismatch path *every iteration*; `RUN_HOM`/`COLLAPSE_HOM`
  hoist it but still divide once per kernel call (millions of calls per run).
- **Proposal:** Add `float mismatch;` set once in the ctor (`= M.ed/M.ee;`) and reuse
  it everywhere (incl. the `g0`/`g1` setup). Exact same value → bit-exact.
- **Impact:** Low–Medium (kills a per-iteration divide in `INIT_HOM`, a per-call
  divide elsewhere). **Effort:** Trivial. **Risk:** None.

### B10. Derive `probSumT` from the SIMD register, not a scalar reload
- **Where:** every kernel does `_mm256_store_ps(&probSumH[0], _sum)` then immediately
  re-reads the 8 floats scalar-wise: `probSumT = probSumH[0]+…+probSumH[7]`.
- **Problem:** The scalar reload right after the vector store incurs a
  store-to-load-forwarding stall, once per kernel call.
- **Proposal:** Keep the store (`probSumH` is consumed later) but compute `probSumT`
  via an in-register horizontal reduction of `_sum` (`_mm256_hadd_ps`×2 +
  `extractf128` + add). **Reorders the 8-way sum → not bit-exact** (Tier 2).
- **Impact:** Low–Medium (per-call stall × call count). **Effort:** Low.
  **Risk:** Low–Medium (validate switch-error).

### B11. Restructure `TRANS_HAP` to load `prob[k*8]` once
- **Where:** `haplotype_segment.h::TRANS_HAP` — loops `h1` (outer, 8) over `k`
  (inner, K), reloading `_beta = prob[k*HAP_NUMBER]` for every `h1`.
- **Problem:** The whole `prob` array is read **8×** (once per `h1`). Called once per
  segment boundary in `backward()`.
- **Proposal:** Interchange to `k` outer: load `_beta` once, then FMA into 8
  accumulators `acc[h1]` with a broadcast of `Alpha[…][k*8+h1]*fact1+fact2`. Reads
  `prob` K times instead of 8K. **Reorders sums → not bit-exact** (Tier 2).
- **Impact:** Medium (memory traffic on the transition path). **Effort:** Medium.
  **Risk:** Medium (validate switch-error).
- **Prior attempt:** implemented once (k-outer + `fmadd`) then reverted — the
  benchmark showed no measurable gain, but that measurement is suspect (possible
  benchmark error). `TRANS_HAP` is per-segment-boundary (not per-locus) and `prob`
  is L1-resident, so the win may genuinely be small; re-measure before re-committing.

> **Compiler check (GCC 11.4, `-O3 -mavx2 -mfma`):** all of B1, B2, B8–B11 were
> confirmed **not** already performed by the compiler (the per-iteration `imulq` +
> reloads for B1; the branch for B8; the `vdivsd` for B9; single-accumulator loops
> with no unroll for B2; no reassociation/interchange for B10/B11 — GCC cannot do
> these without `-ffast-math`, which the module avoids for underflow detection). Only
> B7 is already handled by the compiler.

---

## C. Memory usage

### C1. Store Alpha with a compact type / avoid double the footprint
- **Where:** `haplotype_segment_*.h:75-80` — `Alpha`, `AlphaSum`, `AlphaMissing`,
  `AlphaSumMissing`.
- **Problem:** `Alpha` is `n_segments × HAP_NUMBER × n_cond_haps` scalars **held for
  the whole window**. In the double path this doubles to 8 bytes/element. This is
  the dominant per-thread working-set and it is allocated fresh for every job.
- **Proposal:** (a) Default to the single-precision (float) path everywhere and only
  fall back to double per-segment — the double object currently mirrors the entire
  Alpha store in 8-byte elements even though underflow is rare (see the fallback at
  `phaser_algorithm.cpp:67-72`). (b) Investigate storing `AlphaSum`-normalized alphas
  so only the ratio needed by `TRANS_HAP` is retained. (c) Reserve/pool these
  vectors per worker instead of reallocating per job.
- **Impact:** High (peak RSS scales with K × window length × threads).
  **Effort:** Medium–High. **Risk:** Medium.

### C2. `AlphaMissing` stores a full prob vector per missing site
- **Where:** `haplotype_segment_*.h:79`, filled at
  `haplotype_segment_single.cpp:128`.
- **Problem:** `AlphaMissing` = `n_missing × HAP_NUMBER × n_cond_haps` floats. For
  samples with many missing genotypes this rivals `Alpha` in size, yet `IMPUTE`
  (`:412-430`) only needs the per-hap forward mass combined with the stored
  backward alpha.
- **Proposal:** Store only the reduced quantity `IMPUTE` consumes, or compute
  imputation inline during the backward pass without materializing all missing
  alphas. **Impact:** Medium–High (only when missingness is high).
  **Effort:** Medium. **Risk:** Medium.

### C3. Shrink `indexes_pbwt_neighbour` element type
- **Where:** `conditioning_set_header.h:50` — `std::vector<int>`
  `indexes_pbwt_neighbour` of size `depth × sites_pbwt_ngroups × 2·n_ind`.
- **Problem:** 4 bytes per entry to store a haplotype index. For large cohorts this
  is one of the biggest persistent allocations. If `n_hap < 2^32` it fits in `int`,
  but the `-1` sentinel and index range may allow `int32` → fine, though the
  transpose double-buffers it (`transposePBWTneighbours` uses a second half of the
  same array, `:53,64`), effectively 2× during transpose.
- **Proposal:** (a) Confirm the double-buffer is necessary; an in-place blocked
  transpose or a single reused scratch band would halve peak. (b) Consider `uint32`
  with an explicit sentinel to make the width intentional. **Impact:** Medium.
  **Effort:** Medium. **Risk:** Medium.

### C4. `DProbs` is a 4096-element array always resident
- **Where:** `haplotype_segment_*.h:82` —
  `DProbs[HAP_NUMBER^4]` = 4096 doubles = **32 KB** per segment object, even for the
  single-precision class.
- **Problem:** Only `curr_dipcount × prev_dipcount` (≤ 64×64 worst case, usually far
  fewer) entries are ever used (`TRANS_DIP_MULT`, `:379-389`). The full static
  array bloats the object and hurts cache/stack locality.
- **Proposal:** Size to the actual max transitions for the sample, or use the
  already-known `n_max_transitions` bound from `compute_job` (`compute_job.cpp:31`).
  **Impact:** Low–Medium (cache footprint of the hot object). **Effort:** Low.
  **Risk:** Low.

---

## D. Build / compiler flags

### D1. Enable more aggressive optimization flags
- **Where:** `common/makefile_common.mk:10` — `CXXFLAG=-O3 -mavx2 -mfma`;
  `static_exe` uses `-O2` (`:133`).
- **Proposals:**
  - Add `-funroll-loops` (helps the tight `k`-loops) and `-fstrict-aliasing`.
  - Add **LTO** (`-flto`) so the heavily-`inline`d HMM kernels and cross-TU calls
    (e.g. `bitmatrix::get`, `getForwardTransProb`) actually inline across objects.
  - Provide a `-march=native` / `-mtune=native` build variant for on-prem runs, and
    a separate `-mavx512f -mavx512bw` variant to unlock B2.
  - Bump `static_exe` from `-O2` to `-O3` (the shipped static binary is the one most
    users run — it is currently a tier below the local `desktop`/`rgc` builds).
- **Caveat:** Do **not** add `-ffast-math` blindly — the code relies on
  `isnan`/`isinf`/denormal checks for underflow recovery
  (`TRANS_HAP`/`TRANS_DIP_*`); fast-math would break that logic. If desired, apply
  `-fno-math-errno -fno-trapping-math` selectively instead.
- **Impact:** Medium (LTO + unroll on hot loops), low effort.
  **Risk:** Low–Medium (validate numerics per variant).

### D2. Verify `HAP_NUMBER == 8` alignment assumptions with AVX-512
- **Note:** aligned to 32 bytes (`aligned_vector32`, `__attribute__((aligned(32)))`).
  Moving to `__m512` (B2b) requires 64-byte alignment; bump the allocator alignment
  in a guarded build.
- **Post-`simd8<T>` update:** the double path is now nearly free to widen — add an
  `#ifdef __AVX512F__` specialization of `simd8<double>` backed by one `__m512d`
  (replacing the `__m256d` pair); every kernel, `hsum8`, and the `IMPUTE`
  `reinterpret_cast<double*>` keep working unchanged. Still needs the 64-byte
  alignment bump and a `-mavx512f -mavx512bw` build target. The float hot path gains
  nothing from a wider `simd8<float>` (8 floats already fill a `__m256`); its only
  AVX-512 win is 2-haps-per-`__m512`, which does not fit the one-hap `simd8` model.

---

## Suggested sequencing (highest value first)

| # | Item | Payoff | Speed-up | Effort | Risk |
|---|------|--------|----------|--------|------|
| ✅ | B8 — branchless HOM (bit-exact) | Medium speed | Med | Low | Low |
| ✅ | A1 — template-merge single/double HMM (+ `simd8<T>` unification) | Big consolidation | — | Med | Med |
| ✅ | B2(a) — SIMD unroll to break the FMA dep chain | High speed | High | Med | Med |
| 1 | B9 — cached `M.ed/M.ee` mismatch member (bit-exact) | Low-Med speed | Low-Med | Low | Low |
| ✅ | B1(a) — hoist `__restrict` row ptr out of the k-loop (bit-exact) | Low-Med speed | Low-Med | Low-Med | Low-Med |
| ✅ | B4 — Alpha stored via O(1) swap, not deep copy (bit-exact) | High memory + speed | Med | Med | Med |
| 3 | C1 — Alpha compact type / pooled per-worker buffers | High memory | — | Med-High | Med |
| 4 | B2(b)/D2 — AVX-512 (`__m512d` double drop-in; float hap-pair) | Med-High speed | Med-High | Low (double) / Med (float) | Med |
| 5 | B11 — restructure `TRANS_HAP` (load prob once) | Medium speed | Med | Med | Med |
| 6 | D1 — LTO + `-funroll-loops`, `-O3` static | Free-ish speed | Low-Med | Low | Low |
| 7 | C2 — reduce AlphaMissing storage | Memory (high-miss data) | — | Med | Med |
| 8 | B3 — cache transition probs | Speed | Med | Low | Low |
| 9 | B10 — register-side `probSumT` reduction | Low-Med speed | Low-Med | Low | Low-Med |
| 10 | B5/B6 — PBWT scratch reuse, stat atomics | Speed at scale | Med (at scale) | Low-Med | Low-Med |
| 11 | C3/C4 — shrink index & DProbs arrays | Memory | Low | Low-Med | Low |
| 12 | A2/A3/B7 — dedup + dead-code cleanup | Maintainability | — | Low | Low |

> **Bit-exact (no regression check needed):** B8, B9, B1, B4. **Numeric-affecting** (must
> be validated against a switch-error regression + underflow-recovery path): A1,
> B2, B3, B10, B11, C1–C4, D1 — the module deliberately trades precision for speed and
> relies on exact `isnan`/`isinf`/denormal underflow detection.
