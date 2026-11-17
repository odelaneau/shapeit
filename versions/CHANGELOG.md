---
title: CHANGELOG
layout: default
---

# CHANGELOG

All notable changes to this project are documented in this file.

{: .highlight }
Version 5.1.2 is now available [here] (https://github.com/odelaneau/shapeit5/releases) as docker image or static binaries. 

## v5.1.2
	* Apple Silicon (arm64) support: HMM AVX2/FMA kernels ported to ARM NEON, native macOS build via clang + Homebrew (`make macos`)
	* NEON kernels give comparable performance to AVX2 on Linux and ~10x the speed of the scalar code with compiler auto-vectorization on ARM
	* Native arm64 build is ~2.5x faster than the x86_64 binaries run under Rosetta 2, and removes the Rosetta 2 dependency (Rosetta 2 support for Intel apps ends with macOS 27 / macOS 28)

## v5.1.1
	* Bugfix of a nasty segfault, thanks to Ami Williams (23andMe). 
	* Bugfixes when reading multiple input files data	

## v5.1.0
	* Second release / Version used for phasing UKB WGS 200k data
	* Chromosome X support added in phase_[common/rare]
	* Pedigree phasing support added in phase_[common/rare/ligate]
	* Many small optimizations in phase_[common/rare]
	* IBD2 protection improved in phase_common abd added in phase_rare
	* Automatic indexing of VCF/BCF files	

## v5.0.0
	* First release / Version used for the preprint. (https://doi.org/10.1101/2022.10.19.512867)
	