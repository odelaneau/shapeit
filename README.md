# Segmented HAPlotype Estimation and Imputation Tools version 5 (SHAPEIT5)

> **Note:** This repository replaces the original [odelaneau/shapeit5](https://github.com/odelaneau/shapeit5) repo, which has been repeatedly disabled by GitHub. This is now the official home for the SHAPEIT5 source code.

SHAPEIT5 estimates haplotypes in large datasets, with a special focus on rare variants.

- **docker**: all script needed to build a docker file comprising all binaries
- **docs**: documentation in html
- **ligate**: ligate multiple phased BCF/VCF files into a chromosome length file
- **resources**: genetics maps in b37 and b38 amd coordinates for WGS chunks
- **phase_common**: phase common sites, typically SNP array data
- **phase_rare**: phase rare variants onto a scaffold of common variants
- **static_bins**: static binaries of all executables
- **simulate**: simulate simple haplotype datasets
- **switch**: compute switch error rate and genotyping error rate given simulated or trio data
- **tasks**: scripts used to phase large datasets, good base to start pipelining
- **test**: simulated data for first-step testing of the method
- **versions**: versioning
- **xcftools**: tools to handle XCF files [experimental]

## Download

It is strongly advised to prioritize the use of the latest released version (see released section on the right hand side) instead of cloning the latest version of the source code, as doing so is at your own risk.

## Apple Silicon (macOS)

SHAPEIT5 builds and runs natively on Apple Silicon (arm64) since v5.1.2: the AVX2/FMA SIMD kernels of the HMM engine were ported to ARM NEON, and the build system automatically uses clang + Homebrew dependencies on macOS.

```bash
brew install htslib boost zlib bzip2 lzma curl openssl
make macos    # native arm64 build of all binaries
```

### Why native (and why now)

The previous releases ship **x86_64 static binaries** (and a Docker image of the same). On Apple Silicon those only run under the **Rosetta 2** translation layer. Apple has announced that **macOS 27 is the last release to support Rosetta 2** for general Intel apps; starting with **macOS 28 (expected 2027)**, Intel-only binaries will no longer launch. The native arm64 build removes that dependency entirely, and is also faster.

Measured on an M1 Pro (10 cores), `phase_common` on the haploid test dataset (`--thread 8`), native NEON vs the original x86_64 AVX2 binary run under Rosetta 2:

| Build | Wall time (avg of 2) |
|---|---|
| **Native NEON (arm64)** | **~76 s** |
| x86_64 AVX2 under Rosetta 2 | ~188 s |

So the native build is **~2.5x faster** than running the x86_64 binary through Rosetta 2, and is not subject to the Rosetta 2 end of life.

Notes:

- `make` also works on macOS (Darwin is detected automatically); `make macos` is the explicit target.
- For maximum performance, run with `--thread N` where N is the number of physical cores (10-12 on M3).
- The Docker image in `docker/` contains prebuilt **x86_64** static binaries: on a Mac it runs via Rosetta 2 (~2.5x slower, and only until macOS 28). A multi-arch arm64 Linux image is a planned follow-up.
- The `static_exe` targets (Linux x86_64 static binaries) do not apply on macOS.

## Citation
If you use SHAPEIT5 in your research work, please cite the following paper:

Hofmeister RJ, Ribeiro DM, Rubinacci S., Delaneau O. [Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank](https://www.nature.com/articles/s41588-023-01415-w). Nature Genetics (2023)

## Documentation

Documentation, installation instructions and tutorials can be found at:

[https://odelaneau.github.io/shapeit/](https://odelaneau.github.io/shapeit/)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
