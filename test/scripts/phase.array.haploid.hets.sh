#!/bin/bash

# Declares diploid samples as haploid, so the input holds heterozygous genotypes in haploid samples.
# phase_common must set them missing and impute them as homozygous.

set -euo pipefail

INP=array/target.haploid.bcf
SUB=tmp/target.haploid.hets.samples.txt
HAP=tmp/target.haploid.hets.haploids.txt
OUT=tmp/target.haploid.hets.phased.bcf
mkdir -p tmp

bcftools query -l $INP | awk '/^Dip_/ && ++n <= 500' > $SUB
awk 'NR <= 50' $SUB > $HAP
bcftools view -S $SUB -Ob -o tmp/target.haploid.hets.bcf $INP
bcftools index -f tmp/target.haploid.hets.bcf

../phase_common/bin/phase_common --input tmp/target.haploid.hets.bcf --haploids $HAP --region 1 --map info/chr1.gmap.gz --output $OUT --thread 8

NIN=$(bcftools view -S $HAP tmp/target.haploid.hets.bcf | bcftools query -i 'GT="het"' -f '[%GT\n]' | wc -l)
NOUT=$(bcftools view -S $HAP $OUT | bcftools query -i 'GT="het"' -f '[%GT\n]' | wc -l)
echo "Hets in haploid samples: input=$NIN / output=$NOUT"
if [ "$NIN" -eq 0 ]; then echo "FAIL: the input holds no hets in haploid samples"; exit 1; fi
if [ "$NOUT" -ne 0 ]; then echo "FAIL: haploid samples keep hets after phasing"; exit 1; fi
echo "PASS"
