#!/bin/bash
#SBATCH --error logs/join.err
#SBATCH --output logs/join.out

module load freebayes/1.2 bcftools/1.9

zcat regions/r*.vcf.gz | vcffirstheader | bcftools filter -i 'QUAL>=20' -Ou \
    | bcftools sort -Ov -m $sortmem | vcfuniq | bgzip > variants.Q20.vcf.gz
tabix variants.Q20.vcf.gz
