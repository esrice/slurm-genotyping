#!/bin/bash
#SBATCH --error logs/samtools_merge.%a.err
#SBATCH --output logs/samtools_merge.%a.out

slurm_scripts/coverage_to_regions.py $infile > lists/regions.bed
