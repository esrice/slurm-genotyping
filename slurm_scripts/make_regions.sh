#!/bin/bash
#SBATCH --error logs/samtools_merge.%a.err
#SBATCH --output logs/samtools_merge.%a.out

module load freebayes/1.2 bamtools/2.4

bamtools coverage < $infile | coverage_to_regions.py $faidx $num_regions \
    > regions.bed
