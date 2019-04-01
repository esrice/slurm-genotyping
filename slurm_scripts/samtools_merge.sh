#!/bin/bash
#SBATCH --error logs/samtools_merge.%a.err
#SBATCH --output logs/samtools_merge.%a.out

inputs=$(head -n $SLURM_ARRAY_TASK_ID lists/merges.tsv | cut -f1 | tail -1)
output=$(head -n $SLURM_ARRAY_TASK_ID lists/merges.tsv | cut -f2 | tail -1)

module load samtools/1.9

samtools merge -f $output $inputs
