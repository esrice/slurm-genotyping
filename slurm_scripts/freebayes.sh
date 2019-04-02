#!/bin/bash
#SBATCH --error logs/bwa_mem.%a.err
#SBATCH --output logs/bwa_mem.%a.out

region=$(head -n $SLURM_ARRAY_TASK_ID lists/regions.bed | tail -1)
bams_list=$(cut -f3 lists/alignments.tsv | tr '\n' ' ')

module load freebayes/1.2

freebayes --region $region -f $ref $bams_list | gzip > \
    regions/${SLURM_ARRAY_TASK_ID}.vcf.gz
