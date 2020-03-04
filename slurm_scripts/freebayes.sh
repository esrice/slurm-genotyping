#!/bin/bash
#SBATCH --error logs/freebayes.%a.err
#SBATCH --output logs/freebayes.%a.out

module load freebayes/1.3

((chunk_number = SLURM_ARRAY_TASK_ID - 1))
chunk=$(printf "r%04d" $chunk_number)
bams_list=$(cut -f3 lists/alignments.tsv | tr '\n' ' ')

freebayes -t regions/$chunk -f $ref $bams_list | bgzip > regions/${chunk}.vcf.gz
