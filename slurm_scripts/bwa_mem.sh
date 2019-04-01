#!/bin/bash
#SBATCH --error logs/bwa_mem.%a.err
#SBATCH --output logs/bwa_mem.%a.out

sample=$(head -n $SLURM_ARRAY_TASK_ID lists/libraries.tsv | cut -f1 | tail -1)
inputs=$(head -n $SLURM_ARRAY_TASK_ID lists/libraries.tsv | cut -f2 | tail -1)
output=$(head -n $SLURM_ARRAY_TASK_ID lists/libraries.tsv | cut -f3 | tail -1)

module load bwa/0.7 samtools/1.9

bwa mem -t $threads -R "@RG\tID:${sample}\tSM:${sample}" \
    $bwa_ref $(tr ',' ' ' < $inputs) | \
    samtools view -bh - | samtools sort - | samtools rmdup - $output
samtools index $output

