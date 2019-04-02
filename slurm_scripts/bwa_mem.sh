#!/bin/bash
#SBATCH --error logs/bwa_mem.%a.err
#SBATCH --output logs/bwa_mem.%a.out

sample=$(head -n $SLURM_ARRAY_TASK_ID lists/alignments.tsv | cut -f1 | tail -1)
inputs=$(head -n $SLURM_ARRAY_TASK_ID lists/alignments.tsv | cut -f2 | tail -1)
inputs=$(echo $inputs | tr ',' ' ')
output=$(head -n $SLURM_ARRAY_TASK_ID lists/alignments.tsv | cut -f3 | tail -1)

module load bwa/0.7 samtools/1.9

echo "threads: $threads; sample: $sample"
echo "bwa_ref: $bwa_ref"
echo "inputs: $inputs"
echo "output: $output"
bwa mem -t $threads -R "@RG\tID:${sample}\tSM:${sample}" $bwa_ref $inputs | \
    samtools view -bh - | samtools sort - | samtools rmdup - $output
samtools index $output

