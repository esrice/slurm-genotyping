#!/bin/bash
#SBATCH --error logs/make_regions.err
#SBATCH --output logs/make_regions.out

slurm_scripts/coverage_to_regions.py $infile > lists/regions.bed

num_lines=$(wc -l < lists/regions.bed)
((lines_per_file = (num_lines + num_files - 1) / num_files))
split -d -a4 -l $lines_per_file lists/regions.bed regions/r
