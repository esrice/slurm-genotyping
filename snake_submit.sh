#!/bin/bash
#SBATCH -p abg
#SBATCH --mem 16G
#SBATCH -t 7-00:00:00
#SBATCH -o snakemake.out
#SBATCH -e snakemake.err

# generate regions, if they have not yet been generated
if [ ! -d "variant-calls/regions" ]; then
    ./generate_regions.py
fi

# run snakemake
snakemake -j 30 --rerun-incomplete --cluster-config cluster.yaml \
    --cluster ./slurm_submit.py --cluster-status ./slurm_status.py \
    --restart-times 1
