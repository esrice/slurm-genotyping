# slurm-genotyping
Genotyping with snakemake on SLURM. This package is made to be used on the SLURM
cluster 'crane' at the UNL HCC, but should be able to be used on other slurm
clusters with some modification.

## Installation
This pipeline uses the Snakemake workflow manager. First, if you don't have
miniconda, install it:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Now, install Snakemake:
```
conda install -c bioconda -c conda-forge snakemake
```

Finally, navigate to where you want to do your work and download this repository:
```
git clone https://github.com/esrice/slurm-genotyping.git
```
Everything will be in the `slurm-genotyping` subdirectory after that.

## Running
This pipeline does not include read QC, as this is something I like to do
manually. Once you've got your reads in good shape, put them in the `reads`
subdirectory with names `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`, where
`{sample}` is some descriptive identifier such as library name. Symlinks are
fine.

Next, edit `config.yaml` and change the paths to the paths of the fasta and
indices for your reference genome.

Finally, run the pipeline with the following command:
```
sbatch snake_submit.sh
```

This will break the task into small jobs and submit each of them to the cluster.
If something goes wrong, you can fix it and then run the above command again
and it will start back up right from where it left off. You can monitor the
pipeline's progress in `snakemake.out`.

If everything works, you'll have a vcf at `variant-calls/raw.vcf.gz` at the end.
