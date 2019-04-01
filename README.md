# slurm-genotyping
Genotyping on SLURM. This package is made to be used on the SLURM
cluster 'crane' at the UNL HCC, but should be able to be used on other slurm
clusters with some modification.

## Installation
Navigate to where you want to do your work and download this repository:
```
git clone https://github.com/esrice/slurm-genotyping.git
```
Everything will be in the `slurm-genotyping` subdirectory after that.

## Running
This pipeline does not include read QC, as this is something I like to do
manually. Once you've got your reads in good shape, put them in the `reads`
subdirectory. Symlinks are fine.

To tell the pipeline where your reference genome and reads are, edit
`config.yaml` and change the paths to the paths of the fasta and indices for
your reference genome. Also use the config file to tell this pipeline which
fastq files correspond to which individuals. There are detailed instructions in
the comments of this file, and example entries to start with.

Finally, run the pipeline with the following command:
```
sbatch snake_submit.sh
```

This will break the task into small jobs and submit each of them to the cluster.
If something goes wrong, you can fix it and then run the above command again
and it will start back up right from where it left off. You can monitor the
pipeline's progress in `snakemake.out`.

If everything works, you'll have a vcf at `variant-calls/raw.vcf.gz` at the end.
