# slurm-genotyping
Genotyping on SLURM. This package is made to be used on the SLURM
cluster 'crane' at the UNL HCC, but should be able to be used on other slurm
clusters with some modification.

## Installation
Navigate to where you want to do your work and download this repository:
```
git clone https://github.com/esrice/slurm-genotyping.git
```
Everything will be in the `slurm-genotyping` subdirectory after that. You can
rename this directory to something that describes your project better if you
want.

## Running
This pipeline does not include read QC, as this is something I like to do
manually. Once you've got your reads in good shape, put them in the `reads`
subdirectory. Symlinks are fine.

To tell the pipeline where your reference genome and reads are, edit
`config.yaml` and change the paths to the paths of the fasta and indices for
your reference genome. Also use the config file to tell this pipeline which
fastq files correspond to which individuals. There are detailed instructions in
the comments of this file, and example entries to start with.

### Mapping reads
This pipeline uses `bwa mem` to map reads to the reference. To submit mapping
jobs to the cluster, run
```
./genotype.py map
```
If it works, you will get a message from SLURM letting you know that a batch job
has been submitted, and a message from the pipeline telling you to wait until
all of the jobs are done running and then run the next step. After all the jobs
are done, you can check to make sure they worked by reading the logs, which are
in `logs/bwa_mem.[chunk_number].[err|out]`. If you get out of time/memory
errors, you'll need to increase the limits in `config.yaml`.

### Generating regions
To parallelize the task of variant-calling, this pipeline breaks the genome into
small regions. Some parts of the genome may have higher coverage than others, so
in order to even out resource usage among the jobs, we break the genome into
pieces of variable size based on their read coverage.

