#!/usr/bin/env python3

import argparse
import yaml
import subprocess
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser(prog='genotype.py')

    parser.add_argument('step', choices=['map', 'merge', 'call', 'regions'],
            help='step to run (map, merge, or call)')

    return parser.parse_args()

def run_mapping(config):
    """ Run the mapping step of the pipeline. """
    # first, go through the yaml and make a list of
    # alignments that need to be run, formatting it in a
    # way that a SLURM script can easily read:
    # tab-separated format, one alignment job per line,
    # '{fastqs_in}\t{bam_out}', e.g.,
    # 'lib1_R1.fastq.gz,lib1_R1.fastq.gz\tlib1.bam'
    alignments_list = open('lists/alignments.tsv', 'w')
    num_alignments = 0
    for individual in config['individuals']:
        for i, lane in enumerate(individual['reads']):
            input_filenames = 'reads/{}'.format(lane['r1'])
            if 'r2' in lane:
                input_filenames += ',reads/{}'.format(lane['r2'])

            output_filename = 'alignments/{}_{}.bam'.format(
                    individual['name'], i)

            print('\t'.join([individual['name'], input_filenames,
                output_filename]), file=alignments_list)
            num_alignments += 1

    alignments_list.close()

    # next, run sbatch to get the bwa mem jobs going
    command = ['sbatch', '--mem', '{}G'.format(config['bwa_memory_GB']),
        '--cpus-per-task', str(config['bwa_threads']),
        '--time', '{}-00:00:00'.format(config['bwa_time_days']),
        '--array', '1-{}'.format(num_alignments),
        '--export', 'threads={},bwa_ref={}'.format(
            config['bwa_threads'], config['bwa_ref']),
        'slurm_scripts/bwa_mem.sh']
    subprocess.run(command, check=True, stdout=sys.stdout, stderr=sys.stderr)
    print('Submitted read mapping job array. Wait until all jobs are done\n'
            'and then run merge step.', file=sys.stderr)

def run_merge(config):
    """ Run the merge step of the pipeline. """

    # first, make a list of merges to perform that can be
    # easily read by a SLURM script
    merges_list = open('lists/merges.tsv', 'w')
    num_merges = 0
    for individual in config['individuals']:
        if len(individual['reads']) == 1: # no merging necessary
            # just rename the bam file
            os.rename('alignments/{}_0.bam'.format(individual['name']),
                'alignments/{}.bam'.format(individual['name']))
            pass
        else: # merging necessary
            input_filenames = 'alignments/{}_*.bam'.format(individual['name'])
            output_filename = 'alignments/{}.bam'.format(individual['name'])
            print('\t'.join([input_filenames, output_filename]),
                file=merges_list)
            num_merges += 1
    merges_list.close()

    # then, submit the merge jobs to SLURM as an array,
    # but only if there are any.
    if num_merges > 0:
        command = ['sbatch', '--mem', '{}G'.format(config['merge_memory_GB']),
            '--time', '{}-00:00:00'.format(config['merge_time_days']),
            '--array', '1-{}'.format(num_merges),
            'slurm_scripts/samtools_merge.sh']
        subprocess.run(command, check=True, stdout=sys.stdout,
            stderr=sys.stderr)
        print('Submitted merging job array. Wait until all jobs are done\n'
            'and then run regions step.', file=sys.stderr)
    else:
        print('No merges were necessary. You can run the regions step now',
            file=sys.stderr)

def run_regions(config):
    """ Run the regions step of the pipeline. """
    # rather than calculating total coverage across all
    # samples, which would take a long time, we just pick
    # an arbitrary sample (the first) and hope that it is
    # reasonably representative of the others
    infile = 'alignments/{}.bam'.format(config['individuals'][0]['name'])
    command = ['sbatch', '--mem', '{}G'.format(config['regions_memory_GB']),
        '--time', '{}-00:00:00'.format(config['regions_time_days']),
        '--export', 'infile={},faidx={},num_regions={}'.format(
            infile, config['faidx_ref'], config['freebayes_num_regions']),
        'slurm_scripts/make_regions.sh']
    #subprocess.run(command, check=True, stdout=sys.stdout, stderr=sys.stderr)
    print(' '.join(command), file=sys.stderr)
    print('Submitted region-generating script. Wait until complete and then\n'
        'run the next step.')

def main():
    args = parse_args()

    # TODO check that yaml has everything required
    config = yaml.load(open('config.yaml', 'r'))

    if args.step == 'map':
        run_mapping(config)
    elif args.step == 'merge':
        run_merge(config)
    elif args.step == 'regions':
        run_regions(config)
    elif args.step == 'call':
        run_freebayes(config)

if __name__ == '__main__':
    main()
