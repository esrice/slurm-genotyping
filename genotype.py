#!/usr/bin/env python3

import argparse
import yaml
import subprocess
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser(prog='genotype.py')

    parser.add_argument('step', choices=['map', 'regions', 'call', 'join'],
            help='step to run (map, merge, or call)')

    return parser.parse_args()

def check_files_exist(filenames):
    """ Check if all files in a list exist and exit with an error if not. """
    for f in filenames:
        if not os.path.isfile(f):
            print('FATAL: Missing file: {}'.format(f), file=sys.stderr)
            sys.exit(1)

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
            input_filenames = ['reads/{}'.format(lane['r1'])]
            if 'r2' in lane:
                input_filenames.append('reads/{}'.format(lane['r2']))
            check_files_exist(input_filenames)

            output_filename = 'alignments/{}_{}.bam'.format(
                    individual['name'], i)

            print('\t'.join(map(str, [
                  individual['name'],
                  ','.join(input_filenames),
                  output_filename,
                ])), file=alignments_list)
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
            'and then run regions step.', file=sys.stderr)

def run_regions(config):
    """ Run the regions step of the pipeline. """
    # rather than calculating total coverage across all
    # samples, which would take a long time, we just pick
    # an arbitrary sample (the first) and hope that it is
    # reasonably representative of the others
    infile = 'alignments/{}_0.bam'.format(config['individuals'][0]['name'])
    command = ['sbatch', '--mem', '{}G'.format(config['regions_memory_GB']),
        '--time', '{}-00:00:00'.format(config['regions_time_days']),
        '--export', 'infile={},num_files={}'.format(infile,
            config['freebayes_num_jobs']),
        'slurm_scripts/make_regions.sh']
    subprocess.run(command, check=True, stdout=sys.stdout, stderr=sys.stderr)
    print('Submitted region-generating script. Wait until complete and then\n'
        'run the call step.')

def run_freebayes(config):
    """ Run the variant-calling step of the pipeline. """
    command = ['sbatch', '--mem', '{}G'.format(config['freebayes_memory_GB']),
        '--time', '{}-00:00:00'.format(config['freebayes_time_days']),
        '--array', '1-{}'.format(config['freebayes_num_jobs']),
        '--export', 'ref={}'.format(config['faidx_ref']),
        'slurm_scripts/freebayes.sh']
    subprocess.run(command, check=True, stdout=sys.stdout, stderr=sys.stderr)
    print('Submitted variant-calling job array. Wait until all jobs are done\n'
            'and then run join step.', file=sys.stderr)

def run_join(config):
    """ Run the join step of the pipeline. """
    command = ['sbatch', '--mem', '{}G'.format(config['join_memory_GB']),
        '--time', '{}-00:00:00'.format(config['join_time_days']),
        '--export', 'sortmem={}G'.format(config['join_memory_GB'] - 4),
        'slurm_scripts/join.sh']
    subprocess.run(command, check=True, stdout=sys.stdout, stderr=sys.stderr)
    print('Submitted join job. When it is done your final variant calls will\n'
            'be in variants.Q20.vcf.gz.', file=sys.stderr)

def main():
    args = parse_args()

    # TODO check that yaml has everything required
    config = yaml.load(open('config.yaml', 'r'))

    # TODO before each step, check that required files are
    # present so that the user gets an error from this
    # script and no jobs are submitted, rather than
    # submitting a bunch of jobs that will all fail and
    # making the user look through log files to figure out why.
    if args.step == 'map':
        run_mapping(config)
    elif args.step == 'regions':
        run_regions(config)
    elif args.step == 'call':
        run_freebayes(config)
    elif args.step == 'join':
        run_join(config)

if __name__ == '__main__':
    main()
