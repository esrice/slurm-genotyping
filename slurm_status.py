#!/usr/bin/env python3

# this script is based on
# https://groups.google.com/forum/#!topic/snakemake/7cyqAIfgeq4

import sys
import subprocess

def parse_key_value(scontrol_output):
    """
    Given the output of scontrol, parse the key/value pairs
    into a dict and return it
    """
    params = {}
    for key_value_pair in scontrol_output.split():
        name, var = key_value_pair.partition("=")[::2]
        params[name.strip()] = var
    return params

def main():
    # this script gets run by snakemake to check the job status,
    # with the single argument of the slurm job ID.
    jobid = sys.argv[1]

    # run scontrol to get information about the given job
    out = subprocess.run(['scontrol', 'show', 'jobid', jobid],
            stdout=subprocess.PIPE).stdout.decode('utf-8')

    # get the job state from the scontrol output
    try:
        state = parse_key_value(out)['JobState']

        map_state = {"PENDING": 'running',
                     "RUNNING": 'running',
                     "SUSPENDED": 'running',
                     "CANCELLED": 'failed',
                     "COMPLETING": 'running',
                     "COMPLETED": 'success',
                     "CONFIGURING": 'running',
                     "FAILED": 'failed',
                     "TIMEOUT": 'failed',
                     "PREEMPTED": 'failed',
                     "NODE_FAIL": 'failed',
                     "REVOKED": 'failed',
                     "SPECIAL_EXIT": 'failed',
                     "": 'success'}

        # output whether the job succeeded, failed, or is still
        # running, based on the JobState string
        print(map_state[state])
    except: # if something didn't work, assume job failed
        print('failed')

if __name__ == '__main__':
    main()
