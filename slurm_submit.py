#!/usr/bin/env python3

# this script is based on
# https://groups.google.com/forum/#!topic/snakemake/7cyqAIfgeq4

import os
import sys
from snakemake.utils import read_job_properties

job_script = sys.argv[1]
job_properties = read_job_properties(job_script)

cluster_params = {}

job_resources = job_properties["resources"]

cluster_params['hours'] = job_resources.get('hours', 24)
cluster_params['partition'] = job_properties['cluster']['partition']
cluster_params['mem'] = int(job_resources.get("mem", 12))
cluster_params['rule'] = job_properties['rule']
cluster_params['out'] = job_properties['cluster']['out']
cluster_params['err'] = job_properties['cluster']['err']
cluster_params['threads'] = job_properties.get('threads', 1)

os.system("sbatch -t {hours:02d}:00:00 -p {partition} --mem {mem}G -J {rule} \
        --output {out} --error {err} --ntasks-per-node {threads} {script}"
        .format(script=job_script, **cluster_params))
