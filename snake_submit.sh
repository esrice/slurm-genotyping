#snakemake -j 30 --cluster-config cluster.yaml \
#    --cluster "sbatch -t {cluster.time} -p {cluster.partition} \
#    --mem {cluster.mem} -J {rule} --output {cluster.out} \
#    --error {cluster.err} --ntasks-per-node {threads}"

snakemake -j 30 --rerun-incomplete --cluster-config cluster.yaml \
    --cluster ./slurm_submit.py --cluster-status ./slurm_status.py \
    --restart-times 1
