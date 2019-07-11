#!/bin/bash
#conda activate ribopod
snakemake --snakefile Snakefile\
    --config config_path=configs/$1.py\
    --js $PWD/jobscript.sh\
    --printshellcmds\
    --cluster-config $PWD/cluster.yaml\
    --jobname "{rulename}.{jobid}.$1"\
    --keep-going\
    --stats $PWD/$1.stats\
    --rerun-incomplete\
    -j 200\
    --cluster 'sbatch --partition={cluster.partition} --ntasks={cluster.cores} --mem={cluster.mem} --time={cluster.time} -o {cluster.logout} -e {cluster.logerror}'
