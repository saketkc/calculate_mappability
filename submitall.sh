#!/bin/bash
snakemake --snakefile Snakefile\
    --config config=$1\
    --js $PWD/jobscript.sh\
    --printshellcmds\
    --cluster-config $PWD/cluster.yaml\
    --jobname "$1.{jobid}.{rulename}"\
    --keep-going\
    --stats $PWD/$1.stats\
    --timestamp\
    --rerun-incomplete\
    --resources mem_mb=720000\
    -j 100\
    --cluster 'qsub -q cmb -S /bin/bash -V -l walltime={cluster.time} -l mem={cluster.mem} -l vmem={cluster.mem} -l pmem={cluster.mem} -l nodes=1:ppn={cluster.cores} -o {cluster.logdir} -e {cluster.logdir}'
