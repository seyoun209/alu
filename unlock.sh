#! /bin/bash -login
  
## Exit if any command fails
set -e

## Load required modules
module load python/3.9.6

## Activate virtual environment
source env/bin/activate

case $1 in

    '-h' | '--help' | ' ' | '')
            echo -e '\e[31mSpecify which workflow to unlock (i.e. run_RNAprocessing)'
            exit 2
            ;;
        'rna' | 'run_rnaPipe')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/rna_process.smk --configfile "config/rna_process.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
        'signal' | 'run_signalPipe')
            ## Unlock snakemake workflow
            snakemake -s snakefiles/signal_process.smk --configfile "config/rna_process.yaml" --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
        'vcfpreprocess' | 'run_vcfpreprocess')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/run_vcfpreprocess --configfile "config/rna_prcoess.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
        'signal' | 'run_RNAsignal')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/RNA_signal.snakefile --configfile "config/rna_prcoess.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
        'sqtl' | 'run_sqtl')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/sqtl.snakefile --configfile "config/rna_prcoess.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;

esac

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun."                                                              
