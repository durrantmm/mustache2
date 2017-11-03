#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=res.txt
#
#SBATCH --ntasks=1
#SBATCH --time=120-00:00:
#SBATCH --mem-per-cpu=100
#SBATCH --export=ALL
#$ -cwd

snakemake -p -j 999 --cluster-config slurm_cluster.json --cluster "sbatch --export=ALL -A {cluster.account} -n {cluster.n} -t {cluster.time} -w 'sgiuv300-srcf-d10-01'"
