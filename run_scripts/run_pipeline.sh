#!/bin/bash
#SBATCH -J sm
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 9000
#SBATCH --mem=10000

source activate snakemake
snakemake --snakefile workflow/Snakefile --profile ./profiles/slurm
