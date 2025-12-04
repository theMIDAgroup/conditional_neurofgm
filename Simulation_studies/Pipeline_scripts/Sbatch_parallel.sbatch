#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alessia.mapelli@fht.org
#SBATCH --partition=cpuq
#SBATCH --mem=1GB
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name=Brain_sim
#SBATCH --time=00:20:00

source /center/healthds/singularity_functions
cd $HOME

input=$1
r_script=$2
conf_file=$3

Rscript "${r_script}" "${conf_file}" "${input}"
