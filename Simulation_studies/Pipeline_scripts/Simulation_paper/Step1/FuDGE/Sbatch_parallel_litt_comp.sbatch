#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alessia.mapelli@fht.org
#SBATCH --partition=cpuq
#SBATCH --mem=3GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=Brain_sim_litt_comp
#SBATCH --time=10:00:00

source /center/healthds/singularity_functions
cd $HOME

r_script=$1
conf_file=$2

Rscript "${r_script}" "${conf_file}" 
