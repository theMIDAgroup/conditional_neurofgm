#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alessia.mapelli@fht.org
#SBATCH --output="/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Model A/A050S_2g_new_scores_3/results_of_screening_procedure/logs/75_unbalanced/%j_Node.log"
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

Rscript "${r_script}" "${input}"
