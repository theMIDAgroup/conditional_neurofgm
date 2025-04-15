#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alessia-mapelli@fht.org
#SBATCH --output=/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1/batch_Bayes_parall.log
#SBATCH --partition=cpuq
#SBATCH --mem=16GB
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name=Brain_sim_2
#SBATCH --time=20:00:00

source /center/healthds/singularity_functions
cd $HOME
Rscript /group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Sim_1_batch.R



