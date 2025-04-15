#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alessia-mapelli@fht.org
#SBATCH --output=/group/diangelantonio/users/alessia_mapelli/Brain_simulations/Sim_1/logs/sequential_run_lambda_100_1404.log
#SBATCH --partition=cpuq
#SBATCH --mem=16GB
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --job-name=Brain_sim_1
#SBATCH --time=07:00:00

source /center/healthds/singularity_functions
cd $HOME
Rscript /group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Sim_1c_seq.R



