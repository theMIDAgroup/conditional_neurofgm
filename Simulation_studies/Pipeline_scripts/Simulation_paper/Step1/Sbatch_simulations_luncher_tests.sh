#!/bin/bash

# Sbatch_simulations_launcher_tests.sh - Systematic simulation launcher
source /center/healthds/singularity_functions

# Change to the working directory
cd /group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Pipeline_scripts/Simulation_paper/Step1/

# Configuration
TEMPLATE_CONFIG="/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Pipeline_scripts/Simulation_paper/Step1/config_template.yaml"
YAML_MODIFIER="/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Pipeline_scripts/Simulation_paper/Step1/modify_yaml_config.py"
SAVE_PATH="/group/diangelantonio/users/alessia_mapelli/Brain_simulations/systematic_simulations/Step1/"
JOBS_LIMIT=100
tot_iteration=10

# Simulation parameter arrays (proper bash syntax)
p_seq=(10 15 25 50)
n_g1_seq=(50 75 100 150 200)
n_g2_seq=(50 75 100 150 200)  # Assuming same values as n_g1

for p in "${p_seq[@]}"; do
  for n_g1 in "${n_g1_seq[@]}"; do
    #for n_g2 in "${n_g2_seq[@]}"; do
      n_g2=$n_g1
      red_number=$((p / 3))
      simulation_name="p${p}_n${n_g1}_n${n_g2}"

      echo "=========================================="
      echo "Starting simulation: $simulation_name"
      echo "p=$p, n_g1=$n_g1, n_g2=$n_g2, red_number=$red_number"
      echo "=========================================="

      for iteration in $(seq 1 "$tot_iteration"); do
        echo "Processing iteration $iteration/$tot_iteration for $simulation_name"
        config_dir=$SAVE_PATH$simulation_name/seed_$iteration
        mkdir -p "$config_dir"
        current_config="$config_dir/config_iter_${iteration}.yaml"
        echo "Step 1: Creating config file: $current_config"
        python3 "$YAML_MODIFIER" "$TEMPLATE_CONFIG" "$current_config" \
                    "p=$p" \
                    "save_path=$SAVE_PATH" \
                    "n_g1=$n_g1" \
                    "n_g2=$n_g2" \
                    "red_number=$red_number" \
                    "simulation_name=$simulation_name" \
                    "iteration=$iteration" \
                    "tot_iteration=$tot_iteration"
        
        echo "  Step 2: Running data simulator (iteration $iteration)"

        Rscript Data_simulator_tests.R $current_config

        while [ "$(squeue -u $USER | wc -l)" -ge "$JOBS_LIMIT" ]; do
          echo "Jobs limit reached, sleeping for a while..."
          sleep 60
        done

        echo "Step 2: Launching parallel screening jobs (iteration $iteration)"
        bash Sbatch_parallel_luncher_tests.sh $current_config
      done

      echo "All iterations submitted for $simulation_name"
      echo "Waiting for all screening jobs to complete before running results checker..."
      
      while [ "$(squeue -u $USER | wc -l)" -ge 2 ]; do
        echo "Jobs limit reached, sleeping for a while..."
        sleep 240
      done
      Rscript Check_results_screening_procedure_tests.R $current_config

    #done
  done
done


