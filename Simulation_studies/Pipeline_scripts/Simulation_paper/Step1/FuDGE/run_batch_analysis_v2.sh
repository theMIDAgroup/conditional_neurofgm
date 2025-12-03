#!/bin/bash

# Configuration
BASE_PATH="/group/diangelantonio/users/alessia_mapelli/Brain_simulations/systematic_simulations/Step1"
SCRIPT1="/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Pipeline_scripts/Simulation_paper/Step1/FuDGE/Litterature_comparison.R"
SCRIPT2="/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Pipeline_scripts/Simulation_paper/Step1/Check_results_screening_procedure_tests_v2.R"
SCRIPT3="/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Pipeline_scripts/Simulation_paper/Step1/FuDGE/Check_results_litt_comp.R"
JOBS_LIMIT=200

# Change to base directory
cd "$BASE_PATH"

# Main processing
echo "Starting batch analysis..."
echo "Base path: $BASE_PATH"
echo "Script 1: $SCRIPT1"
echo "Script 2: $SCRIPT2"
echo "Script 3: $SCRIPT3"
echo "========================================"

# Find all simulation directories (p*_ng1*_ng2* pattern)
for results_dir in "$BASE_PATH"/p*_n*_n*; do

    dir_name="$(basename "$results_dir")"
    echo "Processing directory: $dir_name"
    
    # Find all seed directories within this simulation
    for seed_dir in "$results_dir"/seed_*; do
        while [ "$(squeue -u $USER | wc -l)" -ge "$JOBS_LIMIT" ]; do
            echo "Jobs limit reached, sleeping for a while..."
            sleep 240
        done
        seed_name="$(basename "$seed_dir")"
        echo "  Processing seed: $seed_name"
        
        # Process iteration config files
        for config_file in "$seed_dir"/config_iter_*.yaml; do
            LOG_FILE="${seed_dir}/results_lit_comparison/logs/Lit_comparison.log"
            RES=$(sbatch --parsable --output="$LOG_FILE" "/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Pipeline_scripts/Simulation_paper/Step1/FuDGE/Sbatch_parallel_litt_comp.sbatch" "$SCRIPT1" "$config_file")
            echo "Running job ID: $RES"
        done
        
        echo "  Completed seed: $seed_name"
    done

    echo "Completed directory: $dir_name"
    echo ""
done

#while [ "$(squeue -u $USER | wc -l)" -ge 2 ]; do
#    echo "Jobs limit reached, sleeping for a while..."
#    sleep 240
#done


#for results_dir in "$BASE_PATH"/p*_n*_n*; do
#
#    dir_name="$(basename "$results_dir")"
#    echo "Processing directory: $dir_name"
#
#    seed_dir_vec=$"$results_dir"/seed_*
#    seed_dir=$seed_dir_vec[1]
#    config_file="$seed_dir"/config_iter_*.yaml
    
#    Rscript $SCRIPT2 $config_file
#    Rscript $SCRIPT3 $config_file
#
#    echo "Completed directory: $dir_name"
#    echo ""
# done
