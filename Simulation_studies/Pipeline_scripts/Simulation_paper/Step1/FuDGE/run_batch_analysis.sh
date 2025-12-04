#!/bin/bash

# Configuration
BASE_PATH="/group/diangelantonio/users/alessia_mapelli/Brain_simulations/systematic_simulations/Step1"
SCRIPT1="/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Pipeline_scripts/Simulation_paper/Step1/FuDGE/Litterature_comparison.R"
SCRIPT2="/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Pipeline_scripts/Simulation_paper/Step1/Check_results_screening_procedure_tests_v2.R"
SCRIPT3="/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Pipeline_scripts/Simulation_paper/Step1/FuDGE/Check_results_litt_comp.R"

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
        seed_name="$(basename "$seed_dir")"
        echo "  Processing seed: $seed_name"
        
        # Find all config files in this seed directory
        config_found=false
        
        # Process iteration config files
        for config_file in "$seed_dir"/config_iter_*.yaml; do
            Rscript $SCRIPT1 $config_file
        done
        
        echo "  Completed seed: $seed_name"
    done

    Rscript $SCRIPT2 $config_file
    Rscript $SCRIPT3 $config_file

    echo "Completed directory: $dir_name"
    echo ""
done

