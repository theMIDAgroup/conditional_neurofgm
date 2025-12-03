#!/bin/bash
source /center/healthds/singularity_functions

# Accept config file path as the first argument
CONFIG_FILE="$1"
JOBS_LIMIT=70

# Extract parameters from config
NODES=$(grep '^p:' "$CONFIG_FILE" | awk '{print $2}')
SAVE_PATH=$(grep '^save_path:' "$CONFIG_FILE" | awk '{print $2}')
SIM_NAME=$(grep '^simulation_name:' "$CONFIG_FILE" | awk '{print $2}')
ITER=$(grep '^iteration:' "$CONFIG_FILE" | awk '{print $2}')
R_SCRIPT=$(grep '^r_script:' "$CONFIG_FILE" | awk '{print $2}')

LOG_DIR="${SAVE_PATH}${SIM_NAME}/seed_${ITER}/results/logs"
echo ${LOG_DIR}

# Submit jobs
for i in $(seq 1 "$NODES"); do
  while [ "$(squeue -u $USER | wc -l)" -ge "$JOBS_LIMIT" ]; do
    echo "Jobs limit reached, sleeping for a while..."
    sleep 240
  done

  echo "Processing: Node ${i}"
  INPUT=${i}
  LOG_FILE="${LOG_DIR}/Node_${i}.log"
  RES=$(sbatch --parsable --output="$LOG_FILE" "Sbatch_parallel_tests.sbatch" "$INPUT" "$R_SCRIPT" "$CONFIG_FILE")
  echo "Running job ID: $RES"
done



  

