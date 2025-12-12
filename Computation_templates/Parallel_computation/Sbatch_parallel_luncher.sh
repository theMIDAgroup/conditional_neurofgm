#!/bin/bash

CONFIG_FILE="$1"

CURRENT_FOLDER=$(pwd)

R_SCRIPT="Script_sbatch_parallel.R"
JOBS_LIMIT=70

# Extract parameters from config
NODES=$(grep '^n_nodes:' "$CONFIG_FILE" | awk '{print $2}')
SAVE_PATH=$(grep '^output_path:' "$CONFIG_FILE" | awk '{print $2}')

LOG_DIR="${SAVE_PATH}/logs"
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
  RES=$(sbatch --parsable --output="$LOG_FILE" "Sbatch_parallel.sbatch" "$INPUT" "$R_SCRIPT" "$CONFIG_FILE" "$CURRENT_FOLDER")
  echo "Running job ID: $RES"
done



  

