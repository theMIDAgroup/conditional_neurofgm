#!/bin/bash

JOBS_LIMIT=70
r_script=Script_sbatch_parallel.R

current_folder=$(pwd)
LOG_DIR="${current_folder}/logs"


for i in $(seq 1 64); do
  while [ "$(squeue -u $USER |wc -l)" -ge "${JOBS_LIMIT}" ]; do
		echo "Jobs limit reached, I sleep for a while";
		sleep 240
	done
  echo "Processing: Node ${i}"
  input=${i}
  LOG_FILE="${LOG_DIR}/Node_${i}.log"
  RES=$(sbatch --parsable --output="$LOG_FILE" "Sbatch_parallel.sbatch" "${input}" "${r_script}" )
  echo "running job id: ${RES}"
done