#!/bin/bash

JOBS_LIMIT=70
r_script=/group/diangelantonio/users/alessia_mapelli/conditional_neurofgm/Simulation_studies/Sim_1_sbatch_parallel.R


for i in $(seq 1 64); do
  while [ "$(squeue -u $USER |wc -l)" -ge "${JOBS_LIMIT}" ]; do
		echo "Jobs limit reached, I sleep for a while";
		sleep 240
	done
  echo "Processing: Node ${i}"
  input=${i}
  RES=$(sbatch --parsable "Sbatch_parallel.sbatch" "${input}" "${r_script}" )
  echo "running job id: ${RES}"
done