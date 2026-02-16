#!/bin/bash

#SBATCH --job-name=canid
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=1608

beast -statefile canids_srfbd_${SLURM_ARRAY_TASK_ID}.state  -resume canids_rho_${SLURM_ARRAY_TASK_ID}.xml
