#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=1608

beast -statefile canids_fbd_fixed_$SLURM_ARRAY_TASK_ID.state  canids_fbd_fixed.xml
