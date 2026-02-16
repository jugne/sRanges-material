#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=1608

beast -resume penguins_inf_morph_at_start_$SLURM_ARRAY_TASK_ID.xml
