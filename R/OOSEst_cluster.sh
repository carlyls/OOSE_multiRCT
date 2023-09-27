#!/bin/bash

#SBATCH --array=1-7200                     # how many tasks in the array
#SBATCH --output=outputs/o.%a.out   # file to collect standard output
#SBATCH --error=errors/err.%a.log    # file to collect standard output

#SBATCH --cpus-per-task=3       # number of cores
#SBATCH --nodes=1               # number of nodes
#SBATCH --mem=50GB               # memory per __node__

module load conda_R

R CMD BATCH --no-save Simulations_OOSEst.R log/e.$SLURM_ARRAY_TASK_ID
exit 0
