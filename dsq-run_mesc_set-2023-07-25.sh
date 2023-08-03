#!/bin/bash
#SBATCH --output dsq-run_mesc_set-%A_%1a-%N.out
#SBATCH --array 0-6
#SBATCH --job-name dsq-run_mesc_set
#SBATCH --partition day,pi_zhao --mem-per-cpu 30g -t 0-23:59:00 --mail-type ALL

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/run_mesc_set.txt --status-dir /gpfs/gibbs/pi/zhao/cl2384/MESC_Genes

