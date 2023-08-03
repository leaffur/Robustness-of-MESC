#!/bin/bash
#SBATCH --output /dev/null
#SBATCH --job-name dsq-exprscore_mesc
#SBATCH --partition day,pi_zhao,scavenge --mem-per-cpu 10g -t 0-00:30:00 --mail-type ALL

module load PLINK/1.9b_6.21-x86_64; module load R/4.2.0-foss-2020b; Rscript step03_subsample.R
