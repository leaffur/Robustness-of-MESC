args = (commandArgs(TRUE))
chr = as.numeric(args[[1]])
simul = as.numeric(args[[2]])
h_g = as.numeric(args[[3]])
# N = as.numeric(args[[4]])
# Ne = as.numeric(args[[5]])
dir = args[[4]]
h_med_fix = as.numeric(args[[5]])
G = 100

.libPaths("/gpfs/gibbs/project/zhao/cl2384/R/4.2")

mesc <- "/gpfs/gibbs/pi/zhao/cl2384/CONDA/conda_envs/mesc/run_mesc.py"

# step 04: generate mesc / ldsc commands


h_med = 0.2


sink(paste0(dir, "mesc.sh"))
for (tt in 1:G){
  for (h_med in (c(0, 0.2, 0.4, 0.6, 0.8) * h_g)){
    step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ", 
                    dir, "exprun --ref-ld ", dir, "wholeLD --w-ld ", 
                    dir, "wholeLD --frqfile ", dir, "gtexgeno.frq --out ", dir, "temp/unh2med_hm", h_med, "_t", tt,
                    " --chisq-max 500")
    cat(step4, "\n")
    
    step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ", 
                    dir, "expr --ref-ld ", dir, "baselineLD --w-ld ", 
                    dir, "wholeLD --frqfile ", dir, "gtexgeno.frq --out ", 
                    dir, "temp/h2med_hm", h_med, "_t", tt, " --chisq-max 500")
    cat(step4, "\n")
    
    step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ", 
                    dir, "exprun_mesc.", chr, " --ref-ld ", dir, "wholeLD --w-ld ", 
                    dir, "wholeLD --frqfile ", dir, "gtexgeno.frq --out ", dir, "/temp/unh2med_hm", h_med, "_t", tt,
                    "_mesc --chisq-max 500")
    cat(step4, "\n")
    
    step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ", 
                    dir, "expr_mesc.", chr, " --ref-ld ", dir, "baselineLD --w-ld ", 
                    dir, "wholeLD --frqfile ", dir, "gtexgeno.frq --out ", 
                    dir, "temp/h2med_hm", h_med, "_t", tt, "_mesc --chisq-max 500")
    cat(step4, "\n")
  }
}
sink()

# dsq --job-file mesc.txt --partition day,pi_zhao,scavenge --mem-per-cpu 2g -t 0-01:00:00 --mail-type ALL --output /dev/null