# args = (commandArgs(TRUE))
# chr = as.numeric(args[[1]])
# simul = as.numeric(args[[2]])
# # h_g = as.numeric(args[[3]])
# # N = as.numeric(args[[4]])
# # Ne = as.numeric(args[[5]])
# dir = args[[3]]
# h_med_fix = as.numeric(args[[4]])
G = 100

.libPaths("/gpfs/gibbs/project/zhao/cl2384/R/4.2")

chr = 1
# h_g = 0.5


dir = "/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/chr1/s1_samplesize/"
ldsc <- "module load miniconda; conda activate ldsc; /gpfs/gibbs/pi/zhao/cl2384/CONDA/conda_envs/ldsc/ldsc.py"
mesc <- "module load miniconda; conda activate mesc; /gpfs/gibbs/pi/zhao/cl2384/CONDA/conda_envs/mesc/run_mesc.py"
dir2 <- "/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/s1_samplesize/"

# dir2 <- "/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/s1_samplesize/"


# dir <- paste0("/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/")
# dir1 <- paste0("/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/chr", chr, "/genotype/")
# dir2 <- paste0(dir, "chr", chr, "/")


# step 04: generate mesc / ldsc commands


h_med = 0.2
simul = 1




sink(paste0(dir2, "mesc.txt"))
for (tt in 1:G){
  # step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ", 
  #                 dir, "exprun --ref-ld ", dir, "wholeLD --w-ld ", 
  #                 dir, "wholeLD --frqfile ", dir, "gtexgeno.frq --out ", dir, "temp/unh2med_hm", h_med, "_t", tt,
  #                 " --chisq-max 500")
  # cat(step4, "\n")
  # 
  # step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ", 
  #                 dir, "expr --ref-ld ", dir, "baselineLD --w-ld ", 
  #                 dir, "wholeLD --frqfile ", dir, "gtexgeno.frq --out ", 
  #                 dir, "temp/h2med_hm", h_med, "_t", tt, " --chisq-max 500")
  # cat(step4, "\n")
  
  for (i in c(2500)){


    step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ",
                    dir, "exprun_mesc_", i, ".", chr, " --ref-ld ", dir, "wholeLD --w-ld ",
                    dir, "wholeLD --frqfile ", dir, "gtexgeno_", i, ".frq --out ", dir, "/temp/unh2med_hm", h_med, "_t", tt,
                    "_mesc_", i, " --chisq-max 500")
    cat(step4, "\n")

    step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ",
                    dir, "expr_mesc_", i, ".", chr, " --ref-ld ", dir, "baselineLD --w-ld ",
                    dir, "wholeLD --frqfile ", dir, "gtexgeno_", i, ".frq --out ",
                    dir, "temp/h2med_hm", h_med, "_t", tt, "_mesc_", i, " --chisq-max 500")
    cat(step4, "\n")
  }
}
sink()

# dsq --job-file mesc.txt --partition day,pi_zhao,scavenge --mem-per-cpu 2g -t 0-01:00:00 --mail-type ALL --output /dev/null


# step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ", 
#                 dir, "exprun --ref-ld ", dir, "wholeLD --w-ld ", 
#                 dir, "wholeLD --frqfile ", dir, "gtexgeno.frq --out ", dir, "temp/unh2med_hm", h_med, "_t", tt,
#                 " --chisq-max 500")
# write.table(step4, file = paste0(dir2, "mesc.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# 
# step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ", 
#                 dir, "expr --ref-ld ", dir, "baselineLD --w-ld ", 
#                 dir, "wholeLD --frqfile ", dir, "gtexgeno.frq --out ", 
#                 dir, "temp/h2med_hm", h_med, "_t", tt, " --chisq-max 500")
# write.table(step4, file = paste0(dir2, "smesc.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# 
# 
# step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ", 
#                 dir, "exprun_mesc.", chr, " --ref-ld ", dir, "wholeLD --w-ld ", 
#                 dir, "wholeLD --frqfile ", dir, "gtexgeno.frq --out ", dir, "/temp/unh2med_hm", h_med, "_t", tt,
#                 "_mesc --chisq-max 500")
# write.table(step4, file = paste0(dir2, "mesc1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
# 
# step4 <- paste0(mesc, " --h2med ", dir, "temp/sumstat_hm", h_med, "_t", tt, ".sumstats.gz --exp ", 
#                 dir, "expr_mesc.", chr, " --ref-ld ", dir, "baselineLD --w-ld ", 
#                 dir, "wholeLD --frqfile ", dir, "gtexgeno.frq --out ", 
#                 dir, "temp/h2med_hm", h_med, "_t", tt, "_mesc --chisq-max 500")
# write.table(step4, file = paste0(dir2, "smesc1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)


