chr = 1
simul = 1 # simul = 1:5
dir = "/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/chr1/s1_samplesize/"
dir2 <- "/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/s1_samplesize/"

ldsc <- "module load miniconda; conda activate ldsc; /gpfs/gibbs/pi/zhao/cl2384/CONDA/conda_envs/ldsc/ldsc.py"
mesc <- "module load miniconda; conda activate mesc; /gpfs/gibbs/pi/zhao/cl2384/CONDA/conda_envs/mesc/run_mesc.py"
gunzip <- "gunzip -f "
gzip <- "gzip -f "

gtex_id <- read.table(paste0(dir, "gtex_id.txt"), header = TRUE)
gtexexpr <- read.table(paste0(dir, "gtexexpr.txt"), header = TRUE)

c1 <- match(gtex_id[, 1], names(gtexexpr)[-(1:3)])


for (i in c(500, 1000, 1500, 2000, 2500)){
  system(paste0("mkdir -p ", dir, "exp_tmp_", i, ""))
  system(paste0("mkdir -p ", dir, "expun_tmp_", i, ""))
  
  write.table(gtex_id[1:i, ],
              file = paste0(dir, "gtex_id_", i, ".txt"),
              row.names = FALSE, quote = FALSE, col.names = TRUE)
  
  system(paste0("plink --bfile ", dir, "gtexgeno --chr ", chr, " --keep ", dir, "gtex_id_", i, ".txt --make-bed --out ", dir, "gtexgeno_", i))
  # freq
  system(paste0("plink -bfile ", dir, "gtexgeno_", i, " --freq --out ", dir, "gtexgeno_", i))
  
  write.table(gtexexpr[, c(1,2,3, 3 + c1[1:i])], file  = paste0(dir, "gtexexpr_", i, ".txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
}
# 
# # expscore_MESC
# sink(paste0(dir2, "exprscore_mesc_sub.txt"))
# for (i in c(500, 1000, 1500, 2000, 2500)){
#   sent <- paste0(mesc, " --compute-expscore-indiv --plink-path /vast/palmer/apps/avx2/software/PLINK/1.9b_6.21-x86_64/plink --expression-matrix ",
#          dir, "gtexexpr_", i, ".txt --exp-bfile ", dir, "gtexgeno_", i, " --geno-bfile ",
#          dir, "gtexgeno_", i, " --chr ", chr, " --out ", dir, "expr_mesc_", i, " --tmp ",
#          dir, "exp_tmp_", i, "")
#   cat(sent, "\n")
# 
#   sent <- paste0(mesc, " --compute-expscore-indiv --plink-path /vast/palmer/apps/avx2/software/PLINK/1.9b_6.21-x86_64/plink --expression-matrix ",
#          dir, "gtexexpr_", i, ".txt --exp-bfile ", dir, "gtexgeno_", i, " --geno-bfile ",
#          dir, "gtexgeno_", i, " --chr ", chr, " --out ", dir, "exprun_mesc_", i, " --tmp ",
#          dir, "expun_tmp_", i, "", " --num-bins 1")
#   cat(sent, "\n")
# }
# sink()