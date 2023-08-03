args = (commandArgs(TRUE))
chr = as.numeric(args[[1]])
# simul = as.numeric(args[[2]])
# h_g = as.numeric(args[[3]])
# N = as.numeric(args[[4]])
# Ne = as.numeric(args[[5]])
dir = args[[2]]

.libPaths("/gpfs/gibbs/project/zhao/cl2384/R/4.2")


# chr = 1
# simul = 1 # simul = 1:5
# dir = "/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/chr1/s1_samplesize/"


# step 03
# module load miniconda; conda activate ldsc; 
ldsc <- "/gpfs/gibbs/pi/zhao/cl2384/CONDA/conda_envs/ldsc/ldsc.py"
mesc <- "/gpfs/gibbs/pi/zhao/cl2384/CONDA/conda_envs/mesc/run_mesc.py"
dir1 <- paste0("/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/chr", chr, "/genotype/")
dir2 = dir
# dir2 <- "/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/s1_samplesize/"
gunzip <- "gunzip -f "
gzip <- "gzip -f "

# LD scores

LDSCORE <- c(paste0(ldsc, " --bfile ", dir, "gtexgeno --l2 --annot ", dir1, "baselineLD.annot.gz --ld-wind-cm 1 --out ", dir, "baselineLD"),
             paste0(ldsc, " --bfile ", dir, "gtexgeno --l2 --ld-wind-cm 1 --out ", dir, "wholeLD"))
write.table(LDSCORE, file = paste0(dir2, "ldscore.sh"), quote = FALSE, row.names = FALSE, col.names = FALSE)


# expscore_MESC

system(paste0("mkdir -p ", dir, "exp_tmp"))
system(paste0("mkdir -p ", dir, "expun_tmp"))

EXPSCORE <- c(paste0(mesc, " --compute-expscore-indiv --plink-path /vast/palmer/apps/avx2/software/PLINK/1.9b_6.21-x86_64/plink --expression-matrix ",
                     dir, "gtexexpr.txt --exp-bfile ", dir, "gtexgeno --geno-bfile ",
                     dir, "gtexgeno --chr ", chr, " --out ", dir, "expr_mesc --tmp ",
                     dir, "exp_tmp"),
              paste0(mesc, " --compute-expscore-indiv --plink-path /vast/palmer/apps/avx2/software/PLINK/1.9b_6.21-x86_64/plink --expression-matrix ",
                     dir, "gtexexpr.txt --exp-bfile ", dir, "gtexgeno --geno-bfile ",
                     dir, "gtexgeno --chr ", chr, " --out ", dir, "exprun_mesc --tmp ",
                     dir, "expun_tmp", " --num-bins 1"))
write.table(EXPSCORE, file = paste0(dir2, "exprscore_mesc.sh"), quote = FALSE, row.names = FALSE, col.names = FALSE)

# expscore

step3 <- c(paste0(ldsc, " --bfile ", dir, "gtexgeno --l2 --annot ", dir, 
                  "exprun.annot.gz --thin-annot --ld-wind-cm 1 --out ", 
                  dir, "exprun.expscore; cd ", dir, "; ", 
                  gunzip, "exprun.expscore.l2.ldscore.gz;", 
                  " sed -i -e '1s/L2/Cis_herit_bin_1/' exprun.expscore.l2.ldscore; mv exprun.expscore.l2.ldscore exprun.expscore; ", 
                  gzip, "exprun.expscore"),
           paste0(ldsc, " --bfile ", dir, "gtexgeno --l2 --annot ", dir, 
                  "expr.annot.gz --thin-annot --ld-wind-cm 1 --out ", 
                  dir, "expr.expscore; cd ", dir, "; ",
                  gunzip, "expr.expscore.l2.ldscore.gz;", 
                  " sed -i -e '1s/Cis_herit_bin_1L2/ Cis_herit_bin_1/' -e '1s/Cis_herit_bin_2L2/ Cis_herit_bin_2/' -e '1s/Cis_herit_bin_3L2/ Cis_herit_bin_3/' -e '1s/Cis_herit_bin_4L2/ Cis_herit_bin_4/' -e '1s/Cis_herit_bin_5L2/ Cis_herit_bin_5/' ",
                  "expr.expscore.l2.ldscore; mv expr.expscore.l2.ldscore expr.expscore; ", 
                  gzip, "expr.expscore"))
write.table(step3, file = paste0(dir2, "exprscore_ldsc.sh"), quote = FALSE, row.names = FALSE, col.names = FALSE)


# dsq --job-file exprscore_ldsc.txt --partition day,pi_zhao --mem-per-cpu 5g -t 0-03:30:00 --mail-type ALL --output /dev/null