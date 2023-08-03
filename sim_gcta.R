args = (commandArgs(TRUE))
chr = as.numeric(args[[1]])
simul = as.numeric(args[[2]])
h_g = as.numeric(args[[3]])
N = as.numeric(args[[4]])
Ne = as.numeric(args[[5]])
dir = args[[6]]
h_med_fix = as.numeric(args[[7]])

.libPaths("/gpfs/gibbs/project/zhao/cl2384/R/4.2")
library(genio)

##################################

# dir <- paste0("/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/chr", chr, "/simulation_set", simul, "/")
system(paste0("mkdir -p ", dir))
print(dir)
dir1 <- paste0("/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/chr", chr, "/genotype/")
gcta <- "/gpfs/gibbs/pi/zhao/Softwares/gcta/gcta64"


## gtex sample

set.seed(724)
ukbb_fam <- read_fam(paste0(dir1, "ukbball"))
write.table(ukbb_fam[sample(1:nrow(ukbb_fam), Ne, replace = FALSE), 1:2],
            file = paste0(dir, "gtex_id.txt"),
            row.names = FALSE, quote = FALSE)

system(paste0("plink --bfile ", dir1, "ukbball --chr ", chr, " --keep ", dir, "gtex_id.txt --make-bed --out ", dir, "gtexgeno"))
# freq
system(paste0("plink -bfile ", dir, "gtexgeno --freq --out ", dir, "gtexgeno"))



# read data
gene_list <- read.table(paste0(dir1, "gene_list.txt"), header = TRUE)
snp_list <- read.table(paste0(dir1, "ukbball.snplist"))$V1
M <- length(snp_list)
G <- length(gene_list$GENE)


#### Standardization ####------------
gtex <- read_plink(paste0(dir, "gtexgeno"))
bed <- gtex$X
row_means <-  apply(bed, 1, function(y) y[is.na(y)] <- mean(na.omit(y)))
inds <- which(is.na(bed), arr.ind = TRUE)
bed[inds] <- row_means[inds[, 1]]
bed <- t(scale(t(bed)))

# simulation initialization
cis_position <- matrix("0", G, 5)
eqtl_g <- matrix(0, G, 5)
gtex_expr <- matrix(0, Ne, G)
rownames(gtex_expr) = gtex$fam$id
gtex_expr_true <- gtex_expr

h_cis <- rep(0.1, G)
if (simul %in% c(2, 3, 6, 7) ){
  h_cis <- 2^((1:G) / 200)
  h_cis <- h_cis / sum(h_cis) * G * 0.1
}
if (simul == 4){
  set.seed(152)
  h_cis <- rexp(G, 7) * 0.7
}


# expression imputation
set.seed(152)
for (kk in 1:G){
  # find cis-snps in 1 MB region
  inter <- gtex$bim$id[which((gtex$bim$pos <= gene_list$GENE_COORD[kk] + 5e5) & (gtex$bim$pos >= gene_list$GENE_COORD[kk] - 5e5))]
  cis_position[kk, ] <- sample(inter, 5)
  eqtl_g[kk, ] <- rnorm(5, 0, sqrt(h_cis[kk] / 20))
  eqtl_g[kk, 1] = rnorm(1, 0, sqrt(h_cis[kk] / 5 * 4))
  # h_cise ~ 0.8 * chi1 + 0.05 * chi4
  eqtl_g[kk, ] = scale(eqtl_g[kk, ]) * sqrt(h_cis[kk]) / 2
  
  # X_expr_true[, kk] <- colSums(X[cis_position[kk, ], ] * eqtl[kk, ])
  gtex_expr_true[, kk] <- colSums(bed[cis_position[kk, ], ] * eqtl_g[kk, ])
  gtex_expr[, kk] <- gtex_expr_true[, kk] + rnorm(Ne, 0, sqrt(1 - h_cis[kk]))
}

write.table(cbind(gene_list, t(gtex_expr_true)), file  = paste0(dir, "gtextrue.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(cbind(gene_list, t(gtex_expr)), file  = paste0(dir, "gtexexpr.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(cis_position, file  = paste0(dir, "cis_position"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(h_cis, file  = paste0(dir, "h_cis"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(eqtl_g, file  = paste0(dir, "eqtl_g"), quote = FALSE, row.names = FALSE, col.names = FALSE)



#### expscores ####--------------------
h_cise <- h_cis # + rnorm(G, 0, 0.01)
write.table(h_cise, file  = paste0(dir, "h_cise"), quote = FALSE, row.names = FALSE, col.names = FALSE)

gene_cat = NA
gene_cat[order(h_cise)] <- rep(1:5, each = ceiling(G / 5), length.out = G)

gannot <- matrix(0, G, 5)
rownames(gannot) = gene_list$GENE
colnames(gannot) = paste0("Cis_herit_bin_", 1:5)

annot_expr <- matrix(0, M, 5)
rownames(annot_expr) = snp_list
colnames(annot_expr) = paste0("Cis_herit_bin_", 1:5)

for (i in 1:G){
  annot_expr[cis_position[i, ], gene_cat[i]] = annot_expr[cis_position[i, ], gene_cat[i]] + 
    c(0.8, 0.05, 0.05, 0.05, 0.05) * h_cis[i]
  gannot[i, gene_cat[i]] = 1
}


numberg <- matrix(colSums(gannot), 1, 5)
avehis <- matrix(round(colSums(gannot * as.numeric(h_cise)) / numberg, digits = 4), 1, 5)
gannot <- data.frame(GENE = gene_list$GENE, gannot)

write.table(annot_expr, file  = paste0(dir, "expr.annot"), quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(gannot, file  = paste0(dir, "expr.gannot"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(numberg, file  = paste0(dir, "expr.G"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(avehis, file  = paste0(dir, "expr.ave_h2cis"), quote = FALSE, row.names = FALSE, col.names = FALSE)
system(paste0("gzip -f ", dir, "expr.annot"))

# unstratified

annot_exprun = rowSums(annot_expr)
write.table(annot_exprun, file  = paste0(dir, "exprun.annot"), quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data.frame(GENE = gene_list$GENE, "Cis_herit_bin_1" = rep(1, G)), file  = paste0(dir, "exprun.gannot"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(G, file  = paste0(dir, "exprun.G"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(round(mean(h_cise), digits = 4), file  = paste0(dir, "exprun.ave_h2cis"), quote = FALSE, row.names = FALSE, col.names = FALSE)
system(paste0("gzip -f ", dir, "exprun.annot"))



ukbb_fam <- read_fam(paste0(dir1, "ukbball"))
n.sim = 100
# cis_position <- read.table("/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/chr1/simulation_set1/cis_position")
system(paste0("mkdir -p ", dir, "phe_sum/"))
system(paste0("mkdir -p ", dir, "temp/"))

for (h_med in (c(0, 0.2, 0.4, 0.6, 0.8) * h_g)){
  
  if ((h_med_fix == 1) & (h_med != 0.2)) {
    next()
  }
  
  v_alpha = h_med / G / h_cis
  if (simul == 3){
    for (i in 1:5){
      Gi <- sum(gene_cat == i)
      v_alpha[order(h_cis)[(i * Gi):(i * Gi - Gi + 1)]] <- 2^((1:Gi) / 10)
    }
    v_alpha = h_med / sum(v_alpha * h_cis) * v_alpha
  }
  if (simul == 6){
    v_alpha = sample(v_alpha, G)
    v_alpha = h_med / sum(v_alpha * h_cis) * v_alpha
  }
  if (simul == 7){
    v_alpha = v_alpha[G:1]
    v_alpha = h_med / sum(v_alpha * h_cis) * v_alpha
  }
  
  v_gamma = rep((h_g - h_med) / M, M)
  if (simul == 5){
    v_gamma = rep((h_g - h_med) / M / 2, M)
    idx = which(snp_list %in% unique(c(cis_position)))
    v_gamma[idx] = v_gamma[idx] + (h_g - h_med) / length(idx) / 2
  }
  
  for (i in 1:n.sim){
    set.seed(i)
    # ukbb sample
    write.table(ukbb_fam[sample(1:nrow(ukbb_fam), N, replace = FALSE), 1:2],
                file = paste0(dir, "gwas_id.txt"),
                row.names = FALSE, quote = FALSE)
    system(paste0("plink --bfile ", dir1, "ukbball --chr ", chr, " --keep ", dir, "gwas_id.txt --make-bed --out ", dir, "ukbbgeno"))
    
    # effect sizes
    nonmed_gamma <- rnorm(M, 0, sqrt(v_gamma))
    gene_trt <- rnorm(G, 0, sqrt(v_alpha))
    # env_eps <- rnorm(N, 0, sqrt(1 - h_g))
    
    gamma <- nonmed_gamma
    names(gamma) <- snp_list
    for (kk in 1:G){
      gamma[cis_position[kk, ]] <- gamma[cis_position[kk, ]] + 
        gene_trt[kk] * c(rnorm(1, 0, sqrt(h_cis[kk] / 5 * 4)), rnorm(4, 0, sqrt(h_cis[kk] / 20)))
    }
    write.table(gamma, file  = paste0(dir, "phe_sum/beta_hm", h_med, "_t", i, ".snplist"), row.names = TRUE, col.names = FALSE, quote = FALSE)
    
    # pheno
    system(paste0(gcta, " --bfile ", dir, "ukbbgeno --simu-qt --simu-causal-loci ", 
                  dir, "phe_sum/beta_hm", h_med, "_t", i, ".snplist --simu-hsq 0.5 --out ", 
                  dir, "phe_sum/pheno_hm", h_med, "_t", i))
    
    # gwas sum
    system(paste0("plink --bfile ", dir, "ukbbgeno --pheno ", 
                  dir, "phe_sum/pheno_hm", h_med, "_t", i, ".phen --assoc --allow-no-sex --out ",
                  dir, "phe_sum/sum_hm", h_med, "_t", i))
    
    p3 <- read.table(paste0(dir, "phe_sum/sum_hm", h_med, "_t", i, ".qassoc"), header = TRUE)
    z1 = round(p3$BETA / p3$SE, digits = 2)
    
    sumstat <- data.frame(SNP = gtex$bim$id, A1 = gtex$bim$ref, A2 = gtex$bim$alt, N = rep(N, M), Z = z1)
    write.table(sumstat, file  = paste0(dir, "temp/sumstat_hm", h_med, "_t", i, ".sumstats"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    system(paste0("gzip -f ", dir, "temp/sumstat_hm", h_med, "_t", i, ".sumstats"))
  }
}

