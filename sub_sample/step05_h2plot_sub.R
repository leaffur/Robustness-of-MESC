library(stringr)
library(ggplot2)
library(ggpubr)
library(grid)
chr = 1
G = 100
simul = 1
dir = "/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/chr1/s1_samplesize/"
dir2 <- "/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/s1_samplesize/"

system(paste0("mkdir -p ", paste0(dir, "figures")))
# step 05: read h2 and analysis

est <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("h_med", "simul_set", "Model", "estimand", "h_est", "h_se", "Sample"))

for (h_med in c(0.2)){
  for (tt in 1:G){
    h2mesc <- paste0(dir, "temp/unh2med_hm", h_med, "_t", tt, ".all.h2med")
    # h2smesc <- paste0(dir, "temp/h2med_hm", h_med, "_t", tt, ".all.h2med")
    A <- read.table(h2mesc, header = TRUE, row.names = 1)
    est <- rbind(est, 
                 cbind(c(h_med, h_med, h_med), 
                       rep(simul, 3),
                       c("MESC", "MESC", "MESC"),
                       c("h_med", "h_nmed", "h_g"),
                       A$Estimate,
                       A$SE.Estimate.,
                       c("True", "True", "True")))
    
    # # smesc
    # A <- read.table(h2smesc, header = TRUE, row.names = 1)
    # est <- rbind(est, 
    #              cbind(c(h_med, h_med, h_med),
    #                    rep(simul, 3),
    #                    c("SMESC", "SMESC", "SMESC"),
    #                    c("h_med", "h_nmed", "h_g"),
    #                    A$Estimate,
    #                    A$SE.Estimate.,
    #                    c("True", "True", "True")))

    for (i in c(500, 1000, 1500, 2000, 2500)){
      h2mesc_m <- paste0(dir, "temp/unh2med_hm", h_med, "_t", tt, "_mesc_", i, ".all.h2med")
      # h2smesc_m <- paste0(dir, "temp/h2med_hm", h_med, "_t", tt, "_mesc_", i, ".all.h2med")
      
      A <- read.table(h2mesc_m, header = TRUE, row.names = 1)
      est <- rbind(est,
                   cbind(c(h_med, h_med, h_med), 
                         rep(simul, 3),
                         c(paste0("MESC"), paste0("MESC"), paste0("MESC")),
                         c("h_med", "h_nmed", "h_g"),
                         A$Estimate,
                         A$SE.Estimate.,
                         rep(as.character(i), 3)))
      
      # # smesc
      # A <- read.table(h2smesc_m, header = TRUE, row.names = 1)
      # est <- rbind(est,
      #              cbind(c(h_med, h_med, h_med), 
      #                    rep(simul, 3),
      #                    c(paste0("SMESC"), paste0("SMESC"), paste0("SMESC")),
      #                    c("h_med", "h_nmed", "h_g"),
      #                    A$Estimate,
      #                    A$SE.Estimate.,
      #                    rep(as.character(i), 3)))
    }
  }
}

colnames(est) <- c("h_med", "simul_set", "Model", "estimand", "h_est", "h_se", "Sample")
est$h_est <- as.numeric(est$h_est)
est$h_se <- as.numeric(est$h_se)
est$Sample <- factor(est$Sample, levels = c("True", "2500", "2000", "1500", "1000", "500"), ordered = TRUE)

lines = data.frame(Model = c("MESC", "SMESC"), proportion = rep(0.2, 2))
lines2 = data.frame(Model = c("MESC", "SMESC"), proportion = rep(0.5, 2))
lines3 = data.frame(Model = c("MESC", "SMESC"), proportion = rep(0.3, 2))

## new 3.2 ### ------


g1 = ggplot(est[est$estimand == "h_med" & est$Model == "MESC", ], aes(x=Model, y=h_est, fill = Sample)) + 
  ylab(expression("Estimated h"[med]^2)) + xlab(expression("Model")) +
  geom_boxplot(linewidth = 1.2) + #facet_grid(~Model, scale='free_x') +
  stat_summary(fun = mean, geom = "point", shape=5, size=6, col = "red", position = position_dodge(0.75), show.legend = FALSE) +
  geom_hline(aes(yintercept = proportion), lines, lty = 2, col = "grey", linewidth = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.spacing.x = unit(0, "lines"), strip.text.x = element_blank(), legend.position="bottom",
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     text = element_text(size = 15))

g2 = ggplot(est[est$estimand == "h_nmed" & est$Model == "MESC", ], aes(x=Model, y=h_est, fill = Sample)) + 
  ylab(expression("Estimated h"[nmed]^2)) + xlab(expression("Model")) +
  geom_boxplot(linewidth = 1.2) + #facet_grid(~Model, scale='free_x') +
  stat_summary(fun = mean, geom = "point", shape=5, size=6, col = "red", position = position_dodge(0.75), show.legend = FALSE) +
  geom_hline(aes(yintercept = proportion), lines3, lty = 2, col = "grey", linewidth = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.spacing.x = unit(0, "lines"), strip.text.x = element_blank(), legend.position="bottom",
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15))

g3 = ggplot(est[est$estimand == "h_g" & est$Model == "MESC", ], aes(x=Model, y=h_est, fill = Sample)) + 
  ylab(expression("Estimated h"[total]^2)) + xlab(expression("Model")) +
  geom_boxplot(linewidth = 1.2) + #facet_grid(~Model, scale='free_x') +
  stat_summary(fun = mean, geom = "point", shape=5, size=6, col = "red", position = position_dodge(0.75), show.legend = FALSE) +
  geom_hline(aes(yintercept = proportion), lines2, lty = 2, col = "grey", linewidth = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.spacing.x = unit(0, "lines"), strip.text.x = element_blank(), legend.position="bottom",
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15))

top_row = ggarrange(g1, g2, ncol = 2, labels = c("a", "b"), legend = "top", common.legend = TRUE, 
                    font.label = list(size = 15, face = "bold"))
bottom_row = ggarrange(NULL, g3, NULL, ncol = 3, labels = c("", "c", ""), legend = "none", common.legend = TRUE, 
                       widths = c(1,2,1), font.label = list(size = 15, face = "bold"))
final_plot = ggarrange(top_row, bottom_row, ncol = 1)
p1 <- ggarrange(top_row, bottom_row, ncol = 1, heights = c(1, 0.9), legend = "top", common.legend = TRUE, 
                font.label = list(size = 15, face = "bold"))

ggexport(p1, filename = paste0(dir, "figures/ME_sub.eps"))
