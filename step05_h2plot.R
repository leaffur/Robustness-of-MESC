library(stringr)
library(ggplot2)
library(ggpubr)
library(grid)
chr = 1
G = 100
simul = 6
dir = "/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/chr1/simulation_set6/"

# step 05: read h2 and analysis

system(paste0("mkdir -p ", paste0(dir, "figures")))
# step 05: read h2 and analysis

est <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("h_med", "simul_set", "Model", "estimand", "h_est", "h_se"))

for (h_med in c(0, 0.1, 0.2, 0.3, 0.4)){
  for (tt in 1:G){
    h2mesc <- paste0(dir, "temp/unh2med_hm", h_med, "_t", tt, ".all.h2med")
    h2smesc <- paste0(dir, "temp/h2med_hm", h_med, "_t", tt, ".all.h2med")
    A <- read.table(h2mesc, header = TRUE, row.names = 1)
    est <- rbind(est, 
                 cbind(c(h_med, h_med, h_med), 
                       rep(simul, 3),
                       c("MESC", "MESC", "MESC"),
                       c("h_med", "h_nmed", "h_g"),
                       A$Estimate,
                       A$SE.Estimate.))
    
    # smesc
    A <- read.table(h2smesc, header = TRUE, row.names = 1)
    est <- rbind(est, 
                 cbind(c(h_med, h_med, h_med),
                       rep(simul, 3),
                       c("SMESC", "SMESC", "SMESC"),
                       c("h_med", "h_nmed", "h_g"),
                       A$Estimate,
                       A$SE.Estimate.))

    h2mesc_m <- paste0(dir, "temp/unh2med_hm", h_med, "_t", tt, "_mesc.all.h2med")
    h2smesc_m <- paste0(dir, "temp/h2med_hm", h_med, "_t", tt, "_mesc.all.h2med")
    
    A <- read.table(h2mesc_m, header = TRUE, row.names = 1)
    est <- rbind(est,
                 cbind(c(h_med, h_med, h_med), 
                       rep(simul, 3),
                       c(paste0("MESC_est"), paste0("MESC_est"), paste0("MESC_est")),
                       c("h_med", "h_nmed", "h_g"),
                       A$Estimate,
                       A$SE.Estimate.))
    
    # smesc
    A <- read.table(h2smesc_m, header = TRUE, row.names = 1)
    est <- rbind(est,
                 cbind(c(h_med, h_med, h_med), 
                       rep(simul, 3),
                       c(paste0("SMESC_est"), paste0("SMESC_est"), paste0("SMESC_est")),
                       c("h_med", "h_nmed", "h_g"),
                       A$Estimate,
                       A$SE.Estimate.))
  }
}

colnames(est) <- c("h_med", "simul_set", "Model", "estimand", "h_est", "h_se")
est$h_est <- as.numeric(est$h_est)
est$h_se <- as.numeric(est$h_se)
est$Model1 <- paste0(est$Model, est$simul_set)

lines = data.frame(h_med = paste0(c(0, 0.1, 0.2, 0.3, 0.4)), proportion = c(0, 0.1, 0.2, 0.3, 0.4))
lines2 = data.frame(h_med = paste0(c(0, 0.1, 0.2, 0.3, 0.4)), proportion = c(0.5, 0.5, 0.5, 0.5, 0.5))
lines3 = data.frame(h_med = paste0(c(0, 0.1, 0.2, 0.3, 0.4)), proportion = c(0.5, 0.4, 0.3, 0.2, 0.1))

## new 3.2 ### ------

kk = as.character(simul)
g1 = ggplot(est[est$estimand == "h_med" & est$simul_set == kk, ], aes(x=h_med, y=h_est, fill = Model)) + 
  ylab(expression("Estimated h"[med]^2)) + xlab(expression("True h"[med]^2)) +
  geom_boxplot(linewidth = 1.2) + facet_grid(~h_med, scale='free_x') +
  stat_summary(fun = mean, geom = "point", shape=5, size=6, col = "red", position = position_dodge(0.75), show.legend = FALSE) +
  geom_hline(aes(yintercept = proportion), lines, lty = 2, col = "grey", linewidth = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.spacing.x = unit(0, "lines"), strip.text.x = element_blank(), legend.position="bottom",
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     text = element_text(size = 15))

g2 = ggplot(est[est$estimand == "h_nmed" & est$simul_set == kk, ], aes(x=h_med, y=h_est, fill = Model)) + 
  ylab(expression("Estimated h"[nmed]^2)) + xlab(expression("True h"[med]^2)) +
  geom_boxplot(linewidth = 1.2) + facet_grid(~h_med, scale='free_x') +
  stat_summary(fun = mean, geom = "point", shape=5, size=6, col = "red", position = position_dodge(0.75), show.legend = FALSE) +
  geom_hline(aes(yintercept = proportion), lines3, lty = 2, col = "grey", linewidth = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.spacing.x = unit(0, "lines"), strip.text.x = element_blank(), legend.position="bottom",
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15))

g3 = ggplot(est[est$estimand == "h_g" & est$simul_set == kk, ], aes(x=h_med, y=h_est, fill = Model)) + 
  ylab(expression("Estimated h"[total]^2)) + xlab(expression("True h"[med]^2)) +
  geom_boxplot(linewidth = 1.2) + facet_grid(~h_med, scale='free_x') +
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
ggexport(p1, filename = paste0(dir, "figures/ME_set", kk, ".eps"))
  



#   3.3 ### ----------------

kk = as.character(simul)
g1 = ggplot(est[est$estimand == "h_med" & est$simul_set == kk & (est$Model == "MESC" | est$Model == "SMESC"), ], aes(x=h_med, y=h_est, fill = Model)) + 
  ylab(expression("Estimated h"[med]^2)) + xlab(expression("True h"[med]^2)) +
  geom_boxplot(linewidth = 1.2) + facet_grid(~h_med, scale='free_x') +
  stat_summary(fun = mean, geom = "point", shape=5, size=6, col = "red", position = position_dodge(0.75), show.legend = FALSE) +
  geom_hline(aes(yintercept = proportion), lines, lty = 2, col = "grey", linewidth = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.spacing.x = unit(0, "lines"), strip.text.x = element_blank(), legend.position="bottom",
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
                     text = element_text(size = 15))

g2 = ggplot(est[est$estimand == "h_nmed" & est$simul_set == kk & (est$Model == "MESC" | est$Model == "SMESC"), ], aes(x=h_med, y=h_est, fill = Model)) + 
  ylab(expression("Estimated h"[nmed]^2)) + xlab(expression("True h"[med]^2)) +
  geom_boxplot(linewidth = 1.2) + facet_grid(~h_med, scale='free_x') +
  stat_summary(fun = mean, geom = "point", shape=5, size=6, col = "red", position = position_dodge(0.75), show.legend = FALSE) +
  geom_hline(aes(yintercept = proportion), lines3, lty = 2, col = "grey", linewidth = 2) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.spacing.x = unit(0, "lines"), strip.text.x = element_blank(), legend.position="bottom",
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 15))

g3 = ggplot(est[est$estimand == "h_g" & est$simul_set == kk & (est$Model == "MESC" | est$Model == "SMESC"), ], aes(x=h_med, y=h_est, fill = Model)) + 
  ylab(expression("Estimated h"[total]^2)) + xlab(expression("True h"[med]^2)) +
  geom_boxplot(linewidth = 1.2) + facet_grid(~h_med, scale='free_x') +
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
ggexport(p1, filename = paste0(dir, "figures/us_set", kk, ".eps"))
