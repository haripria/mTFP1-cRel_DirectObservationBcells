library(ggplot2)
library(dplyr)
library(ggpubr)
library(moments)
library(EnvStats)
# preStimulusMetabolites1 <- read.csv(file="/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/cRel\ Abundance\ Project/results_distribute_all/allPreStimulusMetabolites.CSV", header = TRUE)
# preStimulusMetabolites2 <- read.csv(file="/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/cRel\ Abundance\ Project/results_IkBeKO/allPreStimulusMetabolites.CSV", header = TRUE)
preStimulusMetabolites1 <- read.csv(file="/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/cRel\ Abundance\ Project/results_cRel_steady/WT+delta_cRel5_adjusted.CSV", header = TRUE)
preStimulusMetabolites2 <- read.csv(file="/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/cRel\ Abundance\ Project/results_cRel_steady/epsilonKO+delta_cRel5_adjusted.CSV", header = TRUE)
# WT
preStimulusMetabolites1$total_RelA <- preStimulusMetabolites1$RelA + preStimulusMetabolites1$RelAn + 
  preStimulusMetabolites1$AA + preStimulusMetabolites1$AAn + 
  preStimulusMetabolites1$IkBaAA + preStimulusMetabolites1$IkBaAAn + preStimulusMetabolites1$IkBbAA + preStimulusMetabolites1$IkBbAAn + 
  preStimulusMetabolites1$IkBeAA + preStimulusMetabolites1$IkBeAAn + preStimulusMetabolites1$IkBdAA + preStimulusMetabolites1$IkBdAAn + 
  preStimulusMetabolites1$A50 + preStimulusMetabolites1$A50n + 
  preStimulusMetabolites1$IkBaA50 + preStimulusMetabolites1$IkBaA50n + preStimulusMetabolites1$IkBbA50 + preStimulusMetabolites1$IkBbA50n + 
  preStimulusMetabolites1$IkBeA50 + preStimulusMetabolites1$IkBeA50n + preStimulusMetabolites1$IkBdA50 + preStimulusMetabolites1$IkBdA50n
preStimulusMetabolites1$total_cRel <- preStimulusMetabolites1$cRel + preStimulusMetabolites1$cReln + 
  preStimulusMetabolites1$C50 + preStimulusMetabolites1$C50n + 
  preStimulusMetabolites1$IkBaC50 + preStimulusMetabolites1$IkBaC50n + preStimulusMetabolites1$IkBbC50 + preStimulusMetabolites1$IkBbC50n + 
  preStimulusMetabolites1$IkBeC50 + preStimulusMetabolites1$IkBeC50n + preStimulusMetabolites1$IkBdC50 + preStimulusMetabolites1$IkBdC50n

# epslion-KO
preStimulusMetabolites2$total_RelA <- preStimulusMetabolites2$RelA + preStimulusMetabolites2$RelAn + 
  preStimulusMetabolites2$AA + preStimulusMetabolites2$AAn + 
  preStimulusMetabolites2$IkBaAA + preStimulusMetabolites2$IkBaAAn + preStimulusMetabolites2$IkBbAA + preStimulusMetabolites2$IkBbAAn + 
  preStimulusMetabolites2$IkBeAA + preStimulusMetabolites2$IkBeAAn + preStimulusMetabolites2$IkBdAA + preStimulusMetabolites2$IkBdAAn + 
  preStimulusMetabolites2$A50 + preStimulusMetabolites2$A50n + 
  preStimulusMetabolites2$IkBaA50 + preStimulusMetabolites2$IkBaA50n + preStimulusMetabolites2$IkBbA50 + preStimulusMetabolites2$IkBbA50n + 
  preStimulusMetabolites2$IkBeA50 + preStimulusMetabolites2$IkBeA50n + preStimulusMetabolites2$IkBdA50 + preStimulusMetabolites2$IkBdA50n
preStimulusMetabolites2$total_cRel <- preStimulusMetabolites2$cRel + preStimulusMetabolites2$cReln + 
  preStimulusMetabolites2$C50 + preStimulusMetabolites2$C50n + 
  preStimulusMetabolites2$IkBaC50 + preStimulusMetabolites2$IkBaC50n + preStimulusMetabolites2$IkBbC50 + preStimulusMetabolites2$IkBbC50n + 
  preStimulusMetabolites2$IkBeC50 + preStimulusMetabolites2$IkBeC50n + preStimulusMetabolites2$IkBdC50 + preStimulusMetabolites2$IkBdC50n

preStimulusMetabolites <- rbind(preStimulusMetabolites1, preStimulusMetabolites2)
preStimulusMetabolites$condition <- rep(c("WT", "EKO"), each=2000)
preStimulusMetabolites$condition <- factor(preStimulusMetabolites$condition , levels=c("WT", "EKO") )

# Initial abundance distribution of cRel
ggplot(preStimulusMetabolites, aes(x=total_cRel, fill=condition, color=condition, linetype=condition)) +
  geom_density(alpha=.2, size=0.7, adjust = 2) +
  scale_color_manual(values=c("grey30", "#FF0066")) + scale_fill_manual(values=c("grey30", "#FF0066")) + 
  scale_linetype_manual(values=c("dotted", "solid")) + scale_x_log10(limits = c(20,200)) + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("density")


write.csv(preStimulusMetabolites1$total_cRel, "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/cRel\ Abundance\ Project/results_cRel_steady/WTcRel.CSV", row.names=FALSE, col.names=TRUE)
write.csv(preStimulusMetabolites1$total_cRel, "/Users/helenhuang/Documents/1st\ Year\ PhD/Hoffmann\ Lab/cRel\ Abundance\ Project/results_cRel_steady/KOcRel.CSV", row.names=FALSE, col.names=TRUE)


skewness(preStimulusMetabolites$total_RelA)
kurtosis(preStimulusMetabolites$total_RelA)
skewness(preStimulusMetabolites$total_cRel)
kurtosis(preStimulusMetabolites$total_cRel)
mean(preStimulusMetabolites$total_cRel)
sd(preStimulusMetabolites$total_cRel)
