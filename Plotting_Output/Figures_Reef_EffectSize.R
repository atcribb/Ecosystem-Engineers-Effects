#Authors: Alison Cribb, Will Gearty
#Summary: Final analyses and output plots for reef effect sizes 

# clear old data
rm(list = ls())

#==== load packages ====#
library(divDyn)
data("stages", package="divDyn")
stage_names <- stages$stage[4:95]
stage_mids <- stages$mid[4:95]
period_names <- unique(stages[which(stages$stage %in% stage_names), 'system'])
period.cols <- unique(stages[which(stages$stage %in% stage_names), 'systemCol'])
library(ggplot2)
library(deeptime)
library(metafor)
library(tidyr)
library(RColorBrewer)

effectsize_theme <- theme(
  panel.border=element_rect(fill=NA, colour='black'),
  legend.position="inside",
  legend.position.inside=c(0.87, 0.8),
  legend.title=element_blank(),
  legend.text=element_text(size=8),
  legend.background=element_rect(fill=NA, colour=NA),
  plot.title=element_text(hjust=0.5),
  axis.text = element_text(color = "black", size=8),
  axis.title = element_text(size=10),
  axis.line.x = element_blank())

#==== load data ====#
load('Output/effectsizes_reefs_occsub.RData')
results_df <- reef_results_df

results_df$stage <- factor(results_df$stage, levels=stage_names)
results_df$period <- factor(results_df$period, levels=period_names)

#===========  MAIN EFFECT SIZES FIGURES ===========

#Richness effect size
#restructure M1 vs M2 data
m1_richness <- results_df[,c(1:3,6:7)]
colnames(m1_richness) <- c('period', 'stage', 'mid_ma', 'richness', 'sd')
m1_richness$formations <- 'Reef-builders present'
m2_richness <- results_df[,c(1:3,8:9)]
colnames(m2_richness) <- c('period', 'stage', 'mid_ma', 'richness', 'sd')
m2_richness$formations <- 'Reef-builders absent'
compare_richness <- rbind(m1_richness, m2_richness)

# these are pretty similar for colorblind folks, might want to change one of them -WG
compare.forms.cols <- c('#A03544', '#B19398')
compare_richness$formations <- factor(compare_richness$formations, 
                                      levels=c('Reef-builders present', 'Reef-builders absent'))

compare_richness_plot <- ggplot(data=compare_richness) +
  #geom_ribbon(aes(x=mid_ma, ymin=richness-sd, ymax=richness+sd, fill=formations), alpha=0.5) +
  geom_errorbar(aes(x=mid_ma, ymin=richness-sd, ymax=richness+sd, 
                    color=formations),
                width=5) +
  geom_line(aes(x=mid_ma, y=richness, colour=formations)) +
  geom_point(aes(x=mid_ma, y=richness, fill=formations), colour='black', size=2.5, pch=21) +
  scale_fill_manual(values=compare.forms.cols) +
  scale_colour_manual(values=compare.forms.cols) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(0,14), name="Generic Richness") +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=530, y=13.5, label="A", size=3) +
  theme_classic() +
  effectsize_theme
compare_richness_plot  

#assign effect size significance thresholds
results_df$genrich_lowerbound <- abs(results_df$HedgesG_genrich) - results_df$g_genrich_sd #find where lower bound does not overlap with the lower significance threshold
results_df$genrich_significance <- rep(NA, nrow(results_df))
results_df[which(abs(results_df$genrich_lowerbound) < 0.2),'genrich_significance'] <- 'no effect'
results_df[which( abs(results_df$genrich_lowerbound) >= 0.2 & abs(results_df$genrich_lowerbound) < 0.5),'genrich_significance'] <- 'weak effect'
results_df[which( abs(results_df$genrich_lowerbound) >= 0.5 & abs(results_df$genrich_lowerbound) < 0.8),'genrich_significance'] <- 'moderate effect'
results_df[which( abs(results_df$genrich_lowerbound) >= 0.8),'genrich_significance'] <- 'strong effect'
results_df$genrich_significance <- factor(results_df$genrich_significance, levels=c('no effect', 'weak effect', 'moderate effect', 'strong effect'))


effectsize_richness_plot <- ggplot(data=subset(results_df, !is.na(HedgesG_genrich))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_genrich-g_genrich_sd,
                    ymax=HedgesG_genrich+g_genrich_sd),
                width=8) +
  geom_line(aes(x=mid_ma, y=HedgesG_genrich), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_genrich, size=genrich_significance), fill='white', shape=21) +
  geom_point(aes(x=mid_ma, y=HedgesG_genrich, fill=period, alpha=genrich_significance, size=genrich_significance), shape=21) +
  scale_size_manual(values=c(1, 1.5, 2, 2.75)) +
  scale_alpha_manual(values=c(0.1, 0.3, 0.7, 1.0)) +
  scale_fill_manual(values=period.cols) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,2.5), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill='none') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=530, y=2.3, label="B", size=3) +
  theme_classic() +
  effectsize_theme
effectsize_richness_plot  

Figure3 <- ggarrange2(compare_richness_plot, effectsize_richness_plot, ncol=1)
ggsave(filename='Draft/Figures/Figure3_reef_effectsize_richness.pdf', Figure3,
       width=166, height=140, unit='mm')


#Shannon's Diversity
#restructure M1 vs M2 data
m1_H <- results_df[,c(1:3,12:13)]
colnames(m1_H) <- c('period', 'stage', 'mid_ma', 'H', 'sd')
m1_H$formations <- 'Reef-builders present'
m2_H <- results_df[,c(1:3,14:15)]
colnames(m2_H) <- c('period', 'stage', 'mid_ma', 'H', 'sd')
m2_H$formations <- 'Reef-builders absent'
compare_H <- rbind(m1_H, m2_H)

compare.forms.cols <- c('#A03544', '#B19398')
compare_H$formations <- factor(compare_H$formations, 
                               levels=c('Reef-builders present', 'Reef-builders absent'))

compare_H_plot <- ggplot(data=compare_H) +
  geom_errorbar(aes(x=mid_ma, ymin=H-sd, ymax=H+sd, 
                    color=formations),
                width=5) +
  geom_line(aes(x=mid_ma, y=H, colour=formations)) +
  geom_point(aes(x=mid_ma, y=H, fill=formations), colour='black', size=2.5, pch=21) +
  scale_fill_manual(values=compare.forms.cols) +
  scale_colour_manual(values=compare.forms.cols) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(0,2.5), name="Shannon's Diversity (H)") +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=530, y=2.4, label="A", size=3) +
  theme_classic() +
  effectsize_theme
compare_H_plot  


#Assign effect size significance thresholds
results_df$H_lowerbound <- abs(results_df$HedgesG_H) - results_df$g_H_sd
results_df$H_significance <- rep(NA, nrow(results_df))
results_df[which(abs(results_df$H_lowerbound) < 0.2),'H_significance'] <- 'no effect'
results_df[which( abs(results_df$H_lowerbound) >= 0.2 & abs(results_df$H_lowerbound) < 0.5),'H_significance'] <- 'weak effect'
results_df[which( abs(results_df$H_lowerbound) >= 0.5 & abs(results_df$H_lowerbound) < 0.8),'H_significance'] <- 'moderate effect'
results_df[which( abs(results_df$H_lowerbound) >= 0.8),'H_significance'] <- 'strong effect'
results_df$H_significance <- factor(results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))

effectsize_H_plot <- ggplot(data=subset(results_df, !is.na(HedgesG_H))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=8) +
  geom_line(aes(x=mid_ma, y=HedgesG_H), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, size=H_significance), fill='white', shape=21) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=period, alpha=H_significance, size=H_significance), shape=21) +
  scale_size_manual(values=c(1, 1.5, 2, 2.75)) +
  scale_alpha_manual(values=c(0.1, 0.3, 0.7, 1.0)) +
  scale_fill_manual(values=period.cols) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,2.5), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill='none') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=530, y=2.2, label="B", size=3) +
  theme_classic() +
  effectsize_theme
effectsize_H_plot  

Figure4 <- ggarrange2(compare_H_plot, effectsize_H_plot, ncol=1)
ggsave(filename='Draft/Figures/Figure4_reef_effectsize_H.pdf', Figure4,
       width=166, height=140, unit='mm')
