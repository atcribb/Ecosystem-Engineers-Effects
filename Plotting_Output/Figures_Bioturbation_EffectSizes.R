#Authors: Alison Cribb, Will Gearty
#Summary: Output plots for bioturbation effect sizes:

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
library(tidyr)

effectsize_theme <- theme(
  panel.border=element_rect(fill=NA, colour='black'),
  legend.position="inside",
  legend.position.inside=c(0.87, 0.8),
  legend.title=element_blank(),
  legend.text=element_text(size=8),
  legend.background=element_rect(fill='white', colour=NA),
  plot.title=element_text(hjust=0.5),
  axis.text = element_text(color = "black", size=8),
  axis.title = element_text(size=10),
  axis.line.x = element_blank())

#==== load data ====#
load('Output/effectsizes_bioturbation_occsub.RData')
results_df <- bioturbation_results_df

results_df$stage <- factor(results_df$stage, levels=stage_names)
results_df$period <- factor(results_df$period, levels=period_names)

#===========  MAIN EFFECT SIZES FIGURES ===========

#Richness effect size
#restructure M1 vs M2 data
m1_richness <- results_df[,c(1:3,6:7)]
colnames(m1_richness) <- c('period', 'stage', 'mid_ma', 'richness', 'sd')
m1_richness$formations <- 'Bioturbators present'
m2_richness <- results_df[,c(1:3,8:9)]
colnames(m2_richness) <- c('period', 'stage', 'mid_ma', 'richness', 'sd')
m2_richness$formations <- 'Bioturbators absent'
compare_richness <- rbind(m1_richness, m2_richness)

compare.forms.cols <- c('#2759A9', '#8F9AAD')
compare_richness$formations <- factor(compare_richness$formations, 
                                      levels=c('Bioturbators present', 'Bioturbators absent'))

compare_richness_plot <- ggplot(data=compare_richness) +
  geom_rect(data=extinction_boundaries, aes(xmin=bottom, xmax=top, ymin=0, ymax=15),
            fill='#e76f51', alpha=0.15) +
  geom_errorbar(aes(x=mid_ma, ymin=richness-sd, ymax=richness+sd, 
                    color=formations),
                width=5) +
  geom_line(aes(x=mid_ma, y=richness, colour=formations)) +
  geom_point(aes(x=mid_ma, y=richness, fill=formations), colour='black', size=2, pch=21) +
  scale_fill_manual(values=compare.forms.cols) +
  scale_colour_manual(values=compare.forms.cols) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(0,15), name="Generic Richness") +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=520, y=14.5, label="A", size=3) +
  theme_classic() +
  effectsize_theme
compare_richness_plot  

#determine significance of effect size:
results_df$genrich_lowerbound <- results_df$HedgesG_genrich - results_df$g_genrich_sd #get lower bound to assess whether uncertainty overlaps with effect size thresholds  
results_df$genrich_significance <- rep(NA, nrow(results_df))
#assign effect size significance
results_df[which(abs(results_df$genrich_lowerbound) < 0.2),'genrich_significance'] <- 'no effect'
results_df[which( abs(results_df$genrich_lowerbound) >= 0.2 & abs(results_df$genrich_lowerbound) < 0.5),'genrich_significance'] <- 'weak effect'
results_df[which( abs(results_df$genrich_lowerbound) >= 0.5 & abs(results_df$genrich_lowerbound) < 0.8),'genrich_significance'] <- 'moderate effect'
results_df[which( abs(results_df$genrich_lowerbound) >= 0.8),'genrich_significance'] <- 'strong effect'
results_df$genrich_significance <- factor(results_df$genrich_significance, levels=c('no effect', 'weak effect', 'moderate effect', 'strong effect'))


effectsize_richness_plot <- ggplot(data=subset(results_df, !is.na(HedgesG_genrich))) +
  geom_rect(data=extinction_boundaries, aes(xmin=bottom, xmax=top, ymin=-1.2, ymax=3.2),
            fill='#e76f51', alpha=0.15) +
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
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill='none', size=guide_legend(byrow=TRUE)) +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=520, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme
effectsize_richness_plot  

Figure1 <- ggarrange2(compare_richness_plot, effectsize_richness_plot, ncol=1)
ggsave(filename='Draft/Figures/Figure1_bioturbation_effectsize_richness.pdf', Figure1,
       width=166, height=140, unit='mm')


#Shannon's Diversity
#restructure M1 vs M2 data
m1_H <- results_df[,c(1:3,12:13)]
colnames(m1_H) <- c('period', 'stage', 'mid_ma', 'H', 'sd')
m1_H$formations <- 'Bioturbators present'
m2_H <- results_df[,c(1:3,14:15)]
colnames(m2_H) <- c('period', 'stage', 'mid_ma', 'H', 'sd')
m2_H$formations <- 'Bioturbators absent'
compare_H <- rbind(m1_H, m2_H)

compare.forms.cols <- c('#2759A9', '#8F9AAD')
compare_H$formations <- factor(compare_H$formations, 
                                      levels=c('Bioturbators present', 'Bioturbators absent'))

compare_H_plot <- ggplot(data=compare_H) +
  geom_rect(data=extinction_boundaries, aes(xmin=bottom, xmax=top, ymin=0, ymax=2.8),
            fill='#e76f51', alpha=0.1) +
    geom_errorbar(aes(x=mid_ma, ymin=H-sd, ymax=H+sd, 
                    color=formations),
                width=5) +
  geom_line(aes(x=mid_ma, y=H, colour=formations)) +
  geom_point(aes(x=mid_ma, y=H, fill=formations), colour='black', size=2, pch=21) +
  scale_fill_manual(values=compare.forms.cols) +
  scale_colour_manual(values=compare.forms.cols) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(0,2.8), name="Shannon's Diversity (H)") +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=2.65, label="A", size=3) +
  theme_classic() +
  effectsize_theme
compare_H_plot  

#assign effect size significance thresholds
results_df$H_lowerbound <- abs(results_df$HedgesG_H) - results_df$g_H_sd
results_df$H_significance <- rep(NA, nrow(results_df))
results_df[which(abs(results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
results_df[which(abs(results_df$H_lowerbound) < 0.2 & is.na(results_df$H_significance)),'H_significance'] <- 'no effect'
results_df[which( abs(results_df$H_lowerbound) >= 0.2 & abs(results_df$H_lowerbound) < 0.5 & is.na(results_df$H_significance)),'H_significance'] <- 'weak effect'
results_df[which( abs(results_df$H_lowerbound) >= 0.5 & abs(results_df$H_lowerbound) < 0.8 & is.na(results_df$H_significance)),'H_significance'] <- 'moderate effect'
results_df[which( abs(results_df$H_lowerbound) >= 0.8 & is.na(results_df$H_significance)),'H_significance'] <- 'strong effect'
results_df$H_significance <- factor(results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))

effectsize_H_plot <- ggplot(data=results_df) +
  geom_rect(data=extinction_boundaries, aes(xmin=bottom, xmax=top, ymin=-1.2, ymax=3.2),
            fill='#e76f51', alpha=0.1) +
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
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill='none') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme
effectsize_H_plot  

Figure2_effects <- ggarrange2(compare_H_plot, effectsize_H_plot, ncol=1)

#====== Mass extinction subsets =====#
results_df$H_significance <- factor(results_df$H_significance, levels=c('strong effect', 'moderate effect', 'weak effect', 'no effect'))
endOrdovician_effects <- subset(results_df, stage %in% endOrdovician_stages)
endDevonian_effects <- subset(results_df, stage %in% endDevonian_stages)
endPermian_effects <- subset(results_df, stage %in% endPermian_stages)
endTriassic_effects <- subset(results_df, stage %in% endTriassic_stages)
KPg_effects <- subset(results_df, stage %in% KPg_stages)

EOMEfx <- ggplot(data=endOrdovician_effects) +
  geom_vline(aes(xintercept=440), col='#e76f51', size=0.5) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=0.25) +
  geom_line(aes(x=mid_ma, y=HedgesG_H), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H), colour='black', fill='gray40', shape=21, size=2) +
  scale_x_reverse(limits=c(endOrdovician_ages[1],endOrdovician_ages[2]), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill='none') +
  coord_geo(pos='bottom', dat='stages', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme
EOMEfx
ggsave(filename='Draft/Figures/assembly/Fig2_EOME.pdf', EOMEfx, width=75, height=60, unit='mm')

EDMEfx <- ggplot(data=endDevonian_effects) +
  geom_vline(aes(xintercept=365), col='#e76f51', size=0.5) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=0.25) +
  geom_line(aes(x=mid_ma, y=HedgesG_H), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H), colour='black', fill='gray40', shape=21, size=2) +  scale_x_reverse(limits=c(endDevonian_ages[1],endDevonian_ages[2]), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill='none') +
  coord_geo(pos='bottom', dat='stages', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme
EDMEfx
ggsave(filename='Draft/Figures/assembly/Fig2_EDME.pdf', EDMEfx, width=75, height=60, unit='mm')


EPMEfx <- ggplot(data=endPermian_effects) +
  geom_vline(aes(xintercept=252), col='#e76f51', size=0.5) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=0.25) +
  geom_line(aes(x=mid_ma, y=HedgesG_H), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H), colour='black', fill='gray40', shape=21, size=2) +  scale_x_reverse(limits=c(endDevonian_ages[1],endDevonian_ages[2]), name='Time (mya)') +
  scale_x_reverse(limits=c(endPermian_ages[1],endPermian_ages[2]), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill='none') +
  coord_geo(pos='bottom', dat='stages', size=2.8, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme
EPMEfx
ggsave(filename='Draft/Figures/assembly/Fig2_EPME.pdf', EPMEfx, width=75, height=60, unit='mm')


ETMEfx <- ggplot(data=endTriassic_effects) +
  geom_vline(aes(xintercept=201.3), col='#e76f51', size=0.5) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=0.25) +
  geom_line(aes(x=mid_ma, y=HedgesG_H), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H), colour='black', fill='gray40', shape=21, size=2) +  scale_x_reverse(limits=c(endDevonian_ages[1],endDevonian_ages[2]), name='Time (mya)') +
  scale_x_reverse(limits=c(endTriassic_ages[1],endTriassic_ages[2]), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill='none') +
  coord_geo(pos='bottom', dat='stages', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme
ETMEfx
ggsave(filename='Draft/Figures/assembly/Fig2_ETME.pdf', ETMEfx, width=75, height=60, unit='mm')


KPgMEfx <- ggplot(data=KPg_effects) +
  geom_vline(aes(xintercept=66), col='#e76f51', size=0.5) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=0.25) +
  geom_line(aes(x=mid_ma, y=HedgesG_H), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H), colour='black', fill='gray40', shape=21, size=2) +  scale_x_reverse(limits=c(endDevonian_ages[1],endDevonian_ages[2]), name='Time (mya)') +
  scale_x_reverse(limits=c(KPg_ages[1],KPg_ages[2]), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill='none') +
  coord_geo(pos='bottom', dat='stages', size=3, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme
KPgMEfx
ggsave(filename='Draft/Figures/assembly/Fig2_KPgMEfx.pdf', KPgMEfx, width=75, height=60, unit='mm')



