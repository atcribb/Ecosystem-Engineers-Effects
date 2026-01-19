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
  legend.title=element_blank(),
  legend.text=element_text(size=6),
  legend.background=element_rect(fill='white', colour=NA),
  axis.title.x=element_text(size=6),
  axis.title.y=element_text(size=6),
  plot.title=element_text(hjust=0.5),
  axis.text = element_text(color = "black", size=6),
  axis.line.x = element_blank())

#mass extinction intervals
endOrdovician_ages <- c(subset(stages, stage=='Katian')$bottom, subset(stages, stage=='Telychian')$top)
endOrdovician_stages <- stages$stage[19:23]
endDevonian_ages   <- c(subset(stages, stage=='Givetian')$bottom, subset(stages, stage=='Visean')$top)
endDevonian_stages <- stages$stage[33:37]
endPermian_ages    <- c(subset(stages, stage=='Capitanian')$bottom, subset(stages, stage=='Ladinian')$top)
endPermian_stages <- stages$stage[49:55]
endTriassic_ages    <- c(subset(stages, stage=='Rhaetian')$bottom, subset(stages, stage=='Pliensbachian')$top)
endTriassic_stages  <- stages$stage[58:61]
KPg_ages           <- c(subset(stages, stage=='Campanian')$bottom, subset(stages, stage=='Selandian-Thanetian')$top)
KPg_stages        <- stages$stage[80:83]
extinction_boundaries <- as.data.frame(rbind(endOrdovician_ages, endDevonian_ages, endPermian_ages, endTriassic_ages, KPg_ages))
colnames(extinction_boundaries) <- c('bottom', 'top')

#==== load data ====#
load('Output/effectsizes_bioturbation_occsub.RData')
load('Output/effectsizes_reefs_occsub.RData')

bioturbation_results_df$stage <- factor(results_df$stage, levels=stage_names)
reefs_results_df$period <- factor(results_df$period, levels=period_names)

#combine
biots_subset <- as.data.frame(cbind(bioturbation_results_df$stage, as.numeric(bioturbation_results_df$mid_ma),
                      as.numeric(bioturbation_results_df$HedgesG_H), as.numeric(bioturbation_results_df$g_H_sd),
                      rep('bioturbators', nrow(bioturbation_results_df))))
colnames(biots_subset) <- c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd', 'EE')
reefs_subset <- as.data.frame(cbind(reefs_results_df$stage, as.numeric(reefs_results_df$mid_ma),
                              as.numeric(reefs_results_df$HedgesG_H), as.numeric(reefs_results_df$g_H_sd),
                              rep('reefs', nrow(reefs_results_df))))
colnames(reefs_subset) <- colnames(biots_subset)                      
MEfx_results <- as.data.frame(rbind(biots_subset, reefs_subset))
MEfx_results$mid_ma <- as.numeric(MEfx_results$mid_ma)
MEfx_results$HedgesG_H <- as.numeric(MEfx_results$HedgesG_H)
MEfx_results$g_H_sd <- as.numeric(MEfx_results$g_H_sd)


#====== Mass extinction FX =====#
endOrdovician_effects <- subset(MEfx_results, stage %in% endOrdovician_stages)
endDevonian_effects <- subset(MEfx_results, stage %in% endDevonian_stages)
endPermian_effects <- subset(MEfx_results, stage %in% endPermian_stages)
endTriassic_effects <- subset(MEfx_results, stage %in% endTriassic_stages)
KPg_effects <- subset(MEfx_results, stage %in% KPg_stages)

EOMEfx <- ggplot(data=endOrdovician_effects) +
  geom_vline(aes(xintercept=444), col='#e76f51', size=0.5) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.25, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.25, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.25, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=1) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=EE), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=EE), colour='black', shape=21, size=1.5) +
  scale_fill_manual(values=c('#1491FF', '#D56482')) +
  scale_color_manual(values=c('#1491FF', '#D56482')) +
  scale_x_reverse(limits=c(endOrdovician_ages[1],endOrdovician_ages[2]), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  coord_geo(pos='bottom', dat='stages', size=2, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme +
  theme(
    legend.position="inside",
    legend.position.inside=c(0.2, 0.8))
EOMEfx

EDMEfx <- ggplot(data=endDevonian_effects) +
  geom_vline(aes(xintercept=371), col='#e76f51', size=0.5) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.25, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.25, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.25, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=2.5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=EE), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=EE), colour='black', shape=21, size=1.5) +
  scale_fill_manual(values=c('#1491FF', '#D56482')) +
  scale_color_manual(values=c('#1491FF', '#D56482')) +
  scale_x_reverse(limits=c(endDevonian_ages[1],endDevonian_ages[2]), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  coord_geo(pos='bottom', dat='stages', size=2, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme +
  theme(legend.position='none')
EDMEfx

EPMEfx <- ggplot(data=endPermian_effects) +
  geom_vline(aes(xintercept=252), col='#e76f51', size=0.5) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.25, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.25, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.25, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=1) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=EE), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=EE), colour='black', shape=21, size=1.5) +
  scale_fill_manual(values=c('#1491FF', '#D56482')) +
  scale_color_manual(values=c('#1491FF', '#D56482')) +
  scale_x_reverse(limits=c(endPermian_ages[1],endPermian_ages[2]), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  coord_geo(pos='bottom', dat='stages', size=1.8, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme  +
  theme(legend.position='none')
EPMEfx

ETMEfx <- ggplot(data=endTriassic_effects) +
  geom_vline(aes(xintercept=201.3), col='#e76f51', size=0.5) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.25, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.25, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.25, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=1) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=EE), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=EE), colour='black', shape=21, size=1.5) +
  scale_fill_manual(values=c('#1491FF', '#D56482')) +
  scale_color_manual(values=c('#1491FF', '#D56482')) +
  scale_x_reverse(limits=c(endTriassic_ages[1],endTriassic_ages[2]), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  coord_geo(pos='bottom', dat='stages', size=2, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme  +
  theme(legend.position='none')
ETMEfx

KPgMEfx <- ggplot(data=KPg_effects) +
  geom_vline(aes(xintercept=66), col='#e76f51', size=0.5) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.25, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.25, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.25, color='gray30') +
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=1) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=EE), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=EE), colour='black', shape=21, size=1.5) +
  scale_fill_manual(values=c('#1491FF', '#D56482')) +
  scale_color_manual(values=c('#1491FF', '#D56482')) +
  scale_x_reverse(limits=c(KPg_ages[1],KPg_ages[2]), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  coord_geo(pos='bottom', dat='stages', size=2, abbrv=TRUE, height=unit(1,'line')) +
  annotate("text", x=525, y=3, label="B", size=3) +
  theme_classic() +
  effectsize_theme +
  theme(legend.position='none')
KPgMEfx

ME_effects <- ggarrange(EOMEfx, EDMEfx, EPMEfx, ETMEfx, KPgMEfx, ncol=2)


