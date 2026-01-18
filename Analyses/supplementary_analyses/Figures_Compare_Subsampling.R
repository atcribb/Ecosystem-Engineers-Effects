#Summary: Compare subsampling results from formation-based and equal area grid cell-based spatial subsampling approaches
#Author: Alison T. Cribb
#Email: A.T.Cribb@soton.ac.uk

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

#Load data
load("Output/bioturbators_spatialsubsample.RData") #Equal area: Bioturbators
load("Output/reefbuilders_spatialsubsample.RData") #Equal area: Reef-builders 
load("Output/effectsizes_bioturbators_occsub.RData") #Formations: Bioturbators; change to collsub or noformsub if that's what you're using
load("Output/effectsizes_reefs_occsub.RData") #Formations: Reef-builders; change to collsub or noformsub if that's what you're using


#Bioturbation figure 
#determine significance of effect size:
#assign effect size significance thresholds
bioturbation_occsub_results_df$H_lowerbound <- abs(bioturbation_occsub_results_df$HedgesG_H) - bioturbation_occsub_results_df$g_H_sd
bioturbation_occsub_results_df$H_significance <- rep(NA, nrow(bioturbation_occsub_results_df))
bioturbation_occsub_results_df[which(abs(bioturbation_occsub_results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
bioturbation_occsub_results_df[which(abs(bioturbation_occsub_results_df$H_lowerbound) < 0.2 & is.na(bioturbation_occsub_results_df$H_significance)),'H_significance'] <- 'no effect'
bioturbation_occsub_results_df[which( abs(bioturbation_occsub_results_df$H_lowerbound) >= 0.2 & abs(bioturbation_occsub_results_df$H_lowerbound) < 0.5 & is.na(bioturbation_occsub_results_df$H_significance)),'H_significance'] <- 'weak effect'
bioturbation_occsub_results_df[which( abs(bioturbation_occsub_results_df$H_lowerbound) >= 0.5 & abs(bioturbation_occsub_results_df$H_lowerbound) < 0.8 & is.na(bioturbation_occsub_results_df$H_significance)),'H_significance'] <- 'moderate effect'
bioturbation_occsub_results_df[which( abs(bioturbation_occsub_results_df$H_lowerbound) >= 0.8 & is.na(bioturbation_occsub_results_df$H_significance)),'H_significance'] <- 'strong effect'
bioturbation_occsub_results_df$H_significance <- factor(bioturbation_occsub_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
bioturbation_occsub_results_df$method <- c('formations')

bioturbators_spatial_results_df$H_lowerbound <- abs(bioturbators_spatial_results_df$HedgesG_H) - bioturbators_spatial_results_df$g_H_sd
bioturbators_spatial_results_df$H_significance <- rep(NA, nrow(bioturbators_spatial_results_df))
bioturbators_spatial_results_df[which(abs(bioturbators_spatial_results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
bioturbators_spatial_results_df[which(abs(bioturbators_spatial_results_df$H_lowerbound) < 0.2 & is.na(bioturbators_spatial_results_df$H_significance)),'H_significance'] <- 'no effect'
bioturbators_spatial_results_df[which( abs(bioturbators_spatial_results_df$H_lowerbound) >= 0.2 & abs(bioturbators_spatial_results_df$H_lowerbound) < 0.5 & is.na(bioturbators_spatial_results_df$H_significance)),'H_significance'] <- 'weak effect'
bioturbators_spatial_results_df[which( abs(bioturbators_spatial_results_df$H_lowerbound) >= 0.5 & abs(bioturbators_spatial_results_df$H_lowerbound) < 0.8 & is.na(bioturbators_spatial_results_df$H_significance)),'H_significance'] <- 'moderate effect'
bioturbators_spatial_results_df[which( abs(bioturbators_spatial_results_df$H_lowerbound) >= 0.8 & is.na(bioturbators_spatial_results_df$H_significance)),'H_significance'] <- 'strong effect'
bioturbators_spatial_results_df$H_significance <- factor(bioturbators_spatial_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
bioturbators_spatial_results_df$method <- c('equal area grid cells')

bioturbators_compare <- rbind(bioturbation_occsub_results_df, bioturbators_spatial_results_df)

#determine significance of effect size:
reefs_occsub_results_df$H_lowerbound <- abs(reefs_occsub_results_df$HedgesG_H) - reefs_occsub_results_df$g_H_sd
reefs_occsub_results_df$H_significance <- rep(NA, nrow(reefs_occsub_results_df))
reefs_occsub_results_df[which(abs(reefs_occsub_results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
reefs_occsub_results_df[which(abs(reefs_occsub_results_df$H_lowerbound) < 0.2 & is.na(reefs_occsub_results_df$H_significance)),'H_significance'] <- 'no effect'
reefs_occsub_results_df[which( abs(reefs_occsub_results_df$H_lowerbound) >= 0.2 & abs(reefs_occsub_results_df$H_lowerbound) < 0.5 & is.na(reefs_occsub_results_df$H_significance)),'H_significance'] <- 'weak effect'
reefs_occsub_results_df[which( abs(reefs_occsub_results_df$H_lowerbound) >= 0.5 & abs(reefs_occsub_results_df$H_lowerbound) < 0.8 & is.na(reefs_occsub_results_df$H_significance)),'H_significance'] <- 'moderate effect'
reefs_occsub_results_df[which( abs(reefs_occsub_results_df$H_lowerbound) >= 0.8 & is.na(reefs_occsub_results_df$H_significance)),'H_significance'] <- 'strong effect'
reefs_occsub_results_df$H_significance <- factor(reefs_occsub_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
reefs_occsub_results_df$method <- c('formations')

reefbuilders_spatial_results_df$H_lowerbound <- abs(reefbuilders_spatial_results_df$HedgesG_H) - reefbuilders_spatial_results_df$g_H_sd
reefbuilders_spatial_results_df$H_significance <- rep(NA, nrow(reefbuilders_spatial_results_df))
reefbuilders_spatial_results_df[which(abs(reefbuilders_spatial_results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
reefbuilders_spatial_results_df[which(abs(reefbuilders_spatial_results_df$H_lowerbound) < 0.2 & is.na(reefbuilders_spatial_results_df$H_significance)),'H_significance'] <- 'no effect'
reefbuilders_spatial_results_df[which( abs(reefbuilders_spatial_results_df$H_lowerbound) >= 0.2 & abs(reefbuilders_spatial_results_df$H_lowerbound) < 0.5 & is.na(reefbuilders_spatial_results_df$H_significance)),'H_significance'] <- 'weak effect'
reefbuilders_spatial_results_df[which( abs(reefbuilders_spatial_results_df$H_lowerbound) >= 0.5 & abs(reefbuilders_spatial_results_df$H_lowerbound) < 0.8 & is.na(reefbuilders_spatial_results_df$H_significance)),'H_significance'] <- 'moderate effect'
reefbuilders_spatial_results_df[which( abs(reefbuilders_spatial_results_df$H_lowerbound) >= 0.8 & is.na(reefbuilders_spatial_results_df$H_significance)),'H_significance'] <- 'strong effect'
reefbuilders_spatial_results_df$H_significance <- factor(reefbuilders_spatial_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
reefbuilders_spatial_results_df$method <- c('equal area grid cells')

reefs_compare <- rbind(reefs_occsub_results_df, reefbuilders_spatial_results_df)



fig_compare_bioturbators <- ggplot(data=bioturbators_compare) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=method),
                alpha=0.8) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, color=method), linewidth=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, size=H_significance, fill=method), alpha=0.85, shape=21) +
  scale_color_manual(values=c('black', 'black')) +
  scale_fill_manual(values=c('#e07a5f', '#81b29a')) +
  scale_size_manual(values=c(1, 1.5, 2, 2.75)) +
  scale_alpha_manual(values=c(0.1, 0.3, 0.7, 1.0)) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  annotate("text", x=500, y=3, label="Bioturbators", size=3) +
  guides(fill=guide_legend(byrow=TRUE), size='none') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme
fig_compare_bioturbators

#Reef figure


fig_compare_reefs <- ggplot(data=reefs_compare) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=method),
              alpha=0.8) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, color=method), linewidth=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, size=H_significance, fill=method), alpha=0.85, shape=21) +
  annotate("text", x=500, y=3, label="Reef-builders", size=3) +
  scale_color_manual(values=c('black', 'black')) +
  scale_fill_manual(values=c('#3d405b', '#f2cc8f')) +
  scale_size_manual(values=c(1, 1.5, 2, 2.75)) +
  scale_alpha_manual(values=c(0.1, 0.3, 0.7, 1.0)) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill=guide_legend(byrow=TRUE), size='none') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme
fig_compare_reefs

compare_subsampling <- ggarrange(fig_compare_bioturbators, fig_compare_reefs, ncol=1)
compare_subsampling



#=== effect size and comparison of H results (Supplementary Figures) ====#
#restructure M1 vs M2 data
m1_H <- bioturbators_spatial_results_df[,c(1:3,12:13)]
colnames(m1_H) <- c('period', 'stage', 'mid_ma', 'H', 'sd')
m1_H$formations <- 'Bioturbators present'
m2_H <- bioturbators_spatial_results_df[,c(1:3,14:15)]
colnames(m2_H) <- c('period', 'stage', 'mid_ma', 'H', 'sd')
m2_H$formations <- 'Bioturbators absent'
compare_H <- rbind(m1_H, m2_H)

bioturbators_spatial_results_df$stage <- factor(bioturbators_spatial_results_df$stage, levels=stage_names)
bioturbators_spatial_results_df$period <- factor(bioturbators_spatial_results_df$period, levels=period_names)

compare.forms.cols <- c('#2759A9', '#8F9AAD')
compare_H$formations <- factor(compare_H$formations, 
                               levels=c('Bioturbators present', 'Bioturbators absent'))
compare_H_plot <- ggplot(data=compare_H) +
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

effectsize_H_plot <- ggplot(data=bioturbators_spatial_results_df) +
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

bioturbation_H_cells <- ggarrange2(compare_H_plot, effectsize_H_plot, ncol=1)

