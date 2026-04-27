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
  legend.position.inside=c(0.78, 0.8),
  legend.title=element_blank(),
  legend.text=element_text(size=6),
  legend.background=element_rect(fill=NA, colour=NA),
  plot.title=element_text(hjust=5),
  axis.text = element_text(color = "black", size=8),
  axis.title = element_text(size=10),
  axis.line.x = element_blank())

#==== compare spatial-level subsampling for bioturbators ====#
#LOAD DATA HERE 

#assign effect size significance thresholds
bioturbation_occsub_results_df$H_lowerbound <- abs(bioturbation_occsub_results_df$HedgesG_H) - bioturbation_occsub_results_df$g_H_sd
bioturbation_occsub_results_df$H_significance <- rep(NA, nrow(bioturbation_occsub_results_df))
bioturbation_occsub_results_df[which(abs(bioturbation_occsub_results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
bioturbation_occsub_results_df[which(abs(bioturbation_occsub_results_df$H_lowerbound) < 0.2 & is.na(bioturbation_occsub_results_df$H_significance)),'H_significance'] <- 'no effect'
bioturbation_occsub_results_df[which( abs(bioturbation_occsub_results_df$H_lowerbound) >= 0.2 & abs(bioturbation_occsub_results_df$H_lowerbound) < 0.5 & is.na(bioturbation_occsub_results_df$H_significance)),'H_significance'] <- 'weak effect'
bioturbation_occsub_results_df[which( abs(bioturbation_occsub_results_df$H_lowerbound) >= 0.5 & abs(bioturbation_occsub_results_df$H_lowerbound) < 0.8 & is.na(bioturbation_occsub_results_df$H_significance)),'H_significance'] <- 'moderate effect'
bioturbation_occsub_results_df[which( abs(bioturbation_occsub_results_df$H_lowerbound) >= 0.8 & is.na(bioturbation_occsub_results_df$H_significance)),'H_significance'] <- 'strong effect'
bioturbation_occsub_results_df$H_significance <- factor(bioturbation_occsub_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
bioturbation_occsub_results_df$method <- c('formations - no threshold')

bioturbation_threshold_results_df$H_lowerbound <- abs(bioturbation_threshold_results_df$HedgesG_H) - bioturbation_threshold_results_df$g_H_sd
bioturbation_threshold_results_df$H_significance <- rep(NA, nrow(bioturbation_threshold_results_df))
bioturbation_threshold_results_df[which(abs(bioturbation_threshold_results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
bioturbation_threshold_results_df[which(abs(bioturbation_threshold_results_df$H_lowerbound) < 0.2 & is.na(bioturbation_threshold_results_df$H_significance)),'H_significance'] <- 'no effect'
bioturbation_threshold_results_df[which( abs(bioturbation_threshold_results_df$H_lowerbound) >= 0.2 & abs(bioturbation_threshold_results_df$H_lowerbound) < 0.5 & is.na(bioturbation_threshold_results_df$H_significance)),'H_significance'] <- 'weak effect'
bioturbation_threshold_results_df[which( abs(bioturbation_threshold_results_df$H_lowerbound) >= 0.5 & abs(bioturbation_threshold_results_df$H_lowerbound) < 0.8 & is.na(bioturbation_threshold_results_df$H_significance)),'H_significance'] <- 'moderate effect'
bioturbation_threshold_results_df[which( abs(bioturbation_threshold_results_df$H_lowerbound) >= 0.8 & is.na(bioturbation_threshold_results_df$H_significance)),'H_significance'] <- 'strong effect'
bioturbation_threshold_results_df$H_significance <- factor(bioturbation_threshold_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
bioturbation_threshold_results_df$method <- c('formations - threshold')

bioturbators_compare <- rbind(bioturbation_occsub_results_df, bioturbation_threshold_results_df)


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


#LOAD DATA HERE 



#assign effect size significance thresholds
bioturbators_spatial_results_df$H_lowerbound <- abs(bioturbators_spatial_results_df$HedgesG_H) - bioturbators_spatial_results_df$g_H_sd
bioturbators_spatial_results_df$H_significance <- rep(NA, nrow(bioturbators_spatial_results_df))
bioturbators_spatial_results_df[which(abs(bioturbators_spatial_results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
bioturbators_spatial_results_df[which(abs(bioturbators_spatial_results_df$H_lowerbound) < 0.2 & is.na(bioturbators_spatial_results_df$H_significance)),'H_significance'] <- 'no effect'
bioturbators_spatial_results_df[which( abs(bioturbators_spatial_results_df$H_lowerbound) >= 0.2 & abs(bioturbators_spatial_results_df$H_lowerbound) < 0.5 & is.na(bioturbators_spatial_results_df$H_significance)),'H_significance'] <- 'weak effect'
bioturbators_spatial_results_df[which( abs(bioturbators_spatial_results_df$H_lowerbound) >= 0.5 & abs(bioturbators_spatial_results_df$H_lowerbound) < 0.8 & is.na(bioturbators_spatial_results_df$H_significance)),'H_significance'] <- 'moderate effect'
bioturbators_spatial_results_df[which( abs(bioturbators_spatial_results_df$H_lowerbound) >= 0.8 & is.na(bioturbators_spatial_results_df$H_significance)),'H_significance'] <- 'strong effect'
bioturbators_spatial_results_df$H_significance <- factor(bioturbators_spatial_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
bioturbators_spatial_results_df$method <- c('equal area cells - no threshold')

bioturbators_spatial_threshold_results_df$H_lowerbound <- abs(bioturbators_spatial_threshold_results_df$HedgesG_H) - bioturbators_spatial_threshold_results_df$g_H_sd
bioturbators_spatial_threshold_results_df$H_significance <- rep(NA, nrow(bioturbators_spatial_threshold_results_df))
bioturbators_spatial_threshold_results_df[which(abs(bioturbators_spatial_threshold_results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
bioturbators_spatial_threshold_results_df[which(abs(bioturbators_spatial_threshold_results_df$H_lowerbound) < 0.2 & is.na(bioturbators_spatial_threshold_results_df$H_significance)),'H_significance'] <- 'no effect'
bioturbators_spatial_threshold_results_df[which( abs(bioturbators_spatial_threshold_results_df$H_lowerbound) >= 0.2 & abs(bioturbators_spatial_threshold_results_df$H_lowerbound) < 0.5 & is.na(bioturbators_spatial_threshold_results_df$H_significance)),'H_significance'] <- 'weak effect'
bioturbators_spatial_threshold_results_df[which( abs(bioturbators_spatial_threshold_results_df$H_lowerbound) >= 0.5 & abs(bioturbators_spatial_threshold_results_df$H_lowerbound) < 0.8 & is.na(bioturbators_spatial_threshold_results_df$H_significance)),'H_significance'] <- 'moderate effect'
bioturbators_spatial_threshold_results_df[which( abs(bioturbators_spatial_threshold_results_df$H_lowerbound) >= 0.8 & is.na(bioturbators_spatial_threshold_results_df$H_significance)),'H_significance'] <- 'strong effect'
bioturbators_spatial_threshold_results_df$H_significance <- factor(bioturbators_spatial_threshold_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
bioturbators_spatial_threshold_results_df$method <- c('equal area cells - threshold')

bioturbators_spatial_compare <- rbind(bioturbators_spatial_results_df, bioturbators_spatial_threshold_results_df)


fig_compare_spatial_bioturbators <- ggplot(data=bioturbators_spatial_compare) +
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
fig_compare_spatial_bioturbators

#===== compare threshold test for reef-builders =====#
#formation-level
#LOAD DATA HERE 

#assign effect size significance thresholds
reefs_occsub_results_df$H_lowerbound <- abs(reefs_occsub_results_df$HedgesG_H) - reefs_occsub_results_df$g_H_sd
reefs_occsub_results_df$H_significance <- rep(NA, nrow(reefs_occsub_results_df))
reefs_occsub_results_df[which(abs(reefs_occsub_results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
reefs_occsub_results_df[which(abs(reefs_occsub_results_df$H_lowerbound) < 0.2 & is.na(reefs_occsub_results_df$H_significance)),'H_significance'] <- 'no effect'
reefs_occsub_results_df[which( abs(reefs_occsub_results_df$H_lowerbound) >= 0.2 & abs(reefs_occsub_results_df$H_lowerbound) < 0.5 & is.na(reefs_occsub_results_df$H_significance)),'H_significance'] <- 'weak effect'
reefs_occsub_results_df[which( abs(reefs_occsub_results_df$H_lowerbound) >= 0.5 & abs(reefs_occsub_results_df$H_lowerbound) < 0.8 & is.na(reefs_occsub_results_df$H_significance)),'H_significance'] <- 'moderate effect'
reefs_occsub_results_df[which( abs(reefs_occsub_results_df$H_lowerbound) >= 0.8 & is.na(reefs_occsub_results_df$H_significance)),'H_significance'] <- 'strong effect'
reefs_occsub_results_df$H_significance <- factor(reefs_occsub_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
reefs_occsub_results_df$method <- c('formations - no threshold')

reef_results_thresholds_df$H_lowerbound <- abs(reef_results_thresholds_df$HedgesG_H) - reef_results_thresholds_df$g_H_sd
reef_results_thresholds_df$H_significance <- rep(NA, nrow(reef_results_thresholds_df))
reef_results_thresholds_df[which(abs(reef_results_thresholds_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
reef_results_thresholds_df[which(abs(reef_results_thresholds_df$H_lowerbound) < 0.2 & is.na(reef_results_thresholds_df$H_significance)),'H_significance'] <- 'no effect'
reef_results_thresholds_df[which( abs(reef_results_thresholds_df$H_lowerbound) >= 0.2 & abs(reef_results_thresholds_df$H_lowerbound) < 0.5 & is.na(reef_results_thresholds_df$H_significance)),'H_significance'] <- 'weak effect'
reef_results_thresholds_df[which( abs(reef_results_thresholds_df$H_lowerbound) >= 0.5 & abs(reef_results_thresholds_df$H_lowerbound) < 0.8 & is.na(reef_results_thresholds_df$H_significance)),'H_significance'] <- 'moderate effect'
reef_results_thresholds_df[which( abs(reef_results_thresholds_df$H_lowerbound) >= 0.8 & is.na(reef_results_thresholds_df$H_significance)),'H_significance'] <- 'strong effect'
reef_results_thresholds_df$H_significance <- factor(reef_results_thresholds_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
reef_results_thresholds_df$method <- c('formations - threshold')

reefs_compare <- rbind(reefs_occsub_results_df, reef_results_thresholds_df)


fig_compare_reefs <- ggplot(data=reefs_compare) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=method),
              alpha=0.8) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, color=method), linewidth=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, size=H_significance, fill=method), alpha=0.85, shape=21) +
  scale_color_manual(values=c('black', 'black')) +
  scale_fill_manual(values=c('#86BBD8', '#EFA7A7')) +
  scale_size_manual(values=c(1, 1.5, 2, 2.75)) +
  scale_alpha_manual(values=c(0.1, 0.3, 0.7, 1.0)) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  annotate("text", x=500, y=3, label="Reef-builders", size=3) +
  guides(fill=guide_legend(byrow=TRUE), size='none') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme
fig_compare_reefs




#cell-level
#LOAD DATA HERE 

#assign effect size significance thresholds
reefbuilders_spatial_results_df$H_lowerbound <- abs(reefbuilders_spatial_results_df$HedgesG_H) - reefbuilders_spatial_results_df$g_H_sd
reefbuilders_spatial_results_df$H_significance <- rep(NA, nrow(reefbuilders_spatial_results_df))
reefbuilders_spatial_results_df[which(abs(reefbuilders_spatial_results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
reefbuilders_spatial_results_df[which(abs(reefbuilders_spatial_results_df$H_lowerbound) < 0.2 & is.na(reefbuilders_spatial_results_df$H_significance)),'H_significance'] <- 'no effect'
reefbuilders_spatial_results_df[which( abs(reefbuilders_spatial_results_df$H_lowerbound) >= 0.2 & abs(reefbuilders_spatial_results_df$H_lowerbound) < 0.5 & is.na(reefbuilders_spatial_results_df$H_significance)),'H_significance'] <- 'weak effect'
reefbuilders_spatial_results_df[which( abs(reefbuilders_spatial_results_df$H_lowerbound) >= 0.5 & abs(reefbuilders_spatial_results_df$H_lowerbound) < 0.8 & is.na(reefbuilders_spatial_results_df$H_significance)),'H_significance'] <- 'moderate effect'
reefbuilders_spatial_results_df[which( abs(reefbuilders_spatial_results_df$H_lowerbound) >= 0.8 & is.na(reefbuilders_spatial_results_df$H_significance)),'H_significance'] <- 'strong effect'
reefbuilders_spatial_results_df$H_significance <- factor(reefbuilders_spatial_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
reefbuilders_spatial_results_df$method <- c('equal area cells - no threshold')

reefbuilders_spatial_threshold_results_df$H_lowerbound <- abs(reefbuilders_spatial_threshold_results_df$HedgesG_H) - reefbuilders_spatial_threshold_results_df$g_H_sd
reefbuilders_spatial_threshold_results_df$H_significance <- rep(NA, nrow(reefbuilders_spatial_threshold_results_df))
reefbuilders_spatial_threshold_results_df[which(abs(reefbuilders_spatial_threshold_results_df$HedgesG_H) < 0.2),'H_significance'] <- 'no effect'
reefbuilders_spatial_threshold_results_df[which(abs(reefbuilders_spatial_threshold_results_df$H_lowerbound) < 0.2 & is.na(reefbuilders_spatial_threshold_results_df$H_significance)),'H_significance'] <- 'no effect'
reefbuilders_spatial_threshold_results_df[which( abs(reefbuilders_spatial_threshold_results_df$H_lowerbound) >= 0.2 & abs(reefbuilders_spatial_threshold_results_df$H_lowerbound) < 0.5 & is.na(reefbuilders_spatial_threshold_results_df$H_significance)),'H_significance'] <- 'weak effect'
reefbuilders_spatial_threshold_results_df[which( abs(reefbuilders_spatial_threshold_results_df$H_lowerbound) >= 0.5 & abs(reefbuilders_spatial_threshold_results_df$H_lowerbound) < 0.8 & is.na(reefbuilders_spatial_threshold_results_df$H_significance)),'H_significance'] <- 'moderate effect'
reefbuilders_spatial_threshold_results_df[which( abs(reefbuilders_spatial_threshold_results_df$H_lowerbound) >= 0.8 & is.na(reefbuilders_spatial_threshold_results_df$H_significance)),'H_significance'] <- 'strong effect'
reefbuilders_spatial_threshold_results_df$H_significance <- factor(reefbuilders_spatial_threshold_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))
reefbuilders_spatial_threshold_results_df$method <- c('equal area cells - threshold')

reefs_spatial_compare <- rbind(reefbuilders_spatial_results_df, reefbuilders_spatial_threshold_results_df)


fig_compare_spatial_reefs <- ggplot(data=reefs_spatial_compare) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=method),
              alpha=0.8) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, color=method), linewidth=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, size=H_significance, fill=method), alpha=0.85, shape=21) +
  scale_color_manual(values=c('black', 'black')) +
  scale_fill_manual(values=c('#86BBD8', '#EFA7A7')) +
  scale_size_manual(values=c(1, 1.5, 2, 2.75)) +
  scale_alpha_manual(values=c(0.1, 0.3, 0.7, 1.0)) +
  scale_x_reverse(limits=c(538,-5), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  annotate("text", x=500, y=3, label="Reef-builders", size=3) +
  guides(fill=guide_legend(byrow=TRUE), size='none') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme
fig_compare_spatial_reefs


panel_fig <- ggarrange(fig_compare_bioturbators, fig_compare_spatial_bioturbators,
          fig_compare_reefs, fig_compare_spatial_reefs, ncol=2, nrow=2)
