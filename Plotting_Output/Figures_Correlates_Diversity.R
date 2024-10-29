#Author: Alison Cribb
#Summary: Figures to plot diversity correlates

#=== load data ===#
load('Output/reefs_evolution_biodiversity.RData')
load('Output/supplemental/bioturbators_evolution_biodiversity.RData')
bioturbators_div <- ee_div_df
#load in a version of effect size data where ecosystem engineers are removed from biodiversity metrics and Hedges' g
load('Output/effectsizes_reefs_noEEs.RData')
reefs_effects <- reef_noEEs_results_df
load('Output/effectsizes_bioturbators_noEEs.RData')
bioturbators_effects <- subset(bioturbation_noEEs_results_df, stage %in% bioturbators_div$stage)

library(ggplot2)
library(deeptime)

divcorr_theme <- theme(
  panel.border=element_rect(fill=NA, colour='black'),
  axis.text=element_text(size=6),
  axis.title.x=element_text(size=7),
  axis.title.y=element_text(size=6.5)
)

#====== bioturbator figures ========#
b_divcurves <- ggplot(data=bioturbators_div) +
  geom_errorbar(aes(x=mid_ma, ymin=bioturbator_richness-rich_sd, ymax=bioturbator_richness+rich_sd), width=8) +
  geom_line(aes(x=mid_ma, y=bioturbator_richness), size=0.4, linetype='dashed', colour='gray20') +
  geom_point(aes(x=mid_ma, y=bioturbator_richness), size=2.5, shape=21, fill='#1D5587', colour='black') +
  scale_x_reverse(limits=c(538,-5), 'Time (Mya)') +
  scale_y_continuous('Global generic richness \n of bioturbators', limits=c(0,220)) +
  annotate('text', label='A', x=525, y=210, size=3) +
  coord_geo(pos='bottom', dat='periods', size=2, height=unit(1, 'line'), abbrv=TRUE) +
  theme_classic() +
  divcorr_theme

b_newgen <- ggplot(data=bioturbators_div) +
  geom_errorbar(aes(x=mid_ma, ymin=n_new_genera-new_sd, ymax=n_new_genera+new_sd), width=8) +
  geom_line(aes(x=mid_ma, y=n_new_genera), size=0.4, linetype='dashed', colour='gray20') +
  geom_point(aes(x=mid_ma, y=n_new_genera), size=2.5, shape=21, fill='#1D5587', colour='black') +
  scale_x_reverse(limits=c(538,-5), 'Time (Mya)') +
  scale_y_continuous('# genera with FAD in stage', limits=c(0,60)) +
  annotate('text', label='B', x=525, y=58, size=3) +
  coord_geo(pos='bottom', dat='periods', size=2, height=unit(1, 'line'), abbrv=TRUE) +
  theme_classic() +
  divcorr_theme

b_div_v_effects_df <- cbind(bioturbators_div$bioturbator_richness, bioturbators_effects$HedgesG_H)
colnames(b_div_v_effects_df) <- c('bioturbator_richness', 'HedgesG_H')
b_div_v_effectsH <- ggplot(b_div_v_effects_df, aes(x=bioturbator_richness, y=HedgesG_H)) +
  geom_point(size=1.5, shape=21, fill='#1D5587', colour='black') +
  geom_smooth(method='lm', colour='#1491FF', fill='#8AC8FF') +
  scale_x_continuous('Global generic richness of bioturbators') +
  scale_y_continuous('Bioturbator effect size on H') +
  annotate('text', label='C', x=0, y=2.5, size=3) +
  theme_classic() +
  divcorr_theme

b_new_v_effects_df <- cbind(bioturbators_div$n_new_genera, bioturbators_effects$HedgesG_H)
colnames(b_new_v_effects_df) <- c('n_new_genera', 'HedgesG_H')
b_new_v_effectsH <- ggplot(b_new_v_effects_df, aes(x=n_new_genera, y=HedgesG_H)) +
  geom_point(size=1.5, shape=21, fill='#1D5587', colour='black') +
  geom_smooth(method='lm', colour='#1491FF', fill='#8AC8FF') +
  scale_x_continuous('# genera with FAD in stage') +
  scale_y_continuous('Bioturbator effect size on H') +
  annotate('text', label='D', x=0, y=2.5, size=3) +
  theme_classic() +
  divcorr_theme

bioturbator_divcorrelates <- ggarrange(b_divcurves, b_newgen, b_div_v_effectsH, b_new_v_effectsH)

#====== breef figures ========#
r_divcurves <- ggplot(data=reefs_div) +
  geom_errorbar(aes(x=mid_ma, ymin=reefs_richness-rich_sd, ymax=reefs_richness+rich_sd), width=8) +
  geom_line(aes(x=mid_ma, y=reefs_richness), size=0.4) +
  geom_point(aes(x=mid_ma, y=reefs_richness), size=1.5, shape=21, fill='#933B53', colour='black') +
  scale_x_reverse(limits=c(538,-5), 'Time (Mya)') +
  scale_y_continuous('Global generic richness \n of reef-builders', limits=c(0, 125)) +
  annotate('text', label='A', x=525, y=120, size=3) +
  coord_geo(pos='bottom', dat='periods', size=2, height=unit(1, 'line'), abbrv=TRUE) +
  theme_classic() +
  divcorr_theme

r_newgen <- ggplot(data=reefs_div) +
  geom_errorbar(aes(x=mid_ma, ymin=n_new_genera-new_sd, ymax=n_new_genera+new_sd), width=8) +
  geom_line(aes(x=mid_ma, y=n_new_genera), size=0.4) +
  geom_point(aes(x=mid_ma, y=n_new_genera), size=1.5, shape=21, fill='#933B53', colour='black') +
  scale_x_reverse(limits=c(538,-5), 'Time (Mya)') +
  scale_y_continuous('# genera with FAD in stage', limits=c(0,45)) +
  annotate('text', label='B', x=525, y=43, size=3) +
  coord_geo(pos='bottom', dat='periods', size=2, height=unit(1, 'line'), abbrv=TRUE) +
  theme_classic() +
  divcorr_theme

r_div_v_effects_df <- cbind(reefs_div$reefs_richness, reefs_effects$HedgesG_genrich, reefs_effects$HedgesG_H)
colnames(r_div_v_effects_df) <- c('reefs_richness', 'HedgesG_genrich', 'HedgesG_H')
r_div_v_effectsH <- ggplot(r_div_v_effects_df, aes(x=reefs_richness, y=HedgesG_H)) +
  geom_point(size=1.5, shape=21, fill='#933B53', colour='black') +
  geom_smooth(method='lm', colour='#D56482', fill='#FFB9CC') +
  scale_x_continuous('Global generic richness of reef-builders') +
  scale_y_continuous('Reefs effect size on H') +
  annotate('text', label='C', x=0, y=2.5, size=3) +
  theme_classic() +
  divcorr_theme

r_new_v_effects_df <- cbind(reefs_div$n_new_genera, reefs_effects$HedgesG_genrich, reefs_effects$HedgesG_H)
colnames(r_new_v_effects_df) <- c('n_new_genera', 'HedgesG_genrich', 'HedgesG_H')
r_new_v_effectsH <- ggplot(r_new_v_effects_df, aes(x=n_new_genera, y=HedgesG_H)) +
  geom_point(size=1.5, shape=21, fill='#933B53', colour='black') +
  geom_smooth(method='lm', colour='#D56482', fill='#FFB9CC') +
  scale_x_continuous('# genera with FAD in stage') +
  scale_y_continuous('Reefs effect size on H') +
  annotate('text', label='D', x=0, y=2.5, size=3) +
  theme_classic() +
  divcorr_theme

reef_divcorrelates <- ggarrange(r_divcurves, r_newgen, r_div_v_effectsH, r_new_v_effectsH)



