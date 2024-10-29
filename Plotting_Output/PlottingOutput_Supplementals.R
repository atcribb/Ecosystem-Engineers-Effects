#Authors: Alison Cribb, Will Gearty
#Summary: Scripts that will produce supplementary figures 1, 3, and 4

#==== Supplemental Figure 1 compare to results without ecosystem engineers ====#
load('Output/effectsizes_bioturbation_occsub.RData')
results_df <-  bioturbation_occsub_results_df
load('Output/effectsizes_bioturbators_noEEs.RData')
biots_compare_remove <- as.data.frame(cbind(bioturbation_noEEs_results_df$stage,
                                      bioturbation_noEEs_results_df$mid_ma,
                                      bioturbation_noEEs_results_df$HedgesG_H,
                                      bioturbation_noEEs_results_df$g_H_sd,
                                      rep('bioturbators removed', nrow(bioturbation_noEEs_results_df))))
colnames(biots_compare_remove) <- c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd', 'treatment')
biots_compare_keeps <- as.data.frame(cbind(results_df$stage,
                                     results_df$mid_ma,
                                     results_df$HedgesG_H,
                                     results_df$g_H_sd,
                                     rep('bioturbators kept', nrow(results_df))))
colnames(biots_compare_keeps) <- colnames(biots_compare_remove)
biots_compare_results <- as.data.frame(rbind(biots_compare_remove, biots_compare_keeps))
biots_compare_results$mid_ma <- as.numeric(biots_compare_results$mid_ma)
biots_compare_results$HedgesG_H <- as.numeric(biots_compare_results$HedgesG_H)
biots_compare_results$g_H_sd <- as.numeric(biots_compare_results$g_H_sd)
biots_compare_results

biots_compare_H <- ggplot(data=biots_compare_results) +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=treatment), alpha=0.5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=treatment)) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=treatment), shape=21, colour='black', size=3, alpha=0.6) +
  scale_fill_manual(values=c('#9E5FA9', '#E49C00')) +
  scale_colour_manual(values=c('#9E5FA9', '#E49C00')) +
  scale_y_continuous(limits=c(0, 3), name="Hedges' g") +
  scale_x_reverse(limits=c(538,-5), 'Time (mya)') +
  annotate("text", x=490, y=2.9, label="A) Bioturbators") +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  theme(
        legend.position=c(0.85,0.9),
        legend.title=element_blank())
biots_compare_H

load('Output/effect_sizes/effectsizes_reefs_occsub.RData')
load('Output/effectsizes_reefs_noEEs.RData')
reefs_compare_remove <- as.data.frame(cbind(reef_noEEs_results_df$stage,
                                            reef_noEEs_results_df$mid_ma,
                                            reef_noEEs_results_df$HedgesG_H,
                                            reef_noEEs_results_df$g_H_sd,
                                            rep('reef-builders removed', nrow(reef_noEEs_results_df))))
colnames(reefs_compare_remove) <- c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd', 'treatment')
reefs_compare_keeps <- as.data.frame(cbind(reefs_occsub_results_df$stage,
                                            reefs_occsub_results_df$mid_ma,
                                            reefs_occsub_results_df$HedgesG_H,
                                            reefs_occsub_results_df$g_H_sd,
                                           rep('reef-builders kept', nrow(reefs_occsub_results_df))))
colnames(reefs_compare_keeps) <- colnames(reefs_compare_remove)
reefs_compare_results <- as.data.frame(rbind(reefs_compare_remove, reefs_compare_keeps))
reefs_compare_results$mid_ma <- as.numeric(reefs_compare_results$mid_ma)
reefs_compare_results$HedgesG_H <- as.numeric(reefs_compare_results$HedgesG_H)
reefs_compare_results$g_H_sd <- as.numeric(reefs_compare_results$g_H_sd)
reefs_compare_results

reefs_compare_H <- ggplot(data=reefs_compare_results) +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=treatment), alpha=0.5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=treatment)) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=treatment), shape=21, colour='black', size=3, alpha=0.6) +
  scale_fill_manual(values=c('#2a9d8f', '#eaac8b')) +
  scale_colour_manual(values=c('#2a9d8f', '#eaac8b')) +
  scale_y_continuous(limits=c(0, 3), name="Hedges' g") +
  scale_x_reverse(limits=c(538,-5), 'Time (mya)') +
  annotate("text", x=490, y=2.9, label="B) Reef-builders") +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  theme(
        legend.position=c(0.85,0.9),
        legend.title=element_blank())
reefs_compare_H

Fig_S1 <- ggarrange(biots_compare_H, reefs_compare_H, ncol=1)


#===== Supplemental Figure 3 - Effect sizes of bioturbators within same palaeoenvironment ====#
load('Output/supplemental_effectsize_bioturbation_carbonates.RData')
load('Output/supplemental_effectsize_bioturbation_siliciclastics.RData')
load('Outputsupplemental_effectsize_bioturbation_shallow.RData')
load('Output/supplemental_effectsize_bioturbation_deep.RData')

#rework data
biots_carbonate_compare <- cbind(bioturbation_carbonate_results_df[,c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd')],
                           rep('carbonates',nrow(bioturbation_carbonate_results_df)))
colnames(biots_carbonate_compare) <- c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd', 'environment')
biots_siliciclastic_compare <-  cbind(bioturbation_siliciclastic_results_df[,c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd')],
                               rep('siliciclastics',nrow(bioturbation_siliciclastic_results_df)))
colnames(biots_siliciclastic_compare) <- colnames(biots_carbonate_compare)
biots_shallow_compare <- cbind(bioturbation_shallow_results_df[,c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd')],
                         rep('shallow',nrow(bioturbation_shallow_results_df)))
colnames(biots_shallow_compare) <- colnames(biots_carbonate_compare)
biots_deep_compare <- cbind(bioturbation_deep_results_df[,c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd')],
                         rep('deep',nrow(bioturbation_deep_results_df)))
colnames(biots_deep_compare) <- colnames(biots_carbonate_compare)
biots_all_compare <- cbind(bioturbation_occsub_results_df[,c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd')],
                           rep('all',nrow(bioturbation_occsub_results_df)))
colnames(biots_all_compare) <- colnames(biots_carbonate_compare)
biots_compare_env <- rbind(biots_carbonate_compare, biots_siliciclastic_compare, biots_shallow_compare, biots_deep_compare, biots_all_compare)
biots_compare_env$mid_ma <- as.numeric(biots_compare_env$mid_ma)
biots_compare_env$HedgesG_H <- as.numeric(biots_compare_env$HedgesG_H)
biots_compare_env$g_H_sd <- as.numeric(biots_compare_env$g_H_sd)

FigS3_a <- ggplot(subset(biots_compare_env, environment %in% c('all', 'siliciclastics'))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=environment), alpha=0.5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=environment), size=0.5) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=environment), shape=21, colour='black', size=2) +
  scale_fill_manual(values=c('#264653', '#2a9d8f')) +
  scale_colour_manual(values=c('#264653', '#2a9d8f')) +
  scale_y_continuous(limits=c(-1.2, 3.2), name="Hedges' g") +
  scale_x_reverse(limits=c(538,-5), 'Time (mya)') +
  annotate("text", x=450, y=2.9, label='A) siliciclastics') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  theme(
    legend.position=c(0.85,0.9),
    legend.title=element_blank())

FigS3_b <- ggplot(subset(biots_compare_env, environment %in% c('all', 'carbonates'))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=environment), alpha=0.5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=environment), size=0.5) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=environment), shape=21, colour='black', size=2) +
  scale_fill_manual(values=c('#264653', '#e9c46a')) +
  scale_colour_manual(values=c('#264653', '#e9c46a')) +
  scale_y_continuous(limits=c(-1.2, 3.2), name="Hedges' g") +
  scale_x_reverse(limits=c(538,-5), 'Time (mya)') +
  annotate("text", x=450, y=2.9, label='B) carbonates') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  theme(
    legend.position=c(0.85,0.9),
    legend.title=element_blank())

FigS3_c <- ggplot(subset(biots_compare_env, environment %in% c('all', 'shallow'))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=environment), alpha=0.5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=environment), size=0.5) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=environment), shape=21, colour='black', size=2) +
  scale_fill_manual(values=c('#264653', '#e76f51')) +
  scale_colour_manual(values=c('#264653', '#e76f51')) +
  scale_y_continuous(limits=c(-1.2, 3.2), name="Hedges' g") +
  scale_x_reverse(limits=c(538,-5), 'Time (mya)') +
  annotate("text", x=430, y=2.9, label='C) shallow marine') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  theme(
    legend.position=c(0.85,0.9),
    legend.title=element_blank())

FigS3_d <- ggplot(subset(biots_compare_env, environment %in% c('all', 'deep'))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=environment), alpha=0.5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=environment), size=0.5) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=environment), shape=21, colour='black', size=2) +
  scale_fill_manual(values=c('#264653', '#8338ec')) +
  scale_colour_manual(values=c('#264653', '#8338ec')) +
  scale_y_continuous(limits=c(-1.2, 3.2), name="Hedges' g") +
  scale_x_reverse(limits=c(538,-5), 'Time (mya)') +
  annotate("text", x=445, y=2.9, label='D) deep marine') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  theme(
    legend.position=c(0.85,0.9),
    legend.title=element_blank())  
FigS3 <- ggarrange(FigS3_a, FigS3_b, FigS3_c, FigS3_d, ncol=2)
ggsave('Draft/Figures/Supplementary/FigureS3_BioturbationEnvironments.pdf', FigS3,
       width=188, height=150, unit='mm')
  
#===== Supplemental Figure 4 - Effect sizes of reefs within same palaeoenvironment ====#
load('Output/supplemental_effectsize_reefs_carbonate.RData')
load('Output/supplemental_effectsize_reefs_siliciclastic.RData')
load('Output/supplemental_effectsize_reefs_shallow.RData')
load('Output/supplemental_effectsize_reefs_deep.RData')

#rework data
reefs_carbonate_compare <- cbind(reefs_carbonate_results_df[,c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd')],
                                 rep('carbonates',nrow(reefs_carbonate_results_df)))
colnames(reefs_carbonate_compare) <- c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd', 'environment')
reefs_siliciclastic_compare <-  cbind(reefs_siliciclastic_results_df[,c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd')],
                                      rep('siliciclastics',nrow(reefs_siliciclastic_results_df)))
colnames(reefs_siliciclastic_compare) <- colnames(reefs_carbonate_compare)
reefs_shallow_compare <- cbind(reefs_shallow_results_df[,c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd')],
                               rep('shallow',nrow(reefs_shallow_results_df)))
colnames(reefs_shallow_compare) <- colnames(reefs_carbonate_compare)
reefs_deep_compare <- cbind(reefs_deep_results_df[,c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd')],
                            rep('deep',nrow(reefs_deep_results_df)))
colnames(reefs_deep_compare) <- colnames(reefs_carbonate_compare)
reefs_all_compare <- cbind(reefs_occsub_results_df[,c('stage', 'mid_ma', 'HedgesG_H', 'g_H_sd')],
                           rep('all',nrow(reefs_occsub_results_df)))
colnames(reefs_all_compare) <- colnames(reefs_carbonate_compare)
reefs_compare_env <- rbind(reefs_carbonate_compare, reefs_siliciclastic_compare, reefs_shallow_compare, reefs_deep_compare, reefs_all_compare)
reefs_compare_env$mid_ma <- as.numeric(reefs_compare_env$mid_ma)
reefs_compare_env$HedgesG_H <- as.numeric(reefs_compare_env$HedgesG_H)
reefs_compare_env$g_H_sd <- as.numeric(reefs_compare_env$g_H_sd)

FigS5_a <- ggplot(subset(reefs_compare_env, environment %in% c('all', 'siliciclastics'))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=environment), alpha=0.5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=environment), size=0.5) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=environment), shape=21, colour='black', size=2) +
  scale_fill_manual(values=c('#264653', '#2a9d8f')) +
  scale_colour_manual(values=c('#264653', '#2a9d8f')) +
  scale_y_continuous(limits=c(-1.2, 3.2), name="Hedges' g") +
  scale_x_reverse(limits=c(538,-5), 'Time (mya)') +
  annotate("text", x=450, y=2.9, label='A) siliciclastics') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  theme(
    legend.position=c(0.85,0.9),
    legend.title=element_blank())

FigS5_b <- ggplot(subset(reefs_compare_env, environment %in% c('all', 'carbonates'))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=environment), alpha=0.5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=environment), size=0.5) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=environment), shape=21, colour='black', size=2) +
  scale_fill_manual(values=c('#264653', '#e9c46a')) +
  scale_colour_manual(values=c('#264653', '#e9c46a')) +
  scale_y_continuous(limits=c(-1.2, 3.2), name="Hedges' g") +
  scale_x_reverse(limits=c(538,-5), 'Time (mya)') +
  annotate("text", x=450, y=2.9, label='B) carbonates') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  theme(
    legend.position=c(0.85,0.9),
    legend.title=element_blank())

FigS5_c <- ggplot(subset(reefs_compare_env, environment %in% c('all', 'shallow'))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=environment), alpha=0.5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=environment), size=0.5) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=environment), shape=21, colour='black', size=2) +
  scale_fill_manual(values=c('#264653', '#e76f51')) +
  scale_colour_manual(values=c('#264653', '#e76f51')) +
  scale_y_continuous(limits=c(-1.2, 3.2), name="Hedges' g") +
  scale_x_reverse(limits=c(538,-5), 'Time (mya)') +
  annotate("text", x=430, y=2.9, label='C) shallow marine') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  theme(
    legend.position=c(0.85,0.9),
    legend.title=element_blank())

FigS5_d <- ggplot(subset(reefs_compare_env, environment %in% c('all', 'deep'))) +
  geom_hline(yintercept=c(-0.2,0.2), linetype='longdash', linewidth=0.5, color='gray70') +
  geom_hline(yintercept=c(-0.5,0.5), linetype='longdash', linewidth=0.5, color='gray50') +
  geom_hline(yintercept=c(-0.8,0.8), linetype='longdash', linewidth=0.5, color='gray30') +
  geom_ribbon(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd, ymax=HedgesG_H+g_H_sd, fill=environment), alpha=0.5) +
  geom_line(aes(x=mid_ma, y=HedgesG_H, colour=environment), size=0.5) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, fill=environment), shape=21, colour='black', size=2) +
  scale_fill_manual(values=c('#264653', '#8338ec')) +
  scale_colour_manual(values=c('#264653', '#8338ec')) +
  scale_y_continuous(limits=c(-1.2, 3.2), name="Hedges' g") +
  scale_x_reverse(limits=c(538,-5), 'Time (mya)') +
  annotate("text", x=445, y=2.9, label='D) deep marine') +
  coord_geo(pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  theme(
    legend.position=c(0.85,0.9),
    legend.title=element_blank())  
FigS5 <- ggarrange(FigS5_a, FigS5_b, FigS5_c, FigS5_d, ncol=2)



  
  