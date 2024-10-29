#Authors: ALison Cribb, Will Gearty
#Summary: temperature correlates with effect sizes


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
library(wesanderson)
library(fANCOVA)
wes_tempcols <- wes_palette('Zissou1', 100, type='continuous')

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
load('Output/effectsizes_reefs_occsub.RData')
reef_results_df <- reefs_occsub_results_df
reef_results_df$stage <- factor(reef_results_df$stage, levels=stage_names)
reef_results_df$period <- factor(reef_results_df$period, levels=period_names)

reef_results_df$H_lowerbound <- abs(reef_results_df$HedgesG_H) - reef_results_df$g_H_sd
reef_results_df$H_significance <- rep(NA, nrow(reef_results_df))
reef_results_df[which(abs(reef_results_df$H_lowerbound) < 0.2),'H_significance'] <- 'no effect'
reef_results_df[which( abs(reef_results_df$H_lowerbound) >= 0.2 & abs(reef_results_df$H_lowerbound) < 0.5),'H_significance'] <- 'weak effect'
reef_results_df[which( abs(reef_results_df$H_lowerbound) >= 0.5 & abs(reef_results_df$H_lowerbound) < 0.8),'H_significance'] <- 'moderate effect'
reef_results_df[which( abs(reef_results_df$H_lowerbound) >= 0.8),'H_significance'] <- 'strong effect'
reef_results_df$H_significance <- factor(reef_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))


load('Output/effectsizes_bioturbation_occsub.RData')
bio_results_df <- bioturbation_occsub_results_df
bio_results_df$stage <- factor(bio_results_df$stage, levels=stage_names)
bio_results_df$period <- factor(bio_results_df$period, levels=period_names)

bio_results_df$H_lowerbound <- abs(bio_results_df$HedgesG_H) - bio_results_df$g_H_sd
bio_results_df$H_significance <- rep(NA, nrow(bio_results_df))
bio_results_df[which(abs(bio_results_df$H_lowerbound) < 0.2),'H_significance'] <- 'no effect'
bio_results_df[which( abs(bio_results_df$H_lowerbound) >= 0.2 & abs(bio_results_df$H_lowerbound) < 0.5),'H_significance'] <- 'weak effect'
bio_results_df[which( abs(bio_results_df$H_lowerbound) >= 0.5 & abs(bio_results_df$H_lowerbound) < 0.8),'H_significance'] <- 'moderate effect'
bio_results_df[which( abs(bio_results_df$H_lowerbound) >= 0.8),'H_significance'] <- 'strong effect'
bio_results_df$H_significance <- factor(bio_results_df$H_significance, levels=c('no effect', 'weak effect', 'moderate effect' , 'strong effect'))

GATs <- as.data.frame(read.csv('Data/scotese_climate/Scotese-2021_GATs1Myr.csv'))
stage_bottoms <- stages$bottom[4:95]
stage_tops <- stages$top[4:95]
stage_durs <- stages$dur[4:95]

variables <- c('stage', 'age_bottom', 'age_mid', 'age_top', 'dur', 'mean_temp')
stage.temps <- as.data.frame(matrix(NA, nrow=length(stage_names), ncol=length(variables)))
colnames(stage.temps) <- variables
rownames(stage.temps) <- stage_names
stage.temps$stage <- stage_names
stage.temps$age_bottom <- stage_bottoms 
stage.temps$age_mid <- stage_mids
stage.temps$age_top <- stage_tops 
stage.temps$dur <- stage_durs
stage.temps

for(i in 1:nrow(stage.temps)){
  
  this.stage <- rownames(stage.temps)[i]
  this.stage.bottom <- stage.temps[this.stage, 'age_bottom']
  this.stage.top <- stage.temps[this.stage, 'age_top']
  this.stage.temps <- subset(GATs, as.numeric(Age) %in% round(c(this.stage.top:this.stage.bottom)))
  avg.stage.temp <- mean(this.stage.temps$GAT)
  stage.temps[this.stage,'mean_temp'] <- avg.stage.temp
  
}

stage.temps
stage.temps$climate <- rep(NA, nrow(stage.temps))
stage.temps[which(stage.temps$mean_temp<=18),'climate'] <- 'icecaps'
stage.temps[which(stage.temps$mean_temp>18),'climate'] <- 'no icecaps'

stage.temps$d_temp <- rep(NA, nrow(stage.temps))
for(i in 2:nrow(stage.temps)){
  this.stage.temp <- stage.temps$mean_temp[i]
  previous.stage.temp <- stage.temps$mean_temp[i-1]
  stage.temps$d_temp[i] <- this.stage.temp-previous.stage.temp
}
stage.temps


#==== effect size figure with climate bands ====#

climatebands_effects_reefs <- ggplot(data=subset(reef_results_df, !is.na(HedgesG_genrich))) +
  #climate bands
  geom_rect(data=stage.temps, aes(xmin = age_bottom, xmax = age_top, ymin = -1.2, ymax = 2.5, fill = mean_temp, col=mean_temp), size = 0.5) +
  #effect sizes 
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=8) +
  geom_line(aes(x=mid_ma, y=HedgesG_H), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, size=H_significance), fill='#D56482', shape=21) +
  scale_size_manual(values=c(1, 2.2, 2.8, 3.4)) +
  scale_fill_gradientn(colours=wes_tempcols, name=expression(Temp~(degree*C))) +
  scale_color_gradientn(colours=wes_tempcols, name=expression(Temp~(degree*C))) +
  scale_x_reverse(limits=c(538.8,0), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,2.5), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill='none') +
  annotate("text", x=525, y=2.3, label="C", size=4) +
  coord_geo(expand=FALSE, pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme +
  theme(legend.position='none') 
climatebands_effects_reefs

climatebands_effects_bioturbators <- ggplot(data=subset(bio_results_df, !is.na(HedgesG_genrich))) +
  #climate bands
  geom_rect(data=stage.temps, aes(xmin = age_bottom, xmax = age_top, ymin = -1.2, ymax = 3.2, fill = mean_temp, col=mean_temp), size = 0.5) +
  #effect sizes 
  geom_errorbar(aes(x=mid_ma, ymin=HedgesG_H-g_H_sd,
                    ymax=HedgesG_H+g_H_sd),
                width=8) +
  geom_line(aes(x=mid_ma, y=HedgesG_H), size=0.4) +
  geom_point(aes(x=mid_ma, y=HedgesG_H, size=H_significance), fill='#1D5587', shape=21) +
  scale_size_manual(values=c(1, 2.2, 2.8, 3.4)) +
  scale_fill_gradientn(colours=wes_tempcols, name=expression(Temp~(degree*C))) +
  scale_color_gradientn(colours=wes_tempcols, name=expression(Temp~(degree*C))) +
  scale_x_reverse(limits=c(538.8,0), name='Time (mya)') +
  scale_y_continuous(limits=c(-1.2,3.2), name=expression(paste("Hedges' g (\u00b11",sigma,")"))) +
  guides(fill='none') +
  annotate("text", x=525, y=3, label="A", size=4) +
  coord_geo(expand=FALSE, pos='bottom', dat='periods', size=3, abbrv=TRUE, height=unit(1,'line')) +
  theme_classic() +
  effectsize_theme +
  theme(legend.position='none') 
climatebands_effects_bioturbators

climatebands_supp <- ggarrange(climatebands_effects_bioturbators, climatebands_effects_reefs, ncol=1)


#==== loess curves ====#
temp_correlate_reefs <- as.data.frame(cbind(stage.temps$mean_temp, reef_results_df$HedgesG_H))
colnames(temp_correlate_reefs) <- c('GAT', 'HedgesG_H')
temp_correlate_reefs <- subset(temp_correlate_reefs, !is.na(HedgesG_H))

#GATs vs effect sizes - use fANCOVA::loess.as
reefs.loess <- loess.as(x=temp_correlate_reefs$GAT, 
                        y=temp_correlate_reefs$HedgesG_H,
                        degree=1,
                        criterion=('gcv'),
                        family='gaussian')

reefs.yhat.loess <- reefs.loess$fitted
reefs.x.loess <- temp_correlate_reefs$GAT

temp_correlate_reefs$yhat.loess <- reefs.yhat.loess
temp_correlate_reefs$x.loess <- reefs.x.loess

#GAT vs effect sizes - use base R loess function
reefs.loess_base <- loess(HedgesG_H ~ GAT, data=temp_correlate_reefs, span=reefs.loess$pars$span, degree=1) #get span from fANCOVA output
predictions <- predict(reefs.loess_base, se=TRUE) #get predictions with standard error
reefs.yhat.loess_base <- predictions$fit
reefs.se.loess_base <- predictions$se.fit
reefs.upper.ci_base <- reefs.yhat.loess_base + 1.96 * reefs.se.loess_base #calcualate 95 CI assuming normal distribution - upper
reefs.lower.ci_base <- reefs.yhat.loess_base - 1.96 * reefs.se.loess_base #lower 
temp_correlate_reefs$reefs.yhat.loess_base <- reefs.yhat.loess_base
temp_correlate_reefs$reefs.se.loess_base <- reefs.se.loess_base
temp_correlate_reefs$reefs.upper.ci_base <- reefs.upper.ci_base
temp_correlate_reefs$reefs.lower.ci_base <- reefs.lower.ci_base

GAT_v_effects_reefs <- ggplot(data=temp_correlate_reefs, aes(x=GAT, y=HedgesG_H)) +
  geom_ribbon(aes(x=GAT, ymin=reefs.lower.ci_base, ymax=reefs.upper.ci_base), fill='#D56482', alpha=0.3) +
  geom_line(aes(x=GAT, y=reefs.yhat.loess_base), color='#D56482', size=1.5) +
  geom_point(shape=21, fill='#933B53', col='black', size=3) +
  scale_x_continuous('Global Average Temperature', expand=c(0,0)) +
  scale_y_continuous('Effect size of reefs on H', expand=c(0,0), limits=c(0,1.7)) +
  annotate("text", y=1.6, x=10.5, label='D') +
  theme_classic() +
  effectsize_theme
GAT_v_effects_reefs

#==== loess curves ====#
temp_correlate_bioturbators <- as.data.frame(cbind(stage.temps$mean_temp, bio_results_df$HedgesG_H))
colnames(temp_correlate_bioturbators) <- c('GAT', 'HedgesG_H')
temp_correlate_bioturbators <- subset(temp_correlate_bioturbators, !is.na(HedgesG_H))

#GATs vs effect sizes - use fANCOVA::loess.as
biots.loess <- loess.as(x=temp_correlate_bioturbators$GAT, 
                        y=temp_correlate_bioturbators$HedgesG_H,
                        degree=1,
                        criterion=('gcv'),
                        family='gaussian')

biots.yhat.loess <- biots.loess$fitted
biots.x.loess <- temp_correlate_bioturbators$GAT

temp_correlate_bioturbators$yhat.loess <- biots.yhat.loess
temp_correlate_bioturbators$x.loess <- biots.x.loess

#GAT vs effect sizes - use base R loess function
biots.loess_base <- loess(HedgesG_H ~ GAT, data=temp_correlate_bioturbators, span=biots.loess$pars$span, degree=1) #get span from fANCOVA output
biots.predictions <- predict(biots.loess_base, se=TRUE) #get predictions with standard error
biots.yhat.loess_base <- biots.predictions$fit
biots.se.loess_base <- biots.predictions$se.fit
biots.upper.ci_base <- biots.yhat.loess_base + 1.96 * biots.se.loess_base #calcualate 95 CI assuming normal distribution - upper
biots.lower.ci_base <- biots.yhat.loess_base - 1.96 * biots.se.loess_base #lower 
temp_correlate_bioturbators$biots.yhat.loess_base <- biots.yhat.loess_base
temp_correlate_bioturbators$biots.se.loess_base <- biots.se.loess_base
temp_correlate_bioturbators$biots.upper.ci_base <- biots.upper.ci_base
temp_correlate_bioturbators$biots.lower.ci_base <- biots.lower.ci_base

GAT_v_effects_bio <- ggplot(data=temp_correlate_bioturbators, aes(x=GAT, y=HedgesG_H)) +
  geom_ribbon(aes(x=GAT, ymin=biots.lower.ci_base, ymax=biots.upper.ci_base), fill='#1491FF', alpha=0.3) +
  geom_line(aes(x=GAT, y=biots.yhat.loess_base), color='#1491FF', size=1.5) +
  geom_point(shape=21, fill='#1D5587', col='black', size=3) +
  scale_x_continuous('Global Average Temperature', expand=c(0,0)) +
  scale_y_continuous('Effect size of bioturbators on H', expand=c(0,0), limits=c(0,2)) +
  annotate("text", y=1.9, x=10.5, label='B') +
  theme_classic() +
  effectsize_theme
GAT_v_effects_bio

climate_correlates <- ggarrange(climatebands_effects_bioturbators, climatebands_effects_reefs,
                                GAT_v_effects_bio, GAT_v_effects_reefs, ncol=2)
