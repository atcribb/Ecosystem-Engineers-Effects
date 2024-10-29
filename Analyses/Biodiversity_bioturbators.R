#Authors: Alison Cribb
#Summary: New ecosystem engineering taxa, diversity of ecosystem engineers

set.seed(541)

library(divDyn)

#==== Load data ====# 
load('Data/Bioturbators_data.RData') #bioturbators data 

data("stages", package="divDyn")
stage_names <- stages$stage[4:95]
stage_mids <- stages$mid[4:95]

#=== subsampling quota ===#
sampling.counts <- as.data.frame(matrix(NA, nrow=length(stage_names), ncol=4))
colnames(sampling.counts) <- c('stage', 'n_occs', 'n_colls', 'included_in_subsampling')
sampling.counts$stage <- stage_names
for(i in 1:nrow(sampling.counts)){
  
  this.stage<-sampling.counts$stage[i]
  this.stage.data <- subset(bioturbators_data, stage==this.stage)
  sampling.counts$n_occs[i] <- nrow(this.stage.data)
  sampling.counts$n_colls[i] <- length(unique(this.stage.data$collection_no))
  if(sampling.counts$n_occs[i]>150){sampling.counts$included_in_subsampling[i]<-'yes'}
  if(sampling.counts$n_occs[i]<=150){sampling.counts$included_in_subsampling[i]<-'no'}
  
}
#View(sampling.counts)

occs.quota <- 150 #set quota - we only want to use the stages where we have sufficient data
quota.stages <- subset(sampling.counts, n_occs>occs.quota)$stage
quota_data <- subset(bioturbators_data, stage %in% quota.stages)

#=== set up ====#
variables <- c('stage', 'mid_ma', 'bioturbator_richness', 'rich_sd', 'n_new_genera', 'new_sd')
ee_div_df <- as.data.frame(matrix(NA, nrow=length(stage_names), ncol=length(variables)))
colnames(ee_div_df) <- variables
ee_div_df$stage <- stage_names
ee_div_df$mid_ma <- stage_mids
#ee_div_df

iter <- 1000

bioturbator_richness <- as.data.frame(matrix(NA, nrow=iter, ncol=length(stage_names)))
colnames(bioturbator_richness) <- stage_names
#head(bioturbator_richness)

new_bioturbators <- as.data.frame(matrix(NA, nrow=iter, ncol=length(stage_names)))
colnames(new_bioturbators) <- stage_names
#head(new_bioturbators)

cat('starting analyses')

for(i in 1:iter){
  
  #create empty dataframe to rbind together subsampled data from each time bin
  subsampled.data <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(bioturbators_data)))
  colnames(subsampled.data) <- colnames(bioturbators_data)
  
  for(j in 1:length(stage_names)){
    this.stage <- stage_names[j]
    this.stage.data <- subset(bioturbators_data, stage==this.stage)
    max.occs <- nrow(this.stage.data)
    row.idxs <- sample(max.occs, occs.quota, replace=TRUE)
    this.stage.subsampled.data <- this.stage.data[row.idxs,]
    subsampled.data <- rbind(subsampled.data, this.stage.subsampled.data)
  }
  
  #range-through subsampled data
  subsampled.ranges <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(subsampled.data)))
  colnames(subsampled.ranges) <- colnames(subsampled.data)
  genus.list <- unique(subsampled.data$genus)
  for(k in 1:length(genus.list)){
    this.genus <- genus.list[k]
    all.genus.data <- subset(subsampled.data, genus==this.genus)
    FAD.genus.stg <- min(all.genus.data$stg) #oldest stage number
    LAD.genus.stg <- max(all.genus.data$stg) #youngest stage number
    stage.range <- FAD.genus.stg:LAD.genus.stg #stage number ranges 
    genus.occ.entry <- all.genus.data[1,] #get the first occurrence entry to duplicate
    genus.range.data <- rbind(genus.occ.entry,
                              genus.occ.entry[rep(1, length(stage.range)-1),]) #duplicate entries for FAD-LAD range
    genus.range.data$stg <- stage.range #assign each entry a stage number (one per stage)
    for(m in 1:nrow(genus.range.data)){ #assign stage names back to corresponding stage numbers
      stg.no <- genus.range.data$stg[m]
      genus.range.data$stage[m] <-subset(stages, stg==stg.no)$stage
    }
    
    subsampled.ranges <- rbind(subsampled.ranges, genus.range.data)
  }
  
  for(t in 1:ncol(bioturbator_richness)){
    this.stage.rt <- colnames(bioturbator_richness)[t]
    this.stage.data <- subset(subsampled.ranges, stage==this.stage.rt)
    
    #calculate generic richness
    bioturbator_richness[i,this.stage.rt] <- length(unique(this.stage.data$clgen))
    
    #calculate # new genera 
    if(t>1){ #skip the Fortunian for this 
      previous.stages <- colnames(bioturbator_richness)[1:t-1] 
      previous.stages.data <- subset(subsampled.ranges, stage %in% previous.stages)
      already.genera <- unique(previous.stages.data$genus) #list of genera that are already present in previous stages
      this.stage.new.data <- subset(this.stage.data, !(genus %in% already.genera)) #get data of genera that are NOT present in preivous stages
      new_bioturbators[i,this.stage.rt] <- length(unique(this.stage.new.data$genus)) #get number of new genera in that stage
    }
    
  }
  
  cat('finished iteration #',i,'\n')
  
}

ee_div_df$bioturbator_richness <- colMeans(bioturbator_richness)
ee_div_df$rich_sd <- apply(bioturbator_richness,2,sd)
ee_div_df$n_new_genera <- colMeans(new_bioturbators)
ee_div_df$new_sd <- apply(new_bioturbators,2,sd)
ee_div_df <- subset(ee_div_df, stage %in% quota.stages)
#View(ee_div_df)

bioturbators_div <- ee_div_df
save(bioturbators_div, file='Output/bioturbators_evolution_biodiversity.RData')
