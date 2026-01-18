#Authors: Alison Cribb, Will Gearty
#Contact: A.T.Cribb@soton.ac.uk
#Summary: Bioturbators effect size with grid-area spatial subsampling

set.seed(541)

# clear old data
rm(list = ls())

#=== USER INPUTS ===#
#how would you like to treat formation subsampling?
#subsampling.method <- 'none'        #no formation subsampling
#subsampling.method <- 'occurrences' #subsampling occurrences per formation
subsampling.method <- 'collections' #subsampling collections per formation


#==== Packages =====#
library(divDyn)
library(terra)

#===== Data input ======#
load('Data/Phanerozoic_clean_final.RData') #phanerozoic PBDB data (date accessed: 1 November 2023)
all_data <- subset(all_data, !(is.na(formation))) #remove data without formation assignments
all_data <- subset(all_data, !(formation=='')) #remove data without formation assignments
all_data <- subset(all_data, !is.na(p_lat))
all_data <- subset(all_data, !is.na(p_lng))

load('Data/Bioturbators_data.RData')
ecoeng_genera <- unique(bioturbators_data$genus) #get each ecosystem engineering genus name
ecoeng_formations <- unique(bioturbators_data$formation) #get each formation name containing ecosystem engineers

data("stages", package="divDyn") #stage info from divDyn
stage_names <- stages$stage[4:95]
stage_mids <- stages$mid[4:95]

# #===== Get paleolat and paleolng =====# 
bioturbators_data$age <- rep(NA, nrow(bioturbators_data))
all_data$age <- rep(NA, nrow(all_data))
all_reef_builders$age <- rep(NA, nrow(all_reef_builders))

for(i in 1:length(stage_names)){
  this.stage <- stage_names[i]
  this.mid <- stage_mids[i]
  all_data[which(all_data$stage==this.stage),'age'] <-  this.mid
  bioturbators_data[which(bioturbators_data$stage==this.stage), 'age'] <- this.mid
  all_reef_builders[which(all_reef_builders$stage==this.stage),'age'] <- this.mid
}

table(all_data$age)
table(bioturbators_data$age)
table(all_reef_builders$age)

all_data <- all_data %>% filter(lat <= 90) %>% filter(lat>=-90)
bioturbators_data <- bioturbators_data %>% filter(lat <= 90) %>% filter(lat>=-90)
all_reef_builders <- all_reef_builders %>% filter(lat <= 90) %>% filter(lat>=-90)

all_data <- palaeorotate(all_data, lng='lng', lat='lat', age='age', model='PALEOMAP')
bioturbators_data <- palaeorotate(bioturbators_data, lng='lng', lat='lat', age='age', model='PALEOMAP')
all_reef_builders <- palaeorotate(all_reef_builders, lng='lng', lat='lat', age='age', model='PALEOMAP')

#===== Rasterise data set up ======#
#following divvy package guidelines (e.g., Antell et al., 2024 Paleobiology)
rast.world <- rast()
proj <- 'EPSG:8857'
rast.proj <- project(rast.world, proj, res = 200000) # 200,000m is approximately 2 degrees
values(rast.proj) <- 1:ncell(rast.proj)

# #coordinate column names for current and target coordinate refernce system
xyCell <- c('cellX', 'cellY')

#get cell number and centroid coordinates associated with each occurrence in all_data
sv.all_data <- vect(all_data, geom=c('p_lng', 'p_lat'), crs='epsg:4326')
proj.all_data <- project(sv.all_data, proj)
all_data$cell <- cells(rast.proj, proj.all_data)[,'cell']
all_data[,xyCell] <- xyFromCell(rast.proj, all_data$cell)

# #==== Summary of spatial metrics =====#
spat.sum <- divvy::sdSumry(all_data, taxVar='genus', xy=xyCell, crs=proj, collections='collection_no')
variables <- c('stage', 'n_taxa', 'n_colls', 'n_sites', 'lat_range')
spatial_summary <- as.data.frame(matrix(NA, ncol=length(variables), nrow=length(stage_names)))
colnames(spatial_summary) <- variables
spatial_summary$stage <- stage_names
for(i in 1:nrow(spatial_summary)){

  this.stage <- spatial_summary$stage[i]
  this.stage.data <- subset(all_data, stage==this.stage)
  this.spat.sum <- divvy::sdSumry(this.stage.data, taxVar='genus', xy=xyCell, crs=proj, collections='collection_no')
  spatial_summary$n_taxa[i] <- as.data.frame(this.spat.sum)$nOcc
  spatial_summary$n_colls[i] <- as.data.frame(this.spat.sum)$nColl
  spatial_summary$n_sites[i] <- as.data.frame(this.spat.sum)$nLoc
  spatial_summary$lat_range[i] <- as.data.frame(this.spat.sum)$latRange
  print(paste('finished: ',this.stage))

}
View(spatial_summary)

#Summary of spatial subsampling:
# for each stage, subsampled 700 occurrences
# then break up data into equal area grid cells
# eliminate grid cells that do not have at least 20 occurrences
# randomly select 20 occurrences from the remaining grid cells 

#=== set up for analysis ===#
#set subsampling variables
n.quota <- 700 #how many fossils per stage?
occs.n.cells <- 20 #how many fossils per cell?
colls.n.cells <- 10 #how many collections per cell?
iter <- 1000 #how many iterations?

#set up dataframe for output
#genrich = generic richness
#H = Shannon's Diversity
#dom = Simpson's Dominance
variables <- c('period', 'stage', 'mid_ma', 'n_EE_forms', 'n_nonEE_forms',
               'M1_genrich', 'M1_genrich_sd', 'M2_genrich', 'M2_genrich_sd', 'HedgesG_genrich', 'g_genrich_sd',
               'M1_H', 'M1_H_sd', 'M2_H', 'M2_H_sd', 'HedgesG_H', 'g_H_sd',
               'M1_dom', 'M1_dom_sd', 'M2_dom', 'M2_dom_sd', 'HedgesG_Dominance', 'g_dom_sd')
results_df <- as.data.frame(matrix(NA, 
                                   nrow=length(stage_names),
                                   ncol=length(variables)))
colnames(results_df) <- variables

#=========================#
#====     ANALYSIS    ====#
#=========================#
for(i in 1:length(stage_names)){
  
  #time data
  this.stage <- stage_names[i] 
  this.mid <- stages[which(stages$stage==this.stage),'mid']
  this.period <- stages[which(stages$stage==this.stage),'system']
  results_df$stage[i] <- this.stage
  results_df$mid_ma[i] <- this.mid
  results_df$period[i] <- this.period
  
  #get data for each stage
  this.stage.data <- subset(all_data, stage==this.stage)
  
  #get stage equal-area cells 
  sv.this.stage.data <- vect(this.stage.data, geom=c('p_lng', 'p_lat'), crs='epsg:4326')
  proj.this.stage.data <- project(sv.this.stage.data, proj)
  this.stage.data$cell <- cells(rast.proj, proj.this.stage.data)[,'cell']

  #ecosystem engineering data
  all_EE_data <- subset(this.stage.data, formation %in% ecoeng_formations) #presence data for entire stage
  results_df$n_EE_forms[i] <- length(unique(all_EE_data$formation)) #how many formations have ecosystem engineers in this stage?
  all_nonEE_data <- subset(this.stage.data, !(formation %in% ecoeng_formations)) #absence data for entire stage
  results_df$n_nonEE_forms[i] <- length(unique(all_nonEE_data$formation)) #how many formations do not have ecosystem engineers in this stage?
  
  #temporary vectors to save subsampled data 
  M1_genrich_iters <- rep(NA, iter)
  M2_genrich_iters <- rep(NA, iter)
  HedgesG_genrich_iters <- rep(NA, iter)
  
  M1_H_iters <- rep(NA, iter)
  M2_H_iters <- rep(NA, iter)
  HedgesG_H_iters <- rep(NA, iter)
  
  M1_dom_iters <- rep(NA, iter)
  M2_dom_iters <- rep(NA, iter)
  HedgesG_dom_iters <- rep(NA, iter)
  
  #subsampling
  for(j in 1:iter){

    #subsample/bootstrap n.quota occurrences per stage 
    max.samples <- nrow(this.stage.data)
    row_idxs <- sample(max.samples, n.quota, replace=TRUE)
    subbed.data <- this.stage.data[row_idxs,]
    
    #subsample from grid cells 
    cell.nos <- unique(subbed.data$cell)
    
    spatial.subbed.data <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(subbed.data)))
    colnames(spatial.subbed.data) <- colnames(subbed.data)
    for(k in 1:length(cell.nos)){
      this.cell.no <- cell.nos[k]
      this.cell.data <- subset(subbed.data, cell==this.cell.no)
      
      if(subsampling.method=='collections'){
        this.cell.collections <- unique(this.cell.data$collection_no)
        coll.idxs <- sample(1:length(this.cell.collections), colls.n.cells, replace=TRUE)
        rand.colls <- unique(this.cell.collections[coll.idxs])
        subbed.cell.data <- subset(this.cell.data, collection_no %in% rand.colls)
      }
      
      if(subsampling.method=='occurrences'){
        max.cell.n <- nrow(this.cell.data)
        rand.rows <- unique(sample(1:max.cell.n, occ.n.cells, replace=TRUE))
        subbed.cell.data <- this.cell.data[rand.rows,]
      }
      
      spatial.subbed.data <- rbind(spatial.subbed.data, subbed.cell.data)
    }
    
    subbed.data <- spatial.subbed.data 
    
    #present ecosystem engineering data from stage-subsampled occurrences:
    presence_data <- subset(subbed.data, formation %in% ecoeng_formations)
    
    #**** If you only want to consider ecosystem engineer impacts on non-ecosystem engineering taxa, but see manuscript for discussion for why we do not do this in the final analyses, change to:
    #presence_data_all <- subset(subbed.data, formation %in% ecoeng_formations)
    #presence_data <- subset(presence_data_all, !(genus %in% ecoeng_genera))
    
    #absent ecosystem engineering data from stage-subsampled occurrences:
    absence_data <- subset(subbed.data, !(formation %in% ecoeng_formations))
    
    #get lists of presence and absence formations from that subsampled data 
    presence_formations <- unique(presence_data$formation)
    absence_formations <- unique(absence_data$formation)
    
    #=== Collect ecological statistics: ===#
    #n1,2 = sample size (no. formations)
    #x1,2 = mean generic richness/diversity/dominance 
    #s1,2 = standard deviation of means 
    
    n1 <- length(unique(presence_formations)) #number of formations containing bioturbating ecosystem engineers
    n2 <- length(unique(absence_formations)) #number of formations NOT containing bioturbating ecosystem engineers
    
    if(n1>0 & n2>0){ #don't get stuck on time periods where there are 0 ecosystem engineering formations
      #== presence statistics (x1) ==#
      genrich.presence.temp <- rep(NA, n1)
      H.presence.temp <- rep(NA, n1) #set up temporary vector to save Shannon's Diversity for each of the n1 presence formations 
      Dom.presence.temp <- rep(NA, n1) #and for Simpson's dominance 
      for(k in 1:n1){
        this.formation <- presence_formations[k] 
        this.formation.data <- subset(presence_data, formation==this.formation) #get data for each formation
        this.formation.subbed <- this.formation.data
        this.presence_genera <- unique(this.formation.subbed$genus) #get list of all genera in the subsampled formation
        this.presence_abundance_data <- as.data.frame(matrix(NA, 
                                                             nrow=length(this.presence_genera),
                                                             ncol=2)) #set up to collect how many of each genus there are 
        colnames(this.presence_abundance_data) <- c('gen', 'n')
        for(l in 1:nrow(this.presence_abundance_data)){ #for each genus
          this.genus <- this.presence_genera[l] 
          this.presence_abundance_data$gen[l] <- this.genus
          this.genus.data <- subset(this.formation.subbed, genus==this.genus) #get data for that genus in the formation
          this.presence_abundance_data$n[l] <- nrow(this.genus.data) #count how many occurrences there are 
        }
        genrich.presence.temp[k] <- nrow(this.presence_abundance_data) #calculate generic richness from n. taxa (n. rows) from abundance matrix
        
        #Shannon's diversity (H)
        presence_tot <- sum(this.presence_abundance_data$n) #total number of occurrences of all genera (formation size)
        presence_gen_props <- this.presence_abundance_data$n/presence_tot #relative abundance for each genus
        presence_shannon_div <- -sum(presence_gen_props*log(presence_gen_props)) #Shannon's Diversity formula 
        H.presence.temp[k] <- presence_shannon_div #and save 
        
        #and we can use this all to calculate Simpson's dominance
        D.presence <- (sum(presence_gen_props^2)) #sum of squared generic proportions
        Dom.presence.temp[k] <- 1/D.presence #Simpson's dominance=1/D
        
      }
      
      #Generic richness presence statistics
      x1 <- mean(genrich.presence.temp, na.rm=TRUE) #find mean generic richness across all of the presence formations
      s1 <- sd(genrich.presence.temp, na.rm=TRUE) #and the standard deviation
      #save
      M1_genrich_iters[j] <- x1
      
      #Shannon's Diversity presence statistics 
      x1_H <- mean(H.presence.temp, na.rm=TRUE) #mean Shannon's Diversity across all presence formations
      s1_H <- sd(H.presence.temp, na.rm=TRUE) #and the standard deviation
      M1_H_iters[j] <- x1_H #save
      
      #Simpson's Dominance statistics 
      x1_dom <- mean(Dom.presence.temp, na.rm=TRUE)
      s1_dom <- sd(Dom.presence.temp, na.rm=TRUE)
      M1_dom_iters[j] <- x1_dom #save 
      
      
      #== absence statistics (x2) ==#
      genrich.absence.temp <- rep(NA, n2)
      H.absence.temp <- rep(NA, n2) #set up temporary vector to save Shannon's Diversity for each of the n2 absence formations 
      Dom.absence.temp <- rep(NA, n2) #and for Simpson's dominance 
      for(k in 1:n2){
        this.formation <- absence_formations[k]
        this.formation.data <- subset(absence_data, formation==this.formation) #get data for each formation
        this.formation.subbed <- this.formation.data

        this.absence_genera <- unique(this.formation.subbed$genus) #get list of all of the genera in this absence formation
        this.absence_abundance_data <- as.data.frame(matrix(NA, 
                                                            nrow=length(this.absence_genera),
                                                            ncol=2)) #set up to collection how many occurrences of each genera in the formation
        colnames(this.absence_abundance_data) <- c('gen', 'n') 
        for(l in 1:nrow(this.absence_abundance_data)){ #for each genus
          this.genus <- this.absence_genera[l] 
          this.absence_abundance_data$gen[l] <- this.genus
          this.genus.data <- subset(this.formation.subbed, genus==this.genus) #get data for that genus in the formation
          this.absence_abundance_data$n[l] <- nrow(this.genus.data)  #count how many occurrences there are 
        }
        genrich.absence.temp[k] <- nrow(this.absence_abundance_data) #calculate generic richness from n. taxa (n. rows) from abundance matrix 
        
        #Shannon's diversity 
        absence_tot <- sum(this.absence_abundance_data$n) #total number of occurrences of all genera (formation size)
        absence_gen_props <- this.absence_abundance_data$n/absence_tot #relative abundance for each genus
        absence_shannon_div <- -sum(absence_gen_props*log(absence_gen_props)) #Shannon's Diversity formula
        H.absence.temp[k] <- absence_shannon_div
        
        #and we can sue this all to calculate Simpson's dominance
        D.absence <- (sum(absence_gen_props^2)) #sum of squared generic proportions
        Dom.absence.temp[k] <- 1/D.absence #Simpson's dominance=1/D
      }
      
      
      #Generic richness absence statistics
      x2 <- mean(genrich.absence.temp, na.rm=TRUE) #find mean generic richness across all of the absence formations
      s2 <- sd(genrich.absence.temp, na.rm=TRUE) #and the standard deviation 
      M2_genrich_iters[j] <- x2 #save
      
      #Shannon's Diversity absence statistics
      x2_H <- mean(H.absence.temp, na.rm=TRUE) #mean Shannon's Diversity across all absence formations
      s2_H <- sd(H.absence.temp, na.rm=TRUE) #and the standard deviation
      M2_H_iters[j] <- x2_H #save
      
      #Simpson's Dominance absence statistics 
      x2_dom <- mean(Dom.absence.temp, na.rm=TRUE)
      s2_dom <- sd(Dom.absence.temp, na.rm=TRUE)
      M2_dom_iters[j] <- x2_dom #save 
      
      #== Effect sizes ==#
      #Effect size for generic richness
      genrich_g <- (x1-x2)/( sqrt( ( ((n1-1)*(s1^2)) + ((n2-1)*(s2^2)) ) / (n1+n2-2)   )  ) # Hedges G comparing the presence and absence data
      HedgesG_genrich_iters[j] <- genrich_g #save 
      
      #Effect size for Shannon's Diversity 
      shannondiv_g <- (x1_H-x2_H)/( sqrt( ( ((n1-1)*(s1_H^2)) + ((n2-1)*(s2_H^2)) ) / (n1+n2-2)   )  ) #Hedges G comparing Shannon's Diversity of presence and absence data
      HedgesG_H_iters[j] <- shannondiv_g #save 
      
      #Effect size for Simpson's Dominance 
      dominance_g <- (x1_dom-x2_dom)/( sqrt( ( ((n1-1)*(s1_dom^2)) + ((n2-1)*(s2_dom^2)) ) / (n1+n2-2)   )  ) #Hedges G comparing Simpson's Dominance of presence and absence data
      HedgesG_dom_iters[j] <- dominance_g
      
    }
    
  }
  
  #== save data with errors in main results ==#
  #Generic richness
  results_df$M1_genrich[i] <- mean(M1_genrich_iters, na.rm=TRUE)
  results_df$M1_genrich_sd[i] <- sd(M1_genrich_iters, na.rm=TRUE)
  results_df$M2_genrich[i] <- mean(M2_genrich_iters, na.rm=TRUE)
  results_df$M2_genrich_sd[i] <- sd(M2_genrich_iters, na.rm=TRUE)
  #Hedges g: Effects on generic richness
  results_df$HedgesG_genrich[i] <- mean(HedgesG_genrich_iters, na.rm=TRUE)
  results_df$g_genrich_sd[i] <- sd(HedgesG_genrich_iters, na.rm=TRUE)
  
  
  #Shannon's Diversity (H)
  results_df$M1_H[i] <- mean(M1_H_iters, na.rm=TRUE)
  results_df$M1_H_sd[i] <- sd(M1_H_iters, na.rm=TRUE)
  results_df$M2_H[i] <- mean(M2_H_iters, na.rm=TRUE)
  results_df$M2_H_sd[i] <- sd(M2_H_iters, na.rm=TRUE)
  #Hedges g: Effect on Shannon's Diversity 
  results_df$HedgesG_H[i] <- mean(HedgesG_H_iters, na.rm=TRUE)
  results_df$g_H_sd[i] <- sd(HedgesG_H_iters, na.rm=TRUE)
  
  
  #Simpson's Dominance (1/D)
  results_df$M1_dom[i] <- mean(M1_dom_iters, na.rm=TRUE)
  results_df$M1_dom_sd[i] <- sd(M1_dom_iters, na.rm=TRUE)
  results_df$M2_dom[i] <- mean(M2_dom_iters, na.rm=TRUE)
  results_df$M2_dom_sd[i] <- sd(M2_dom_iters, na.rm=TRUE)
  #Hedges g: Effect on Simpson's Dominance 
  results_df$HedgesG_Dominance[i] <- mean(HedgesG_dom_iters, na.rm=TRUE)
  results_df$g_dom_sd[i] <- sd(HedgesG_dom_iters, na.rm=TRUE)
  
  #print that this stage is finished
  print(paste('finished stage:', stage_names[i]))
  
  
}

bioturbators_spatial_results_df <- results_df
save(bioturbators_spatial_results_df, file='Output/bioturbators_spatialsubsample.RData')


