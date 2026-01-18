# Marine Phanerozoic biodiversity increased in presence of ecosystem engineers

<a href="https://doi.org/10.21203/rs.3.rs-5447601/v1"><img src="https://img.shields.io/badge/Download preprint here!-ffb703"></a>

Analyses associated with the manuscript "Marine Phanerozoic biodiversity increased in presence of ecosystem engineers" by AT Cribb, SAF Darroch, and W Gearty. <b>This manuscript is currently in revision and is in the process of being updated in accordance with reviewer suggestions. Do not use the scripts without contacting me first.</b>

Contact: A.T.Cribb@soton.ac.uk

<img src="https://img.shields.io/badge/README%20is%20under%20construction-ff5400"> <img src="https://img.shields.io/badge/unpublished:-in revision-7678ed"> <a href="https://doi.org/10.5281/zenodo.14196644"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14196644.svg" alt="DOI"></a>

This repository contains all data and R scripts neede to reproduce the analyses and figures in the manuscript main text and suppelementary materials. You should run the analyses to generate new Output files in the Output folder created here. If you download this repository, set it as your working directory and everything should run smoothly. Should you have any questions or problems, or if you wish to make modifications for your own publication, please contact me at the email address above.

## Datasets
There are three .RData files in <kbd>[Data](https://github.com/atcribb/Ecosystem-Engineers-Biodiversity/tree/main/Data)</kbd>:
* ``Bioturbators_data.RData`` - contains the bioturbator ecosystem engineering data
* ``Reef_Ecosystem_Engineers_Final.RData`` - contains the reef-builder ecosystem engineering data
* ``Phanerozoic_clean_final.Rdata`` - contains the cleaned and stratigraphically binned (largely following Kocsis et al. 2019 ddPhanero protocol with additional cleaning of formation data) Phanerozoic data. The raw PBDB data was downloaded 1 November 2023. Using this dataset will reproduce the figures in the paper. 

There is also a ``scotese_climate`` folder that contains a .csv file ``Scotese-2021_GATs1Myr.csv``. These are global average temperatures for every one million years over the Phanerozoic extracted from <a href="https://doi.org/10.1016/j.earscirev.2021.103503"><b>Scotese et al. (2021) Phanerozoic paleotemperatures: The earth's changing climate during the last 540 million years</b></a>. Please do not use this data without also appropriately citing Scotese et al. (2021). 

## Analyses
There are two sets of <kbd>[Analyses](https://github.com/atcribb/Ecosystem-Engineers-Biodiversity/tree/main/Data)</kbd>: the main effect text analyses, and the supplemental analyses.

### Effect sizes 
* ``EffectSize_Bioturbation.R`` - calculates Hedges g effect sizes for bioturbation on biodiversity 
* ``EffectSize_Reefs.R`` - calculates Hedges g effect size for reef-buidlers on biodiversity
* Both of these scripts contain three different subsampling protocols based on how or if formations are subsampled. By default, these are set to 20 occurrences per formation beacuse that is what is used in the main analysis, but by un-commenting either of the other two options and the corresponding file save line at the end of the script, you can run the analyses with the other two subsampling methods we test. You will need to do all three for the sampling biases analyses. Save all data to your /Output folder that you created in your working repository.
* Also, note that these scripts calculate richness, Shannon's Diversity (H), and Simpson's Dominance. Only Shannon's Diversity is discussed in the manuscript, because it accounts for both richness and evenness together, but we have left each diversity metric in for any interested users. 
* ``Biodiversity_bioturbators.R`` - calculates 1. a gloabl generic richness curve of bioturbators, and 2. the number of genera whose FAD are in each stage. Data is subsampled to 150 occurrences per formation, which requires the removal of stages that do not meet that occurrence quota. 
* ``Biodiversity_reefs.R`` - the same thing as above, but for reef-builders.

### Supplementary analyses
There are twelve scripts needed to produce the supplementary figures:
* ``EffectSize_Bioturbation_Bathymetry_Deep.R`` and ``EffectSize_Reefs_Bathymetry_Deep.R`` - Effect size analyses for only deep marine palaeoenvironments
* ``EffectSize_Bioturbation_Bathymetry_Shallow.R`` and ``EffectSize_Reefs_Bathymetry_Shallow.R`` - Effect size analyses for only shallow marine palaeoenvironments
* ``EffectSize_Bioturbation_Environment_Carbonate.R`` and ``EffectSize_Reefs_Environment_Carbonate.R`` - Effect size analyses for only carbonate facies
* ``EffectSize_Bioturbation_Environment_Siliciclastic.R`` and ``EffectSize_Reefs_Environment_Siliciclastic.R`` - Effect size analyses for only siliciclastic facies 
* ``EffectSize_Bioturbators_RemoveEEs.R`` and ``EffectSize_Reefs_RemoveEEs.R`` - effect size analyses if ecosystem engineers are removed from calculation of diversity metrics 
* ``Bioturbation_SamplingBiases.R`` - used to determine which formation subsampling method best reduces biases of sampling effort on effect sizes for bioturbators 
* ``Reef_SamplingBiases.R`` - uesd to determine which formation subsampling method best reduces biases of sampling effort on effect sizes for reef-builders
* These analyses will reproduce supplemenetal figures x and y. You need to have run all three subsampling versions of the effect size analyses to use the sampling biases analyses scripts. 

## Plotting outputs 
* By default, the plotting output scripts are set to use the results from subsampling 20 occurrences per formation, as this is what is used in the main text. Change what file you load if you want to see results for different subsampling methods.
* ``Figures_Bioturbation_EffectSizes.R`` - plotting outputs for Figure 1 comparing biodiversity measurements per formation with and without bioturbators in each stage and effect size of the presence of bioturbators in each stage. 
* ``Figures_Reef_EffectSizes.R`` - plotting outputs for Figure 2 comparing biodiversity measurements per formation with and without reef-builders in each stage and effect size of the presence of reef-builders in each stage.
* ``Figures_MassExtinctions.R`` - subsets and replots data for each mass extinction to create Figure 3
* ``Figures_Correlates_Diversity.R`` - uses the effect size outputs and ecosystem engineer diveristy outputs to create Figure 4 and Figure 5
* ``Figures_Correaltes_Climate.R`` - uses the effect size outputs and the Scotese et al. (2021) GATs to assess global temperature as a predictor of effect sizes to create Figure 6. 

# Change log
### 10 February 2025
``` diff
+ starting revision process
```

### 29 October 2024
```diff
+ add diversity and climate correlate analyses
+ add all supplementary effect size analyses variations 
+ add new main text figure plotting output scripts 
+ add supplemental figure plotting output scripts
```

### 1-3 June 2024
```diff
+ Update bioturbation dataset to include a small number of infauna which were not previously included 
+ Update sampling bias analysis scripts to plot but linear and polynomial regressions more easily 
```

### 30 May 2024
```diff
+ Minor edits to plotting aesthetics, and fixing issue where effect sizes strength was not being assessed based on uncertainty bounds
```

### 28 May 2024
```diff
-  remove weighed means/standard deviation because of lack of ability to consistently apply in each stage - with large iter, more likely to deal with n1=1 in stages where EEs are not dominant, where sd and thereby weighted means/sd cannot be calculated. Switching back to unweighted mean and standard deviation to opt for consistency. None of this impacts Hedges g.
```

### 15 May 2024
```diff
+ update effect size analyses to calculate weighted means and standard deviations for generic richness and Shannon's Diversity
- remove AIC tests from sampling biases analyses
```

### 25 April 2024
```diff
+ upload full PBBD dataset
+ update README to reflect new dataset upload
+ update README for Analyses and Plotting_Output
- remove DataClean in lieu of uploading PBDB data
```

### 4 April 2024
```diff
+ writing the README
+ uploading Analyses scripts
```

