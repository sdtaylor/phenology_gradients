## Opportunistically collected photographs can be used to estimate large-scale phenological trends

This is the code and primary data for the above study. 

The analysis develops and tests a method, the Weibull Grid, to estimate flowering onset across large-scales using only observations of flowering presence. It uses a Weibull estimator (from the [phest](https://github.com/willpearse/phest) package) integrated into a spatially explicit ensemble. The method is tested with a simulation study and then applied to actually phenology data derived from iNaturalist. 

Code for the Weibull Grid model and making flowering simulations is in the [flowergrids](https://github.com/sdtaylor/flowergrids) R package. 

Data sources:

[iNaturalist](https://www.inaturalist.org)  
[Phenocam](https://phenocam.sr.unh.edu/webcam/)  
[USA National Phenology Network*](https://www.usanpn.org)  
[iDigBio*](https://www.idigbio.org)  
*These data sources were only used to show historic trends in data collection, not in the main analysis.  

### Repository files

`data/` : All original data used in the analysis. 

`results/` : Results from the simulation analysis

`inat_data_prep/` : Scripts and config files for downloading and scoring the iNaturalist photographs for the two species.

`prior_study_phenology` : Notes and data on prior measurements of onset for the two species analysed.

`install_packages.R` : Install the packages used in this analysis

`config.R` : Setup for the simulation analysis parameters. 

`run_estimators.R` : Runs the simulation analysis. 

`analysis_simulation_errors.R` : Produce the figures and tables related to the simulations. 

`conceptual_figures.R` : Create Figure 1 from the introduction.

`phenocam_site_explore.R` : Mapping out the different phenocam sites for selection.

`process_phenocam_data.R` : Downloading and process the phenocam data used in the analysis

`inaturalist_onset_models.R` : Fit the iNaturalist models for the 2 species and create the figure comparing them to phenocam onset

`supplemental_figure_simulation_example.R` : Create the supplemental figures showing the example simulations.

`other_dataset_comparison.R` : Create Figure 5 showing the sample size of different data sources over time.
