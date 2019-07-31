library(tidyverse)
source('spatial_estimators.R')

#inat_species = c('Rudbeckia hirta')
inat_species = c('Rudbeckia hirta','Maianthemum canadense')
inat_years = c(2016,2017)

# Rudbeckia: long flowering season, very large range = moderate spatial gradient, potentially relatively linear
#               bisq_distance = extent of 30x25, 15 = 30 * 0.5 (0.5 is generally the best for this scenario)
#               weibull_box_size =  extent of 30x25, 10 = mean(c(30,20) * .4) (0.4 box size is best here)
# Maianthemum: short flowering season, small range (mostly NE/Great lakes) = weak spatial gradient, potentially non-linear
#               bisq_distance = extent of 35x15, 7 = 35 * 0.2 (0.2 is generally the best for this scenario)
#               weibull_box_size =  extent of 35x15, 10 = mean(c(25,10) * .4) (0.4 box size is best here)
inat_model_params = tribble(
  ~species,         ~site_type, ~idw_power, ~bisq_distance, ~weibull_grid_box_size, ~weibull_grid_buffer, ~weibull_grid_xlimits, ~weibull_grid_ylimits,
  'Rudbeckia hirta', 'GR',      2,          15,             4,                   0,                  c(-100,-70),            c(25,50),
  'Maianthemum canadense', 'DB',2,    5,              7,                   0,                  c(-95,-60),             c(35,50)
)

inat_data = read_csv('~/data/inat/all_observations.csv') %>%
  filter(species %in% inat_species, open_flower==1, data_source=='inaturalist') %>%
  filter(longitude < -65, longitude > -130, latitude > 20, latitude < 60) 

phenocam_data = read_csv('data/phenocam_transition_dates.csv') %>%
  filter(year<=2017, direction=='rising')

phenocam_site_info = phenocam_data %>%
  select(site, site_type, lat, lon) %>% 
  distinct() %>%
  mutate(x=lon, y=lat)

###################################################################

inat_model_estimates = tibble()

for(this_species in inat_species){
  this_species_params = inat_model_params %>%
    filter(species == this_species)
  
  for(this_year in inat_years){
    this_spp_data = inat_data %>%
      filter(species == this_species, year== this_year) %>%
      mutate(x=longitude, y=latitude)
    
    weibull_grid_model = WeibullGridEstimator(doy_points = this_spp_data,
                                              xlimits = this_species_params$weibull_grid_xlimits[[1]],
                                              ylimits = this_species_params$weibull_grid_ylimits[[1]],
                                              stratum_size_x = 5,
                                              stratum_size_y = 3,
                                              boxes_per_stratum = 40,
                                              box_size = this_species_params$weibull_grid_box_size,
                                              edge_buffer = this_species_params$weibull_grid_buffer)
  
    predictions = phenocam_site_info %>%
      filter(site_type == this_species_params$site_type) %>%
      group_by(site) %>%
      do(predict.WeibullGridEstimator(weibull_grid_model, doy_points = ., type = 'onset', se=T)) %>%
      ungroup() %>%
      rename(doy = estimate, 
             doy_high = estimate_upper,
             doy_low  = estimate_lower) %>%
      mutate(inat_species = this_species, 
             year = this_year, 
             source = 'iNat model')
    
    inat_model_estimates = inat_model_estimates %>%
      bind_rows(predictions)
    
  }
}

#write_csv(inat_model_estimates, 'inaturalist_onset_estimates_at_phenocam_sites.csv')




#############################################

phenocam_onset = phenocam_data %>%
  select(site, doy, doy_high, doy_low, year) %>%
  mutate(source = 'Phenocam Transitions')

combined_estimates = phenocam_onset %>%
  bind_rows(inat_model_estimates) %>%
  #gather(source, onset, phenocam_onset_estimate, inat_onset_estimate) %>%
  left_join(phenocam_site_info, by='site') %>%
  mutate(y_nudge = ifelse(source=='iNat model', 0.2,-0.1))

combined_estimates$site_type = factor(combined_estimates$site_type,
                                      levels = c('DB','GR'),
                                      labels = c('Decid. Broadlead\nMaianthemum canadense',
                                                 'Grassland\nRudbeckia hirta'))

ggplot(combined_estimates, aes(x=doy, y=lat + y_nudge, color=as.factor(site))) + 
  geom_errorbarh(aes(xmax = doy_high, xmin=doy_low, linetype=source), size=1.5, height=0.5) +
  scale_linetype_manual(values=c('dashed','solid')) + 
  geom_point(aes(shape=source), size=6, color='black') +
  geom_point(aes(shape=source), size=4) +
  #gthemes::scale_color_colorblind() + 
  #scale_color_brewer(palette = 'Dark2', direction = -1) + 
  facet_wrap(site_type~year, scales='free_y') +
  theme_dark(25) +
  theme(panel.background = element_rect(fill='grey70'),
        strip.background = element_rect(fill='grey40')) +
  labs(y='Latitude',x='Day Of Year (DOY)', 
       color='Phenocam Site', linetype='Estimate Source', shape = 'Estimate Source')



ggplot(combined_estimates, aes(x=phenocam_onset_estimate , y = inat_onset_estimate, color=as.factor(year))) + 
  geom_point()  +
  geom_abline(slope = 1, intercept = 0) +
  #geom_vline(xintercept = 0) + 
  facet_wrap(~method)


ggplot(combined_estimates, aes(x=phenocam_onset_estimate - inat_onset_estimate, y = lat, color=as.factor(year))) + 
  geom_point()  +
  #geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = 0) + 
  facet_wrap(~method)
