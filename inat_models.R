library(tidyverse)
source('spatial_estimators.R')

inat_species = c('Rudbeckia hirta','Maianthemum canadense')
inat_years = c(2016,2017)

# Rudbeckia: long flowering season, very large range = moderate spatial gradient, potentially relatively linear
#               bisq_distance = extent of 30x25, 15 = 30 * 0.5 (0.5 is generally the best for this scenario)
#               weibull_box_size =  extent of 30x25, 10 = mean(c(30,20) * .4) (0.4 box size is best here)
# Maianthemum: short flowering season, small range (mostly NE/Great lakes) = weak spatial gradient, potentially non-linear
#               bisq_distance = extent of 35x15, 7 = 35 * 0.2 (0.2 is generally the best for this scenario)
#               weibull_box_size =  extent of 35x15, 10 = mean(c(25,10) * .4) (0.4 box size is best here)
inat_model_params = tribble(
  ~species,         ~idw_power, ~bisq_distance, ~weibull_grid_box_size, ~weibull_grid_buffer, ~weibull_grid_xlimits, ~weibull_grid_ylimits,
  'Rudbeckia hirta', 2,          15,             4,                   0,                  c(-100,-70),            c(25,50)
  #'Maianthemum canadense', 2,    5,              7,                   0,                  c(-95,-60),             c(35,50),
)

inat_data = read_csv('~/data/inat/all_observations.csv') %>%
  filter(species %in% inat_species, open_flower==1, data_source=='inaturalist') %>%
  filter(longitude < -65, longitude > -130, latitude > 20, latitude < 60) 

phenocam_dates = read_csv('data/phenocam_transition_dates.csv') %>%
  filter(year!=2019, direction=='rising')

phenocam_site_info = phenocam_dates %>%
  select(site, lat, lon) %>% 
  distinct() %>%
  mutate(x=lon, y=lat)

###################################################################

all_predictions = tibble()

for(this_species in inat_species){
  this_species_params = inat_model_params %>%
    filter(species == this_species)
  
  for(this_year in inat_years){
    this_spp_data = inat_data %>%
      filter(species == this_species, year== this_year) %>%
      mutate(x=longitude, y=latitude)
    
    idw_model = InterpolationEstimator(doy_points = this_spp_data,
                                       percentile = 0.99,
                                       method='idw',
                                       power = this_species_params$idw_power)
    
    bisq_model = InterpolationEstimator(doy_points = this_spp_data,
                                        percentile = 0.99,
                                        method='bisq',
                                        bisq_distance = this_species_params$bisq_distance)
    
    weibull_grid_model = WeibullGridEstimator(doy_points = this_spp_data,
                                              xlimits = this_species_params$weibull_grid_xlimits[[1]],
                                              ylimits = this_species_params$weibull_grid_ylimits[[1]],
                                              stratum_size_x = 5,
                                              stratum_size_y = 3,
                                              boxes_per_stratum = 20,
                                              box_size = this_species_params$weibull_grid_box_size,
                                              edge_buffer = this_species_params$weibull_grid_buffer)
    
    
    predictions = phenocam_site_info %>%
      mutate(idw = predict.InterpolationEstimator(idw_model, doy_points = ., type = 'onset'),
             bisq = predict.InterpolationEstimator(bisq_model, doy_points = ., type = 'onset'),
             weibull = predict.WeibullGridEstimator(weibull_grid_model, doy_points = ., type = 'onset')) %>%
      select(-lat, -lon, -x, -y) %>%
      gather(method, onset_estimate, idw, bisq, weibull) %>%
      mutate(species = this_species, year = this_year)
    
    all_predictions = all_predictions %>%
      bind_rows(predictions)
    
  }
}

write_csv(all_predictions, 'inaturalist_onset_estimates_at_phenocam_sites.csv')