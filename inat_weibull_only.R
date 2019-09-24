library(tidyverse)
library(flowergrids)

#inat_species = c('Rudbeckia hirta')
inat_species = c('Rudbeckia hirta','Maianthemum canadense')
inat_years = c(2016,2017,2018,2019)

# Rudbeckia: long flowering season, very large range = moderate spatial gradient, potentially relatively linear
#               weibull_box_size =  extent of 30x25, 10 = mean(c(30,20) * .4) (0.4 box size is best here)
# Maianthemum: short flowering season, small range (mostly NE/Great lakes) = weak spatial gradient, potentially non-linear
#               weibull_box_size =  extent of 35x15, 10 = mean(c(25,10) * .4) (0.4 box size is best here)
inat_model_params = tribble(
  ~species,         ~site_type, ~idw_power, ~bisq_distance, ~weibull_grid_box_size, ~weibull_grid_buffer, ~weibull_grid_xlimits, ~weibull_grid_ylimits,
  'Rudbeckia hirta', 'GR',      2,          15,             4,                   0,                  c(-100,-70),            c(25,50),
  'Maianthemum canadense', 'DB',2,          5,              5,                   0,                  c(-95,-60),             c(35,50)
)

#inat_data = read_csv('~/data/inat/all_observations.csv') %>%
inat_data = read_csv('data/all_inat_data.csv') %>%
  filter(species %in% inat_species, open_flower==1) %>%
  filter(longitude < -65, longitude > -130, latitude > 20, latitude < 60) 

phenocam_data = read_csv('data/phenocam_transition_dates.csv') %>%
  filter(year %in% inat_years, direction=='rising') 

phenocam_site_info = phenocam_data %>%
  select(site, site_type, lat, lon) %>% 
  distinct() %>%
  mutate(x=lon, y=lat)

###################################################################

inat_model_estimates = tibble()
set.seed(1234)
for(this_species in inat_species){
  this_species_params = inat_model_params %>%
    filter(species == this_species)
  
  for(this_year in inat_years){
    this_spp_data = inat_data %>%
      filter(species == this_species, year== this_year) %>%
      mutate(x=longitude, y=latitude)
    
    weibull_grid_model = weibull_grid(doy_points = this_spp_data,
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
      do(predict.weibull_grid(weibull_grid_model, doy_points = ., type = 'onset', se=T)) %>%
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

# Bound the lower estimate to 1, or else it gets really  negative
# inat_model_estimates = inat_model_estimates %>%
#   mutate(doy_low = ifelse(doy_low < 1, 1, doy_low),
#          doy = ifelse(doy < 1, 1, doy))

#write_csv(inat_model_estimates, 'inaturalist_onset_estimates_at_phenocam_sites.csv')



#############################################
phenocam_onset = phenocam_data %>%
  select(site, doy, doy_high, doy_low, year) %>%
  mutate(source = 'Phenocam Transitions')

combined_estimates = phenocam_onset %>%
  bind_rows(inat_model_estimates) %>%
  #gather(source, onset, phenocam_onset_estimate, inat_onset_estimate) %>%
  left_join(phenocam_site_info, by='site') %>%
  mutate(y_nudge = ifelse(site_type=='DB', 0.3, 0.2)) %>%
  mutate(y_nudge = ifelse(source=='iNat model',y_nudge,y_nudge*-1)) 

combined_estimates$site_type = factor(combined_estimates$site_type,
                                      levels = c('DB','GR'),
                                      labels = c('atop("Decid. Broadleaf",italic("M. canadense"))',
                                                 'atop("Grassland",italic("R. hirta"))'))

# Average delay between phenocam greenup and the species flowering across all years, sites
combined_estimates %>%
  select(site, doy, year, site_type, source) %>%
  spread(source, doy) %>%
  rename(inat = 'iNat model', phenocam = 'Phenocam Transitions') %>%
  mutate(difference = inat - phenocam) %>%
  group_by(site_type) %>%
  summarise(mean_diff = mean(difference, na.rm=T))


# Sample sizes for inat data
inat_yearly_summaries = inat_data %>%
  group_by(species, year) %>%
  summarise(n=n(),
            mean_doy = mean(doy),
            yearly_onset_estimate = phest::weib.limit(doy, k=50)[1]) %>%
  ungroup() %>%
  mutate(sample_size_label_y = 47,
         sample_size_label_x = ifelse(year==2019,40,250)) %>%
  mutate(sample_size_label = paste0('n=',n))

# match up with the site type, species labels for the faceting
inat_yearly_summaries$site_type = factor(inat_yearly_summaries$species,
                                      levels = c('Maianthemum canadense','Rudbeckia hirta'),
                                      labels = c('atop("Decid. Broadleaf",italic("M. canadense"))',
                                                 'atop("Grassland",italic("R. hirta"))'))


site_labels = combined_estimates %>%
  filter(year==2019, source == 'iNat model') %>%
  mutate(doy = 199, 
         y_nudge=case_when( # Nudge these 2 site labels just a bit or else they overlap.
           site=='konza' ~ -0.15,
           site=='ninemileprairie' ~ 0.35,
           TRUE ~ 0
         ))

ggplot(combined_estimates, aes(x=doy, y=lat + y_nudge, color=interaction(site, site_type))) + 
  geom_errorbarh(aes(xmax = doy_high, xmin=doy_low), size=1.5, height=0) +
  scale_linetype_manual(values=c('dashed','solid')) + 
  geom_point(aes(shape=source), size=6) +
  geom_point(aes(shape=source), size=3, color='white') + 
  scale_shape_manual(values = c(16,17)) +
  scale_color_manual(values = c("#56B4E9","#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7","#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")) + 
  geom_label(data = inat_yearly_summaries, 
             aes(x=sample_size_label_x, y=sample_size_label_y, label = sample_size_label), 
             size=5,
             inherit.aes = F) + 
  geom_vline(data = inat_yearly_summaries, aes(xintercept = mean_doy),
             linetype='dotted', size=2) +
  geom_label(data = site_labels, aes(label=site), 
             size=4,hjust=0, fontface='bold') + 
  #ggthemes::scale_color_gdocs() + 
  #scale_color_brewer(palette = 'Dark2', direction = -1) + 
  scale_x_continuous(breaks=c(1,100,200)) + 
  coord_cartesian(xlim = c(1,310)) + 
  facet_grid(site_type~year, scales='free', labeller = label_parsed) +
  theme_bw(25) +
  theme(legend.position = 'bottom') +
  guides(shape = guide_legend(),
         color = FALSE) + 
  labs(y='Latitude',x='Day Of Year (DOY)', 
       color='Phenocam Site', linetype='Estimate Source', shape = 'Estimate Source')
