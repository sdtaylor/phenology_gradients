library(phenocamr)

phenocam_site_info = phenocamr::list_sites() %>%
  select(site, lat,lon)

selected_sites = bind_rows(
  tibble(site_name= c('marena','kansas','konza','uiefprairie','ninemileprairie','oakville','sweetbriargrass','harvardfarmnorth'),
         veg_type = 'GR'),
  tibble(site_name = c('shiningrock','harvard','proctor','sanford','willowcreek','boundarrywaters','shalehillsczo'),
         veg_type = 'DB')
)

# for(i in 1:nrow(selected_sites)){
#   try(
#   download_phenocam(site = selected_sites$site_name[i],
#                     veg_type = selected_sites$veg_type[i],
#                     frequency = 3,
#                     phenophase = T,
#                     out_dir = 'data/phenocam_site_data/')
#   )
# }


#############################
# process the phenocamr derived transition dates

transition_date_files = list.files('data/phenocam_site_data/', pattern = '*transition*', full.names = T)

read_transition_file = function(f){read_csv(f, skip = 16)}
dates = purrr::map_df(transition_date_files, read_transition_file) %>%
  filter(gcc_value=='gcc_90') %>% # this is the GCC 90th percentiel
  select(site, roi_id, direction, transition_10, transition_10_lower_ci, transition_10_upper_ci) %>%
  mutate(doy = lubridate::yday(transition_10),
         doy_high = lubridate::yday(transition_10_upper_ci),
         doy_low = lubridate::yday(transition_10_lower_ci),
         year = lubridate::year(transition_10)) %>%
  left_join(phenocam_site_info, by='site') %>%
  filter(year>=2016)


write_csv(dates, 'data/phenocam_transition_dates.csv')
