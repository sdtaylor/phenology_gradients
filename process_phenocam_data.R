library(phenocamr)
library(tidyverse)

phenocam_site_info = phenocamr::list_sites() %>%
  select(site, lat,lon)

selected_sites = bind_rows(
  tibble(site_name= c('marena','konza','uiefprairie','ninemileprairie','oakville','harvardfarmnorth'),
         veg_type = 'GR'),
  tibble(site_name = c('shiningrock','harvard','proctor','willowcreek','boundarywaters','shalehillsczo'),
         veg_type = 'DB')
)

# for(i in 1:nrow(selected_sites)){
#   try(
#   download_phenocam(site = paste0(selected_sites$site_name[i],'$'),
#                     veg_type = selected_sites$veg_type[i],
#                     roi_id = 1000,
#                     frequency = 3,
#                     phenophase = T,
#                     out_dir = 'data/phenocam_site_data/')
#   )
# }

# All ROI's are = to 1000, except for these which have some others
# download_phenocam(site = 'ninemileprairie$',
#                   veg_type = 'GR',
#                   roi_id = 2000,
#                   frequency = 3,
#                   phenophase = T,
#                   out_dir = 'data/phenocam_site_data/')
# download_phenocam(site = 'boundarywaters$',
#                   veg_type = 'DB',
#                   roi_id = 4000,
#                   frequency = 3,
#                   phenophase = T,
#                   out_dir = 'data/phenocam_site_data/')
# download_phenocam(site = 'shiningrock$',
#                   veg_type = 'DB',
#                   roi_id = 3000,
#                   frequency = 3,
#                   phenophase = T,
#                   out_dir = 'data/phenocam_site_data/')
# download_phenocam(site = 'shalehillsczo$',
#                   veg_type = 'DB',
#                   roi_id = 2000,
#                   frequency = 3,
#                   phenophase = T,
#                   out_dir = 'data/phenocam_site_data/')
#############################
# process the phenocamr derived transition dates

transition_date_files = list.files('data/phenocam_site_data/', pattern = '*transition*', full.names = T)

read_transition_file = function(f){read_csv(f, skip = 16, col_types = cols(roi_id = col_character()))}
dates = purrr::map_df(transition_date_files, read_transition_file) %>%
  filter(gcc_value=='gcc_90') %>% # this is the GCC 90th percentiel
  select(site, site_type = veg_type, roi_id, direction, transition_10, transition_10_lower_ci, transition_10_upper_ci) %>%
  mutate(doy = lubridate::yday(transition_10),
         doy_high = lubridate::yday(transition_10_upper_ci),
         doy_low = lubridate::yday(transition_10_lower_ci),
         year = lubridate::year(transition_10)) %>%
  left_join(phenocam_site_info, by='site') %>%
  filter(year>=2016)

# In 2019 in Konza the low estimate was Dec 13, 2018. So it needs to be a negative DOY to be
# put in the figure correctly. 
dates = dates %>%
  mutate(doy_low = ifelse(site=='konza' & year==2019, -18, doy_low))

write_csv(dates, 'data/phenocam_transition_dates.csv')
