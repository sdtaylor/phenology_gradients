library(tidyverse)
source('config.R')

# Combine the direction annotations from ImageAnt with the full iNaturalist metadata
# to get the iNat data used for analysis. 
inat_metadata = read_csv('./data/raw_inat_data/inaturalist_metadata.csv', col_types = cols('inat_id'='c')) %>%
  mutate(species = word(scientific_name, 1,2)) %>%
  select(-scientific_name)

annotations = read_csv('inat_data_prep/imageant_sessions/inat_rudbeckia.csv', col_types = cols('open_flowers'='c','flowers'='c')) %>%
  bind_rows(read_csv('inat_data_prep/imageant_sessions/inat_maianthemum.csv', col_types = cols('open_flowers'='c','flowers'='c'))) %>%
  filter(open_flowers != 'U') %>%
  mutate(open_flower = as.numeric(open_flowers)) %>%
  select(-time, -open_flowers) %>%
  mutate(inat_id = tools::file_path_sans_ext(file))

inat_data = inat_metadata %>%
  left_join(annotations, by=c('inat_id')) %>%
  mutate(year = lubridate::year(date),
         doy  = lubridate::yday(date))

write_csv(inat_data, 'data/processed_inat_data.csv')