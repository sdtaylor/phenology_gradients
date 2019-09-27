library(tidyverse)
library(snakecase)
library(janitor)
library(jsonlite)
##############
# Produce the discussion figure comparing sample sizes of different data sources over time.
##############

# note this data is from the initial pull from inaturalist.org, thus contains
# *all* research grade observations for the 2 spp
inat_observations = read_csv('./data/raw_inat_data/inaturalist_metadata.csv') %>%
  mutate(year = lubridate::year(date),
         data_source = 'iNaturalist',
         species = stringr::word(scientific_name, 1,2)) %>%
  select(species, year, latitude, longitude, data_source)

npn_observations = read_csv('./data/npn_data/status_intensity_observation_data.csv.zip') %>%
  filter(Phenophase_ID == 501) %>%
  mutate(year = lubridate::year(Observation_Date),
         species = paste(Genus, Species),
         data_source = 'USA-NPN') %>%
  select(species, year, latitude=Latitude, longitude=Longitude, data_source)

idigbio_observations = read_csv('./data/idigbio_data/maianthemum/occurrence.csv.zip') %>%
  bind_rows(read_csv('./data/idigbio_data/rudbeckia/occurrence.csv', col_types = cols('dwc:fieldNumber'='c'))) %>%
  janitor::clean_names() %>% # ugh, there are : in the column names
  filter(dwc_basis_of_record == 'preservedspecimen') %>% # herbarium specimens only
  mutate(species = snakecase::to_sentence_case(gbif_canonical_name),
         year = lubridate::year(idigbio_event_date),
         data_source = 'Herbariums') %>%
  select(species, year, idigbio_geo_point, data_source)


all_data = inat_observations %>%
  bind_rows(npn_observations) %>%
  bind_rows(idigbio_observations) %>%
  count(data_source, species, year) 

all_data$species = factor(all_data$species, levels = c('Rudbeckia hirta','Maianthemum canadense'),
                                            labels = c("bold(A.)~italic('Rudbeckia hirta')","bold(B.)~italic('Maianthemum canadense')"))

discussion_figure = ggplot(all_data, aes(x=year, y=n, color=data_source)) + 
  geom_line(size=2) +
  geom_point(size=4) +
  #scale_y_log10(limits = c(1, 5000)) + 
  scale_y_continuous(breaks = c(0,250,500,1000,1500,2000), minor_breaks = c(125,250, 750, 1250)) + 
  scale_x_continuous(breaks = c(seq(2000,2018,4),c(2019)), limits = c(2000,2019)) +  
  #ggthemes::scale_color_colorblind() +
  scale_color_manual(values = c('grey30', "#009E73", "#E69F00")) + 
  #scale_linetype_manual(values = c('solid','dotted')) + 
  facet_wrap(~species, ncol=1, labeller = label_parsed) + 
  theme_bw() +
  theme(legend.position = c(0.31, 0.85),
        legend.background = element_rect(color='black', size=0.5),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size=16),
        panel.background  = element_rect(color='black', size=0.7),
        axis.text.y = element_text(size=12, color='black'),
        axis.text.x = element_text(size=14, color='black'),
        axis.title = element_text(size=20, color='black')) +
  labs(y='Samples Per Year', x='Year') +
  guides(color=guide_legend(title = 'Data Source', override.aes = list(size=3), keywidth = unit(40,'mm'), reverse = T))

ggsave(plot = discussion_figure, filename = 'manuscript/figs/fig5_data_sample_sizes.png', height=14, width = 14, units = 'cm', dpi = 300)

##########################################################
# Map of all observations, including phenocams

# fix the idigbio lat/longs which are stored in, ugh, a json string...
extract_json_lat_lon = function(single_row){
  if(is.na(single_row$idigbio_geo_point)){
    extracted = list(lat=NA, lon=NA)
  } else {
    extracted = jsonlite::fromJSON(single_row$idigbio_geo_point)
  }
  single_row$latitude = extracted$lat
  single_row$longitude = extracted$lon
  return(single_row)
}

idigbio_observations = idigbio_observations %>%
  mutate(row_id = 1:n()) %>%
  group_split(row_id) %>%
  map_dfr(.f = extract_json_lat_lon) %>%
  select(-idigbio_geo_point, -row_id)
################################


all_observations = idigbio_observations %>%
  bind_rows(inat_observations) %>%
  bind_rows(npn_observations) %>%
  mutate(data_source = factor(data_source, levels = c('iNaturalist','Herbariums','USA-NPN'), ordered = T))

basemap = map_data('world') %>%
  filter(region %in% c('USA','Canada'))

lat_range = c(20,60)
lon_range = c(-150,50)

ggplot() + 
  geom_polygon(data = basemap, aes(x=long, y = lat, group = group), fill=NA, color='black', size=1.5) +
  geom_point(data = all_observations, aes(x=longitude, y=latitude, color=data_source), size=2, alpha=0.6) + 
  scale_color_manual(values = c("#009E73",'grey30', "#E69F00")) + 
  #scale_x_continuous(breaks=seq(lon_range[1],lon_range[2],10), minor_breaks = seq(lon_range[1],lon_range[2],1)) +
  #scale_y_continuous(breaks=seq(lat_range[1],lat_range[2],10), minor_breaks = seq(lat_range[1],lat_range[2],1)) + 
  theme_bw() +
  coord_fixed(1.3, xlim = c(-105, -50), ylim=c(24, 55)) +  
  facet_wrap(species~data_source) + 
  theme(panel.background = element_rect(fill='white'),
        panel.grid.major = element_line(color='#0072B2', size=0.5),
        panel.grid.minor = element_line(color='#56B4E9', size=0.2),
        axis.text = element_text(size=20),
        axis.title = element_text(size=25),
        strip.text = element_text(size=18),
        legend.title = element_text(size=25),
        legend.text = element_text(size=20),
        legend.background = element_rect(size=1, color='black'),
        legend.position = 'bottom') +
  labs(x='Longitude',y='Latitude', color='') +
  guides(color = guide_legend(override.aes = list(size=5, alpha=1)))

