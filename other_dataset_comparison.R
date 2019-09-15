library(tidyverse)
library(snakecase)
library(janitor)
library(ggforce)
##############
# Produce the discssion figure comparing sample sizes of different data sources over time.
##############

inat_observations = read_csv('~/data/phenology_gradients/inaturalist_metadata.csv') %>%
  mutate(year = lubridate::year(date),
         data_source = 'iNaturalist',
         species = stringr::word(scientific_name, 1,2)) %>%
  select(species, year, data_source)

npn_observations = read_csv('~/data/phenology_gradients/npn_data/status_intensity_observation_data.csv') %>%
  filter(Phenophase_ID == 501) %>%
  mutate(year = lubridate::year(Observation_Date),
         species = paste(Genus, Species),
         data_source = 'USA-NPN') %>%
  select(species, year, data_source)

idigbio_observations = read_csv('~/data/phenology_gradients/idigbio_data/maianthemum/occurrence.csv') %>%
  bind_rows(read_csv('~/data/phenology_gradients/idigbio_data/rudbeckia/occurrence.csv', col_types = cols('dwc:fieldNumber'='c'))) %>%
  janitor::clean_names() %>% # ugh, there are : in the column names
  filter(dwc_basis_of_record == 'preservedspecimen') %>% # herbarium specimens only
  mutate(species = snakecase::to_sentence_case(gbif_canonical_name),
         year = lubridate::year(idigbio_event_date),
         data_source = 'Herbariums') %>%
  select(species, year, data_source)


all_data = inat_observations %>%
  bind_rows(npn_observations) %>%
  bind_rows(idigbio_observations) %>%
  count(data_source, species, year)


ggplot(all_data, aes(x=year, y=n, color=data_source, linetype = species)) + 
  geom_line(size=2) +
  #scale_y_log10(limits = c(1, 5000)) + 
  scale_y_continuous(breaks = c(0,100,500,1000,1500,2000), minor_breaks = c(50,250, 750, 1250)) + 
  scale_x_continuous(breaks = c(seq(1950,2009,10),c(2010,2015,2019)), limits = c(1950,2019)) +  
  #ggthemes::scale_color_colorblind() +
  scale_color_manual(values = c('grey30', "#009E73", "#E69F00")) + 
  scale_linetype_manual(values = c('solid','dotted')) + 
  #facet_wrap(~species, ncol=1) + 
  theme_bw() +
  theme(legend.position = c(0.2, 0.6),
        legend.background = element_rect(color='black', size=0.5),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18),
        panel.background  = element_rect(color='black', size=0.7),
        axis.text = element_text(size=16, color='black'),
        axis.title = element_text(size=20, color='black')) +
  labs(y='Samples Per Year', x='') +
  guides(color=guide_legend(title = 'Data Source', override.aes = list(size=3), keywidth = unit(40,'mm')),
         linetype = guide_legend(title = 'Species', override.aes = list(size=2), keywidth = unit(40,'mm')))
          