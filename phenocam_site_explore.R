library(tidyverse)



library(phenor)
library(phenocamr)
library(ggrepel)
download_phenocam(veg_type = "DB",
                  site = "acadia",
                  phenophase = TRUE)

d = format_phenocam()

phenocam_sites = phenocamr::list_sites()
grassland_sites = phenocam_sites %>% 
  filter(primary_veg_type=='GR')

forest_sites = phenocam_sites %>% 
  filter(primary_veg_type%in% c('DB','EN')) %>%
  filter(date_start < '2015-01-01', date_end > '2019-01-01')

basemap = map_data('state')

ggplot() + 
  geom_polygon(data = basemap, aes(x=long, y = lat, group = group), fill=NA, color='black', size=1.5) +
  geom_point(data=forest_sites, aes(x=lon, y=lat), size=2.5) + 
  geom_label_repel(data=forest_sites, aes(x=lon, y=lat,label=site), size=2.5) + 
  theme_bw() +
  #ylim(35,55)+
  #xlim(-100,-65)+
  #coord_fixed(1.3, xlim=c(-120,-65), ylim=c(25,70)) +  
  theme(panel.background = element_rect(fill='white'),
        axis.text = element_text(size=25),
        axis.title = element_text(size=30)) +
  labs(x='Longitude',y='Latitude')


# grassland sites I chose by hand
# 'tuckerprairie' and 'twosfpr' also look good but they dont start till mid 2018
#  Notes on ROI's.
# marena: single ROI of grass
# kansas: single ROI of grass
# konza: singel ROI of grass
# uiefprairie: single ROI of grass
# ninemileprairie: DB_1000 and DB_2000 are trees, GR_1000 and  GR_2000 are the same grass patch
# oakville: single ROI of grass
# sweetbriargrass: single ROI of grass
# harvardfarmnorth: hmm, 2 ROIS, says they're different but they look the same. use GR_1000
selected_gr_sites = c('marena','kansas','konza','uiefprairie','ninemileprairie','oakville','sweetbriargrass','harvardfarmnorth')

# forest sites I chose by hand
# there's a *bunch*, so inital filter is within latitude (35,55), longitude (-100,-65) and has tart dates < 2015 and
# end date > 2019 (ie. up to the present).
# Notes on ROI's. Pick 1000 in all forest cases. 
# shiningrock: all 3 ROIs (1000,2000,3000) are of the same general area. suspect they're changing from moving camera
# harvard: all 2 ROIs (0001, 1000) are of same areas. 
# proctor: single ROI of trees
# sanford: single ROI of trees
# willowcreek: single ROI of trees (and sensor)
# boundarrywaters: all 4 ROIs (1000, 2000, 3000, 4000) are of individual trees, just pick the 1st. 
# shalehillsczo: both are different areas of trees, very similar. 
selected_forest_sites = c('shiningrock','harvard','proctor','sanford','willowcreek','boundarrywaters','shalehillsczo')
selected_sites_info = phenocam_sites %>%
  filter(site %in% selected_sites)
