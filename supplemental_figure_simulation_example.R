library(tidyverse)
library(flowergrids)

##################
# This code creates supplemental figure 1 showing examples of simulated data
##################

# example parameters
fl_length = 30
fl_onset = 90
gradient_strength = 3.36 / 0.1


simulation_seed = 275
cell_size = 0.01

linear_gradient = spatial_flowering_grid(start_doy = fl_onset,
                                         flowering_length = fl_length,
                                         flowering_gradient = gradient_strength,
                                         cell_size = cell_size,
                                         seed = seed,
                                         spatial_gradient_type = 'linear') %>%
  mutate(example_type = 'a')

nonlinear_gradient = spatial_flowering_grid(start_doy = fl_onset,
                                            flowering_length = fl_length,
                                            flowering_gradient = gradient_strength,
                                            cell_size = cell_size,
                                            seed = simulation_seed,
                                            spatial_gradient_type = 'non-linear') %>%
  mutate(example_type = 'b')


example_gradients = linear_gradient %>%
  bind_rows(nonlinear_gradient)

##################

sample_size = 200

non_clustered_samples = flowergrids::spatial_flowering_sampler(n = sample_size,
                                                               flowering_length = 30,
                                                               start_doy = fl_onset,
                                                               flowering_gradient = gradient_strength,
                                                               seed = simulation_seed,
                                                               clustering = FALSE) %>%
  mutate(example_type = 'a')

clustered_samples = flowergrids::spatial_flowering_sampler(n = sample_size,
                                                               flowering_length = 30,
                                                               start_doy = fl_onset,
                                                               flowering_gradient = gradient_strength,
                                                               seed = simulation_seed,
                                                               clustering = TRUE) %>%
  mutate(example_type = 'b')

example_observations = non_clustered_samples %>%
  bind_rows(clustered_samples)

facet_types = c('a','b')
facet_labels = c('A. Linear Gradient\nNon-Clustered Sampling', 'B. Nonlinear Gradient\nClustered Sampling')
example_gradients$example_type = factor(example_gradients$example_type, levels = facet_types, labels = facet_labels)
example_observations$example_type = factor(example_observations$example_type, levels = facet_types, labels = facet_labels)

ggplot(example_gradients, aes(x=x, y=y)) + 
  geom_raster(aes(fill=onset)) +
  geom_point(data = example_observations, size=4, color='black') +
  geom_point(data = example_observations, size=2, color='grey60') +
  scale_fill_viridis_c() +
  facet_wrap(~example_type) +
  theme_bw(20) +
  theme(panel.grid = element_blank()) +
  labs(fill = 'Onset DOY')

########
# Some code to make a *bunch* of nonlinear examples
########

# make_map = function(s){
#   spatial_flowering_grid(start_doy = 90,
#                          flowering_length = 30,
#                          flowering_gradient = 33,
#                          cell_size = 0.01,
#                          seed = s,
#                          spatial_gradient_type = 'non-linear') %>%
#     mutate(seed = s)
# }
# 
# n_maps = 20
# all_maps = map_dfr(sample.int(1000, size=n_maps), .f = make_map)
# 
# 
# ggplot(all_maps, aes(x=x,y=y, fill=onset)) +
#   geom_raster() +
#   scale_fill_viridis_c() +
#   facet_wrap(~seed)
