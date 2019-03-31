library(tidyverse)
source('flowering_gradient_generator.R')
source('spatial_grid_estimator.R')

estimator_metadata = read_csv('results/models/2019-03-26-1533/estimator_metadata.csv')
m = read_rds(estimator_metadata$model_filename[1])


# a grid  of evenly spaced points to predict with
prediction_grid = expand.grid(x=seq(0,1,by=0.05),
                              y=seq(0,1,by=0.05))
prediction_grid$dat_type = 'Modelled'

predictions = prediction_grid %>%
  bind_cols(predict.PhenologyGridEstimator(model = m, doy_points = .)) %>%
  rename(onset=onset_estimate)

true_gradient = spatialFloweringGrid(start_doy = 90,
                                     flowering_gradient = estimator_metadata$flowering_gradients[1])

true_gradient$data_type = 'Simulated Gradient'

both = predictions %>%
  bind_rows(true_gradient)

ggplot(both, aes(x=x,y=y,fill=onset)) + 
  geom_raster() + 
  scale_fill_viridis_c() +
  facet_wrap(~data_type)

