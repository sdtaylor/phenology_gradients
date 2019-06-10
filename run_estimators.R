library(tidyverse)
source('flowering_gradient_generator.R')
source('spatial_grid_estimator.R')


###################################################
# Config

n_cores = 2

# Sampling parameters
sample_sizes = c(300, 600)
flowering_lengths = c(15,30)
flowering_gradients = c(
  # slope of lm(sos~latitude) from Melaas et al. 2018 / scale  of simulated scale
  3.36/0.1, 
  # Half the above, representing a more uniform spatial gradient
  1.68/0.1)
spatial_gradient_types = c('linear','non-linear')
clustering = c(TRUE, FALSE)
n_bootstrap = 100

# Spatial model parameters
#n_boxess = c(200)
box_sizes = c(0.2, 0.4)
edge_buffers = c(0, 0.1, 0.2)

#all_estimates_file = 'results/all_estimates.csv'
current_time =  format(Sys.time(),'%Y-%m-%d-%H%M')
model_dir = paste0('results/models/',current_time,'/')
dir.create(model_dir, recursive = T)

model_run_note = "
using more realistic gradients with landsat start of season estimates from Melaas et al. 2018 data (10.1002/2017GL076933)
add non-linear spatial gradients
fixed non-linear spatial gradient reproducebility
fixed seed issues
changed mean estimate to median
"

cat(model_run_note, file=paste0(model_dir,'notes.txt'))
#####################################################
# extract the slope estimates from a SpatialGridModel
get_slope_estimates = function(spatial_model,
                               prediction_grid){
  prediction_grid = prediction_grid %>%
    bind_cols(predict.PhenologyGridEstimator(spatial_model, prediction_grid)) %>%
    gather(metric, doy_estimate, onset_estimate, end_estimate, peak_estimate)
  
  prediction_grid %>% 
    group_by(metric) %>% 
    do(broom::tidy(lm(.$doy_estimate ~ .$y))) %>%
    ungroup() %>%
    mutate(term = case_when(
      term == '(Intercept)' ~ 'intercept',
      term == '.$y' ~ 'slope' ))
  
}

# a grid  of evenly spaced points to predict with
prediction_grid = expand.grid(x=seq(0.05,0.95,by=0.05),
                              y=seq(0.05,0.95,by=0.05))

#####################################
# Setup parallel processing
library(doParallel)
cl=makeCluster(n_cores)
registerDoParallel(cl)
####################################

all_parameter_combos = expand.grid(sample_size = sample_sizes,
                                   flowering_lengths = flowering_lengths,
                                   flowering_gradients = flowering_gradients,
                                   spatial_gradient_types = spatial_gradient_types,
                                   clustering = clustering,
                                   
                                   box_size = box_sizes,
                                   edge_buffer = edge_buffers,
                                   
                                   bootstrap_i = 1:n_bootstrap)

all_parameter_combos$model_id = 1:nrow(all_parameter_combos)
all_parameter_combos$model_filename = paste0(model_dir,'estimator_',all_parameter_combos$model_id,'.rds')
all_parameter_combos$seed = runif(nrow(all_parameter_combos), 0, 10e8)

write_csv(all_parameter_combos, paste0(model_dir,'estimator_metadata.csv'))

#all_estimates = foreach(iteration_i = 1:nrow(all_parameter_combos), .combine = bind_rows, .errorhandling = 'remove', .packages = c('readr','dplyr','tidyr','broom')) %dopar% {
all_errors = foreach(iteration_i = 1:nrow(all_parameter_combos), .combine = bind_rows, .errorhandling = 'remove', .packages = c('dplyr','tidyr','broom')) %do% {
    
  this_sample_size = all_parameter_combos$sample_size[iteration_i]
  this_length = all_parameter_combos$flowering_lengths[iteration_i]
  this_gradient = all_parameter_combos$flowering_gradients[iteration_i]
  this_gradient_type = all_parameter_combos$spatial_gradient_types[iteration_i]
  do_clustering = all_parameter_combos$clustering[iteration_i]

  this_box_size = all_parameter_combos$box_size[iteration_i]
  this_edge_buffer = all_parameter_combos$edge_buffer[iteration_i]

  this_seed = all_parameter_combos$seed[iteration_i]  
  bootstrap_i = all_parameter_combos$bootstrap_i[iteration_i]

  model_filename = all_parameter_combos$model_filename[iteration_i]
  model_id = all_parameter_combos$model_id[iteration_i]
  
  simulated_sample_data = spatialFloweringSampler(n=this_sample_size,
                                                  distribution_type = 'a',
                                                  start_doy = 90,
                                                  flowering_length = this_length,
                                                  flowering_gradient = this_gradient,
                                                  spatial_gradient_type = this_gradient_type,
                                                  clustering=do_clustering,
                                                  seed = this_seed)


  # Model estimating onset, peak, and end of flowering across space
  spatial_estimator = PhenologyGridEstimator(doy_points = simulated_sample_data,
                                             boxes_per_stratum = 2,
                                             box_size = this_box_size,
                                             edge_buffer = this_edge_buffer)
  
  write_rds(spatial_estimator, model_filename)
  ########################################
  # Now calculate model errors
  
  # a grid  of evenly spaced points to predict with
  prediction_grid = expand.grid(x=seq(0,1,by=0.05),
                                y=seq(0,1,by=0.05))
  predictions = prediction_grid %>%
    bind_cols(predict.PhenologyGridEstimator(model = spatial_estimator, doy_points = .)) %>%
    select(primary_model = onset_estimate, x, y)
  
  # Fit a simple linear model to the data for a "naive" comparison
  linear_fl_model = lm(doy~y, data=spatial_estimator$doy_points)
  predictions$naive_99 = predict(linear_fl_model, newdata = predictions, interval='prediction', level=0.99)[,2]
  
  true_dates  = spatialFloweringGrid(start_doy = 90,
                                     flowering_gradient = this_gradient,
                                     spatial_gradient_type = this_gradient_type,
                                     seed = this_seed)

  errors = predictions %>%
    gather(model_type, onset_estimate, primary_model, naive_99) %>%
    left_join(true_dates, by=c('x','y')) %>%
    mutate(error = onset_estimate - onset)
  
  errors$model_id = model_id
  
  return(errors)
}

write_csv(all_errors, paste0(model_dir,'all_errors.csv'))
