library(tidyverse)
source('flowering_gradient_generator.R')
source('spatial_grid_estimator.R')


###################################################
# Config

n_cores = 2

# Sampling parameters
sample_sizes = c(300, 600)
flowering_lengths = c(30)
flowering_gradients = c(10/0.1, 30/0.1)
clustering = c(TRUE)
n_bootstrap = 4

# Spatial model parameters
#n_boxess = c(200)
box_sizes = c(0.2, 0.4)
edge_buffers = c(0, 0.1, 0.2)

#all_estimates_file = 'results/all_estimates.csv'
current_time =  format(Sys.time(),'%Y-%m-%d-%H%M')
model_dir = paste0('results/models/',current_time,'/')
dir.create(model_dir, recursive = T)

model_run_note = "
running for the first time w/ the grid estimator using stratified random sampling.
not testing non-clustered samples.
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
                                   clustering = clustering,
                                   
                                   box_size = box_sizes,
                                   edge_buffer = edge_buffers,
                                   
                                   bootstrap_i = 1:n_bootstrap)

all_parameter_combos$model_id = 1:nrow(all_parameter_combos)
all_parameter_combos$model_filename = paste0(model_dir,'estimator_',all_parameter_combos$model_id,'.rds')

#all_estimates = foreach(iteration_i = 1:nrow(all_parameter_combos), .combine = bind_rows, .errorhandling = 'remove', .packages = c('dplyr','tidyr','broom')) %dopar% {
all_estimates = foreach(iteration_i = 1:nrow(all_parameter_combos), .combine = bind_rows, .errorhandling = 'remove', .packages = c('dplyr','tidyr','broom')) %do% {
    
  this_iteration_estimates = tibble()
  
  this_sample_size = all_parameter_combos$sample_size[iteration_i]
  this_length = all_parameter_combos$flowering_lengths[iteration_i]
  this_gradient = all_parameter_combos$flowering_gradients[iteration_i]
  do_clustering = all_parameter_combos$clustering[iteration_i]

  this_box_size = all_parameter_combos$box_size[iteration_i]
  this_edge_buffer = all_parameter_combos$edge_buffer[iteration_i]
  
  bootstrap_i = all_parameter_combos$bootstrap_i[iteration_i]

  model_filename = all_parameter_combos$model_filename[iteration_i]
  
  simulated_sample_data = spatialFloweringSampler(n=this_sample_size,
                                                  distribution_type = 'a',
                                                  start_doy = 90,
                                                  flowering_length = this_length,
                                                  flowering_gradient = this_gradient,
                                                  clustering=do_clustering)


  # Model estimating onset, peak, and end of flowering across space
  spatial_estimator = PhenologyGridEstimator(doy_points = simulated_sample_data,
                                             boxes_per_stratum = 2,
                                             box_size = this_box_size,
                                             edge_buffer = this_edge_buffer)
  
  write_rds(spatial_estimator, model_filename)
  # slope_estimates = get_slope_estimates(spatial_estimator,
  #                                       prediction_grid = prediction_grid)
  # 
  # slope_estimates$sample_size = this_sample_size
  # slope_estimates$flowering_length = this_length
  # slope_estimates$flowering_gradient = this_gradient
  # slope_estimates$clustering = do_clustering
  # slope_estimates$bootstrap_i = bootstrap_i
  # slope_estimates$n_boxes = this_n_boxes
  # slope_estimates$box_size = this_box_size
  # slope_estimates$edge_buffer = this_edge_buffer
  # 
  # this_iteration_estimates = this_iteration_estimates %>%
  #   bind_rows(slope_estimates)
  
  return(this_iteration_estimates)
}

#write_csv(all_estimates, all_estimates_file)
write_csv(all_parameter_combos, paste0(model_dir,'estimator_metadata.csv'))