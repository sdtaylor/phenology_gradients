library(tidyverse)
source('flowering_gradient_generator.R')
source('spatial_grid_estimator.R')


###################################################
# Config

n_cores = 10

# Sampling parameters
sample_sizes = c(300, 600)
flowering_lengths = c(15,30)
flowering_gradients = c(10/0.1, 20/0.1, 30/0.1)
clustering = c(TRUE,FALSE)
n_bootstrap = 2

# Spatial model parameters
n_boxess = c(200)
box_sizes = c(0.3, 0.6)
edge_buffers = c(0, 0.1, 0.2)

all_estimates_file = 'results/all_estimates.csv'
###################################################
###################################################
generate_clustered_points = function(n, uniform_random_percent=0.2){
  uniform_n = floor(n * uniform_random_percent)
  clumped_n = n - uniform_n
  point_pattern = spatstat::rThomas(5, 0.05, clumped_n*1.5)
  
  # Spatstat will generate a random number of points, so from
  # those choose the n points requested.
  number_of_random_points = length(point_pattern$x)
  if(number_of_random_points < clumped_n){
    # If not enough random points were generated try again
    return(generate_clustered_points(clumped_n, uniform_random_percent = uniform_random_percent))
  } else {
    uniform_x = runif(n=uniform_n, min=0, max=1)
    uniform_y = runif(n=uniform_n, min=0, max=1)
    random_selection = sample(1:number_of_random_points, size = clumped_n,replace = FALSE)
    x = c(point_pattern$x[random_selection], uniform_x)
    y = c(point_pattern$y[random_selection], uniform_y)
    return(tibble(x=x, y=y))
  }
}


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
                                   bootstrap_i = 1:n_bootstrap)


all_estimates = foreach(iteration_i = 1:nrow(all_parameter_combos), .combine = bind_rows, .errorhandling = 'remove', .packages = c('dplyr','tidyr','broom')) %dopar% {
#all_estimates = foreach(iteration_i = 1:nrow(all_parameter_combos), .combine = bind_rows, .errorhandling = 'remove', .packages = c('dplyr','tidyr','broom')) %do% {
    
  this_iteration_estimates = tibble()
  
  this_sample_size = all_parameter_combos$sample_size[iteration_i]
  this_length = all_parameter_combos$flowering_lengths[iteration_i]
  this_gradient = all_parameter_combos$flowering_gradients[iteration_i]
  do_clustering = all_parameter_combos$clustering[iteration_i]
  bootstrap_i = all_parameter_combos$bootstrap_i[iteration_i]

  # The flowering sampler will generate completely random, and spatially
  # uniform, points by default. If doing clustering this will pass in
  # clusted points
  if(do_clustering){
    clustered_points = generate_clustered_points(this_sample_size)
    x = clustered_points$x
    y = clustered_points$y
  } else {
    x = NULL
    y = NULL
  }
  
  simulated_sample_data = spatialFloweringSampler(n=this_sample_size,
                                                  distribution_type = 'a',
                                                  start_doy = 90,
                                                  flowering_length = this_length,
                                                  flowering_gradient = this_gradient,
                                                  x=x,y=y)

  for(this_n_boxes in n_boxess){
    for(this_box_size in box_sizes){
      for(this_edge_buffer in edge_buffers){
      
        # Model estimating onset, peak, and end of flowering across space
        spatial_estimator = PhenologyGridEstimator(doy_points = simulated_sample_data,
                                                   n_boxes = this_n_boxes,
                                                   max_box_size = this_box_size, min_box_size = this_box_size,
                                                   edge_buffer = this_edge_buffer)
        
        slope_estimates = get_slope_estimates(spatial_estimator,
                                              prediction_grid = prediction_grid)
        
        slope_estimates$sample_size = this_sample_size
        slope_estimates$flowering_length = this_length
        slope_estimates$flowering_gradient = this_gradient
        slope_estimates$clustering = do_clustering
        slope_estimates$bootstrap_i = bootstrap_i
        slope_estimates$n_boxes = this_n_boxes
        slope_estimates$box_size = this_box_size
        slope_estimates$edge_buffer = this_edge_buffer
        
        this_iteration_estimates = this_iteration_estimates %>%
          bind_rows(slope_estimates)
      }
      }
  }
  return(this_iteration_estimates)
}

write_csv(all_estimates, all_estimates_file)