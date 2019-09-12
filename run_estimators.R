library(tidyverse)
#devtools::install_github('sdtaylor/flowergrids')

source('config.R')

###################################################
# Config

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

added other models
"

cat(model_run_note, file=paste0(model_dir,'notes.txt'))

#####################################
# Setup parallel processing
library(doParallel)
cl=makeCluster(n_cores)
registerDoParallel(cl)
####################################

simulation_combos = expand.grid(sample_size = sample_sizes,
                                   flowering_lengths = flowering_lengths,
                                   flowering_gradients = flowering_gradients,
                                   spatial_gradient_types = spatial_gradient_types,
                                   clustering = clustering,
                                   bootstrap_i = 1:n_bootstrap)

simulation_combos$simulation_id = 1:nrow(simulation_combos)
simulation_combos$seed = runif(nrow(simulation_combos), 0, 10e8)

write_csv(simulation_combos, paste0(model_dir,'simulation_metadata.csv'))

required_packages =  c('dplyr','tidyr','broom','spatstat','phest')

#all_estimates = foreach(simulation_i = 1:nrow(all_parameter_combos), .combine = bind_rows, .errorhandling = 'remove', .packages = required_packages)) %dopar% {
all_errors = foreach(simulation_i = 1:nrow(simulation_combos), .combine = bind_rows, .errorhandling = 'remove', .packages =required_packages) %do% {
    
  this_sample_size = simulation_combos$sample_size[simulation_i]
  this_length = simulation_combos$flowering_lengths[simulation_i]
  this_gradient = simulation_combos$flowering_gradients[simulation_i]
  this_gradient_type = simulation_combos$spatial_gradient_types[simulation_i]
  do_clustering = simulation_combos$clustering[simulation_i]

  this_seed = simulation_combos$seed[simulation_i]  
  bootstrap_i = simulation_combos$bootstrap_i[simulation_i]
  
  simulation_id = simulation_combos$simulation_id[simulation_i]
  
  simulated_sample_data = flowergrids::spatial_flowering_sampler(n=this_sample_size,
                                                                 distribution_type = 'a',
                                                                 start_doy = 90,
                                                                 flowering_length = this_length,
                                                                 flowering_gradient = this_gradient,
                                                                 spatial_gradient_type = this_gradient_type,
                                                                 clustering=do_clustering,
                                                                 seed = this_seed)
  
  true_dates  = flowergrids::spatial_flowering_grid(start_doy = 90,
                                                    xlimits = c(0,1),
                                                    ylimits = c(0,1),
                                                    cell_size = 0.05,
                                                    flowering_gradient = this_gradient,
                                                    spatial_gradient_type = this_gradient_type,
                                                    seed = this_seed)
  
  # a grid  of evenly spaced points to predict with. matches true_dates
  prediction_grid = expand.grid(x=seq(0,1,by=0.05),
                                y=seq(0,1,by=0.05))
  
  model_estimates = tibble()
  #############################################
  # The Weibull grid estimator
  #############################################
  # Model estimating onset, peak, and end of flowering across space
  for(model_i in weibull_model_parameters$model_id){
    params = filter(weibull_model_parameters, model_id == model_i)
    weibull_grid_model = flowergrids::weibull_grid(doy_points = simulated_sample_data,
                                                   boxes_per_stratum = 5,
                                                   box_size = params$box_size,
                                                   edge_buffer = params$edge_buffer)
  
    predictions = prediction_grid %>%
      mutate(onset_estimate = flowergrids::predict.weibull_grid(model = weibull_grid_model, doy_points = ., type='onset')) %>%
      mutate(model_id = model_i,
             model = 'weibull_grid')
    
    model_estimates = model_estimates %>%
      bind_rows(predictions)
  }

  #############################################
  # The linear model
  #############################################
  linear_fl_model = lm(doy~y, data=simulated_sample_data)
  predictions = prediction_grid %>%
    mutate(onset_estimate = predict(linear_fl_model, newdata = predictions, interval='prediction', level=0.99)[,2]) %>%
    mutate(model_id = 1,
           model = 'linear')
  
  model_estimates = model_estimates %>%
    bind_rows(predictions)
  #############################################

  # This makes a cool figure showing rasters for all estimators + the true raster
  # true_dates %>%
  #   mutate(model_i = 1, model = 'true_dates') %>%
  #   rename(onset_estimate = onset) %>%
  #   bind_rows(model_estimates) %>%
  # ggplot(aes(x=x,y=y, fill=onset_estimate)) + 
  #   geom_raster() + 
  #   scale_fill_viridis_c() + 
  #   facet_wrap(model~model_id)
  
  errors = model_estimates %>%
    left_join(true_dates, by=c('x','y')) %>%
    group_by(model, model_id) %>%
    summarise(rmse = sqrt(mean((onset_estimate - onset)**2, na.rm=T)),
              bias = mean(onset_estimate, na.rm=T) - mean(onset, na.rm=T),
              n=n(),
              num_na = sum(is.na(onset_estimate)),
              percent_na = num_na/n) %>%
    ungroup()
  
  errors$simulation_id = simulation_id
  
  return(errors)
}

write_csv(all_errors, paste0(model_dir,'all_errors.csv'))
