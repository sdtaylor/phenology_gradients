##############################################################
# This file holds:
# -the paramters for the simulated data, ie the underlying phenology parameters, sample sizes, etc.
# -the different parameters used for each of the methods in the simulation study

####################################
# global config stuff
####################################
photo_folder = '/home/shawn/data/phenology_gradients/photos/'


###################################
# stuff affecting runtim and resources in run_estimators.R
###################################
n_cores = 1
n_bootstrap = 50

###################################
# underlying phenology parameters
###################################
flowering_lengths = c(15,30,45,60)
flowering_gradients = c(
  # slope of lm(sos~latitude) from Melaas et al. 2018 / scale  of simulated scale
  3.36/0.1, 
  # double the above
  6.72/0.1,
  # Half the above, representing a relativly uniform spatial gradient
  1.68/0.1)
spatial_gradient_types = c('linear','non-linear')
###################################
# sampling parameters
###################################
sample_sizes = c(150, 300, 600, 1200)
clustering = c(TRUE, FALSE)

########################################
# Spatial model parameters
######################################
#n_boxess = c(200)
weibull_model_parameters = expand.grid(
  box_size = c(0.2, 0.4),
  num_boxes = c(5,10,20,40),
  stratum_size = c(0.1, 0.2, 0.5)
)
weibull_model_parameters$model_id = 1:nrow(weibull_model_parameters)

# bi-squared weighted distance (bisq) model
bisq_model_parameters = expand.grid(
  b_distance = c(0.2, 0.5, 0.75, 1.0)
)
bisq_model_parameters$model_id = 1:nrow(bisq_model_parameters)


# inverse weighted distance (idw) model
idw_model_parameters = expand.grid(
  idw_power = c(1,2,3)
)
idw_model_parameters$model_id = 1:nrow(idw_model_parameters)


# # linear model 
# linear_model_parameters = expand.grid(
#   quantile = c(0.9, 0.95, 0.99)
# )
# linear_model_parameters$model_id = 1:nrow(linear_model_parameters)
# 
