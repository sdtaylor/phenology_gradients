library(tidyverse)
library(ggrepel)

source('config.R')

simulation_metadata = read_csv('results/models/2019-09-16-1739/simulation_metadata.csv')
all_model_errors = read_csv('results/models/2019-09-16-1739/all_errors.csv')

# For every suite of 100 phenology/sampling simulations, get the average error for a particular
# weibull model parameter set (identified by model id), and then find the best parameter set for each unique simulation. 
# For the best ones, do not consider ones where > 10% of predictions were NA (usually happens when there aren't enough boxes.)
model_avg_errors = all_model_errors %>%
  left_join(simulation_metadata, by='simulation_id') %>%
  group_by(model, model_id, sample_size, clustering, flowering_lengths, flowering_gradients, spatial_gradient_types) %>%
  summarise(rmse = mean(rmse, na.rm=T), 
            bias = mean(bias, na.rm=T),
            percent_na =  mean(percent_na),
            n=n()) %>%
  ungroup() 

model_avg_errors$model = factor(model_avg_errors$model, levels = c('linear','weibull_grid'),
                                                labels = c('Naive Model','Weibull Grid'))

##############################################################
##############################################################
# Primary simulation result figures.
#
# Errors for 1) underlying phenology, and 2) sampling scenarios
##############################################################
phenology_error_facet_labels = tribble(
  ~spatial_gradient_types, ~flowering_gradients, ~facet_label,~facet_order,
  'linear',16.8,     'A.      Weak Linear Gradient', 1, 
  'non-linear',16.8, 'B.    Weak Non-Linear Gradient', 2,
  'linear',33.6,     'C.   Moderate Linear Gradient', 3, 
  'non-linear',33.6, 'D. Moderate Non-Linear Gradient', 4, 
  'linear',67.2,     'E.      Strong Linear Gradient', 5, 
  'non-linear',67.2, 'F.    Strong Non-Linear Gradient', 6
)

# find the top model for all underlying phenology combinations for sample size of 300 and with clustering
phenology_error_means =  model_avg_errors %>%
  filter(sample_size == 300, clustering == T) %>%
  group_by(model, flowering_lengths, flowering_gradients, spatial_gradient_types) %>%
  top_n(1, -rmse) %>%
  ungroup() %>%
  group_by(model, model_id, flowering_lengths, flowering_gradients, spatial_gradient_types) %>%
  summarise(rmse = mean(rmse), n=n(), bias=mean(bias), percent_na = mean(percent_na)) %>%
  ungroup() %>%
  mutate(flowering_gradients = round(flowering_gradients,1)) %>% # round these so it joings correctly in the next line
  left_join(phenology_error_facet_labels, by=c('spatial_gradient_types','flowering_gradients')) %>%
  mutate(facet_label = forcats::fct_reorder(facet_label, facet_order))

phenology_error_model_labels = phenology_error_means %>%
  group_by(model, facet_label) %>%
  filter(flowering_lengths == max(flowering_lengths)) %>%
  ungroup() 

ggplot(phenology_error_means, aes(y=rmse, x=flowering_lengths, color=model)) +
  geom_line(size=2) +
  geom_point(size=3) +
  #scale_color_brewer(palette = 'Dark2') + 
  scale_color_manual(values = c('black','#0072B2','#E69F00')) +
  #scale_color_manual(values = viridis::magma(4, end=0.8)) + 
  #scale_color_viridis_d(end = 0.9) + 
  #scale_color_manual(values = c('black','grey20','grey60','grey80')) + 
  #scale_linetype_manual(values = c('solid','dashed','solid','dashed')) + 
  geom_text_repel(data = phenology_error_model_labels, aes(x=flowering_lengths + 0.5, label=stringr::str_wrap(model, 10)), 
                  size=4,fontface='bold', family='sans',
                  xlim=c(62,70), hjust=0, min.segment.length = 0.1) + 
  scale_x_continuous(breaks=c(15,30,45,60), labels = c(15,30,45,60), limits = c(15,70),
                     minor_breaks = c()) + 
  coord_cartesian(ylim=c(0,16)) + 
  facet_wrap(~facet_label, ncol=2) +
  theme_bw(base_size = 20) +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_rect(color='black', size=0.5),
        legend.key.width = unit(20,'mm'),
        axis.text = element_text(size=18),
        strip.text = element_text(size=20, hjust = 0),
        axis.title = element_text(size=20)) + 
  labs(x='Flowering Length (Days)', y='RMSE', color='', linetype='')

##############################################################

sampling_error_facet_labels = tribble(
  ~spatial_gradient_types, ~clustering, ~facet_label,~facet_order,
  'linear',FALSE,     'A.  Moderate Linear Gradient\n   Non-Clustered Sampling', 1, 
  'linear',TRUE,      'B.  Moderate Linear Gradient\n Clustered Sampling', 2, 
  'non-linear',FALSE, 'C.  Moderate Non-Linear Gradient\n Non-Clustered Sampling', 3,
  'non-linear',TRUE,  'D.  Moderate Non-Linear Gradient\n Clustered Sampling', 4
)

# find the top model for all sampling scenarios, but only within a flowering length of 30 and strength of 33.6
# still incorportating linear/nonlinear gradient here
sampling_error_means = model_avg_errors %>%
  filter(round(flowering_gradients,1) == 33.6, flowering_lengths == 30) %>%
  group_by(model, spatial_gradient_types, sample_size, clustering) %>%
  top_n(1, -rmse) %>%
  ungroup() %>%
  group_by(model,model_id, spatial_gradient_types, sample_size, clustering) %>%
  summarise(rmse = mean(rmse), bias=mean(bias), n=n()) %>%
  ungroup() %>%
  left_join(sampling_error_facet_labels, by=c('spatial_gradient_types','clustering')) %>%
  mutate(facet_label = forcats::fct_reorder(facet_label, facet_order))

sampling_error_model_labels = sampling_error_means %>%
  group_by(model, clustering) %>%
  filter(sample_size == max(sample_size)) %>%
  ungroup()

ggplot(sampling_error_means, aes(x=sample_size, y=rmse, color=model)) + 
  geom_line(size=2) +
  geom_point(size=3) +
  #scale_color_brewer(palette = 'Dark2') + 
  scale_color_manual(values = c('black','#0072B2','#E69F00')) +
  geom_text_repel(data = sampling_error_model_labels, aes(x=sample_size + 0.5, label=stringr::str_wrap(model, 10)),
                  size=5, fontface='bold', family='sans',
                  xlim=c(1200,1500), hjust=0, min.segment.length = 1) + 
  scale_x_continuous(breaks=c(150,300,600,1200), limits = c(150, 1500),
                     minor_breaks = c()) + 
  facet_wrap(~facet_label) +
  theme_bw(base_size = 20) + 
  theme(axis.text = element_text(size=18),
        strip.text = element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = 'none',
        legend.key.width = unit(20,'mm'),
        legend.background = element_rect(color='black'),
        legend.title = element_blank(),
        legend.justification = 'center') + 
  labs(x='Sample Size', y='RMSE', color='')

################################################
# Supplemental figures
#
# Same plots as above, but only show the weibull errors,
# and have labels showing the weibull grid parameters which generated those errors.
# Also a table showing *all* errors.
################################################

#####
# The best weibul grid parameters with regard to flowering length, gradient type & strength

underlying_phenology_best_model_parameters = model_avg_errors %>%
  filter(model == 'Weibull Grid') %>%
  filter(sample_size == 300, clustering == T) %>%
  group_by(model, flowering_lengths, flowering_gradients, spatial_gradient_types) %>%
  top_n(2, -rmse) %>%
  ungroup()  %>%
  left_join(weibull_model_parameters, by='model_id') %>%
  mutate(weibull_label = paste('Box size:',box_size,
                               '\n# Boxes:',num_boxes,
                               '\nGrid Size:',stratum_size,
                               '\nRMSE:',round(rmse,2),
                               '\nBias:',round(bias,2),
                               '\n% NA:',round(percent_na,2))) %>%
  mutate(label_y = ifelse(flowering_lengths%%10==0, 'on_top','on_bottom')) %>%
  mutate(flowering_gradients = round(flowering_gradients,1)) %>% # round these so it joings correctly in the next line
  left_join(phenology_error_facet_labels, by=c('spatial_gradient_types','flowering_gradients'))

ggplot(underlying_phenology_best_model_parameters,aes(y=rmse, x=flowering_lengths)) +
  geom_label_repel(data=filter(underlying_phenology_best_model_parameters, label_y=='on_top'),
                   aes(label=weibull_label), size=2.5,
                   direction = 'x', nudge_y=25) + 
  geom_label_repel(data=filter(underlying_phenology_best_model_parameters, label_y=='on_bottom'),
                   aes(label=weibull_label), size=2.5,
                   direction = 'x', nudge_y = -25) + 
  geom_point( size=3) + 
  scale_x_continuous(breaks=c(15,30,45,60), labels = c(15,30,45,60), limits=c(0,70),
                     minor_breaks = c()) + 
  ylim(-45,50) + 
  #coord_cartesian(ylim=c(0,16)) + 
  facet_wrap(~facet_label, ncol=2) +
  theme_bw(base_size = 20) +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_rect(color='black', size=0.5),
        legend.key.width = unit(20,'mm'),
        axis.text = element_text(size=18),
        strip.text = element_text(size=20, hjust = 0),
        axis.title = element_text(size=20)) + 
  labs(x='Flowering Length (Days)', y='RMSE', color='', linetype='')

############
# The best parameters with regard to sample size and clustering/no clustering
sampling_scenario_best_model_parameters = model_avg_errors %>%
  filter(model == 'Weibull Grid') %>%
  filter(round(flowering_gradients,1) == 33.6, flowering_lengths == 30) %>%
  group_by(model, spatial_gradient_types, sample_size, clustering) %>%
  top_n(2, -rmse) %>%
  ungroup() %>%
  left_join(weibull_model_parameters, by='model_id') %>%
  mutate(weibull_label = paste('Box size:',box_size,
                               '\n# Boxes:',num_boxes,
                               '\nGrid Size:',stratum_size,
                               '\nRMSE:',round(rmse,2),
                               '\nBias:',round(bias,2),
                               '\n% NA:',round(percent_na,2))) %>%
  left_join(sampling_error_facet_labels, by=c('spatial_gradient_types','clustering')) %>%
  mutate(facet_label = forcats::fct_reorder(facet_label, facet_order)) %>%
  mutate(label_y = ifelse(sample_size %in% c(150,600), 'on_top','on_bottom'))

ggplot(sampling_scenario_best_model_parameters,aes(x=sample_size, y=rmse, color=model)) + 
    geom_point(size=3) +
  geom_label_repel(data = filter(sampling_scenario_best_model_parameters, label_y=='on_top'),
                   aes(label = weibull_label),size=3, 
                   vjust=0, direction = 'x', nudge_y = 10) +
  geom_label_repel(data = filter(sampling_scenario_best_model_parameters, label_y=='on_bottom'),
                   aes(label = weibull_label),size=3, 
                   vjust=0, direction = 'x', nudge_y = -15) +
  scale_color_manual(values = c('black','#0072B2','#E69F00')) +
  scale_x_continuous(breaks=c(150,300,600,1200), limits = c(-100, 1500),
                     minor_breaks = c()) + 
  ylim(-10,15) + 
  facet_wrap(~facet_label) +
  theme_bw(base_size = 20) + 
  theme(axis.text = element_text(size=18),
        strip.text = element_text(size=20),
        axis.title = element_text(size=20),
        legend.position = 'none',
        legend.key.width = unit(20,'mm'),
        legend.background = element_rect(color='black'),
        legend.title = element_blank(),
        legend.justification = 'center') + 
  labs(x='Sample Size', y='Bias', color='')
  
##############################################################
# Supplemental table showing the best error for each sampling x phenology scenarios

library(kableExtra)

model_avg_errors %>%
  filter(sample_size == 300, clustering == T) %>%
  group_by(model, flowering_lengths, flowering_gradients, spatial_gradient_types) %>%
  top_n(1, -rmse) %>%
  ungroup()

best_model_error_for_all_simulations = model_avg_errors %>%
  select(-model_id, -bias, -percent_na) %>%
  group_by(model, sample_size, clustering, flowering_lengths, flowering_gradients, spatial_gradient_types) %>%
  top_n(1, -rmse) %>%
  ungroup()


error_table_data = best_model_error_for_all_simulations %>%
  spread(model, rmse) %>%
  rename(linear = 'Naive Model', weibull_grid = 'Weibull Grid') %>%
  mutate(error = paste0(round(weibull_grid,1), ' (',round(linear,1),')')) %>%
  mutate(weibull_is_better = weibull_grid < linear) %>%
  mutate(error = ifelse(weibull_is_better, paste0('\\textbf{',error,'}'), error)) %>% # add latex bold to entries where weibull < linear
  select(clustering, spatial_gradient_types, flowering_gradients, flowering_lengths, sample_size, error) %>%
  spread(sample_size, error) %>%
  arrange(clustering, spatial_gradient_types, flowering_gradients, flowering_lengths)

error_table_data %>%
  select(-clustering) %>%
  kable(format = 'latex', escape=F) %>%
  add_header_above(c(' ' = 3, 'Sample Size' = 4)) %>%
  pack_rows('Clustered Sampling', 1, 24) %>%
  pack_rows('Non-Clustered Sampling', 25,48)
