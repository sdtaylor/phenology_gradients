library(tidyverse)
library(ggrepel)

source('config.R')

simulation_metadata = read_csv('results/models/2019-09-12-2208/simulation_metadata.csv')
all_model_errors = read_csv('results/models/2019-09-12-2208/all_errors.csv')

# For each simulation run, pick the models best performance from their best parameters used
best_model_errors = all_model_errors %>%
  group_by(model, simulation_id) %>%
  top_n(1, -rmse) %>%
  ungroup()

# attach simulation info, such as phenology length, sample size, etc. 
best_model_errors = best_model_errors %>%
  left_join(simulation_metadata, by='simulation_id')

best_model_errors$model = factor(best_model_errors$model, levels = c('linear','weibull_grid'),
                                                labels = c('Naive Model','Weibull Grid'))

##############################################################
##############################################################
# make average error plots centered around the 1) underlying phenology, 2) sampling
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

phenology_error_means =  best_model_errors %>%
  filter(sample_size == 300, clustering == T) %>%
  group_by(model, flowering_lengths, flowering_gradients, spatial_gradient_types) %>%
  summarise(rmse = mean(rmse), n=n()) %>%
  ungroup() %>%
  mutate(flowering_gradients = round(flowering_gradients,1)) %>%
  left_join(phenology_error_facet_labels, by=c('spatial_gradient_types','flowering_gradients')) %>%
  mutate(facet_label = forcats::fct_reorder(facet_label, facet_order))

phenology_error_model_labels = phenology_error_means %>%
  group_by(model, facet_label) %>%
  filter(flowering_lengths == max(flowering_lengths)) %>%
  ungroup() 

ggplot(phenology_error_means, aes(y=rmse, x=flowering_lengths, color=model)) +
  geom_line(size=2) +
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
  facet_wrap(~facet_label, scales='free', ncol=2) +
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
  'non-linear',FALSE, 'B.  Moderate Non-Linear Gradient\n Non-Clustered Sampling', 2,
  'linear',TRUE,      'C.  Moderate Linear Gradient\n Clustered Sampling', 3, 
  'non-linear',TRUE,  'D.  Moderate Non-Linear Gradient\n Clustered Sampling', 4
)

sampling_error_means = best_model_errors %>%
  filter(round(flowering_gradients,1) == 33.6, flowering_lengths == 30) %>%
  group_by(model, spatial_gradient_types, sample_size, clustering) %>%
  summarise(rmse = mean(rmse), n=n()) %>%
  ungroup() %>%
  left_join(sampling_error_facet_labels, by=c('spatial_gradient_types','clustering')) %>%
  mutate(facet_label = forcats::fct_reorder(facet_label, facet_order))

# sampling_error_means$clustering = factor(sampling_error_means$clustering , levels = c(TRUE, FALSE), labels = c('Clustered Sampling','Non-Clustered Sampling'))
# sampling_error_means$spatial_gradient_types = factor(sampling_error_means$spatial_gradient_types, levels = c('linear','non-linear'),
#                                                      labels = c('Moderate Linear Gradient','Moderate Non-Linear Gradient'))

sampling_error_model_labels = sampling_error_means %>%
  group_by(model, clustering) %>%
  filter(sample_size == max(sample_size)) %>%
  ungroup()

ggplot(sampling_error_means, aes(x=sample_size, y=rmse, color=model)) + 
  geom_line(size=3) +
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

#################################################
# Large supplement figure which has boxplots for *all* simulation combinations
library(kableExtra)
large_plot_errors = model_errors %>%
  group_by(model, flowering_lengths, flowering_gradients, spatial_gradient_types, clustering, sample_size) %>%
  summarise(rmse = mean(rmse)) %>%
  ungroup() 
  mutate(flowering_lengths_text = case_when(
    flowering_gradients == 16.8 ~ 'Weak Gradient (16.8 days)',
    flowering_gradients == 33.6 ~ 'Moderage Gradient (33.6 days)',
    flowering_gradients == 67.2 ~ 'Strong Gradient (67.2 days)')) 

# large_plot_errors = large_plot_errors %>%
#   filter(flowering_lengths %in% c(15,45),
#          sample_size %in% c(150,600,1200))


x = large_plot_errors %>%
  mutate(rmse = round(rmse, 2)) %>%
  spread(model, rmse) %>%
  select(spatial_gradient_types, flowering_gradients, flowering_lengths, clustering, sample_size, BiSquared, IDW, `Linear Model`,`Weibull Grid`) %>%
  arrange(spatial_gradient_types, flowering_gradients, flowering_lengths, clustering, sample_size, BiSquared, IDW, `Linear Model`,`Weibull Grid`)

kable(x, 'latex') %>%
  collapse_rows(columns = 1:4, valign = 'top')

ggplot(model_errors, aes(x=))

###################################
# Below is a method to make the plots using a random forest based error models, they end up
# equivalent to just taking the means of all errors like above. 
# model_errors$model = as.factor(model_errors$model)
# model_errors$spatial_gradient_types = as.factor(model_errors$spatial_gradient_types)
# model_errors$clustering = as.factor(model_errors$clustering)
# error_model = randomForest(rmse ~ model + sample_size + clustering + flowering_lengths + flowering_gradients + spatial_gradient_types ,
#                            data=model_errors)
# 
# pdp_phenology_data = partial(error_model, pred.var = c('model','flowering_lengths','flowering_gradients','spatial_gradient_types'), n.trees=500) %>%
#   as_tibble()
# pdp_phenology_facet_labels = tribble(
#   ~spatial_gradient_types, ~flowering_gradients, ~facet_label,~facet_order,
#   'linear',16.8, 'Weak Linear Gradient', 1, 
#   'non-linear',16.8, 'Weak Non-Linear Gradient', 2,
#   'linear',33.6, 'Moderate Linear Gradient', 3, 
#   'non-linear',33.6, 'Moderate Non-Linear Gradient', 4, 
#   'linear',67.2, 'Strong Linear Gradient', 5, 
#   'non-linear',67.2, 'Strong Non-Linear Gradient', 6
# )
# x = pdp_phenology_data %>%
#   mutate(flowering_gradients = round(flowering_gradients,1)) %>%
#   left_join(pdp_phenology_facet_labels, by=c('spatial_gradient_types','flowering_gradients')) %>%
#   mutate(facet_label = forcats::fct_reorder(facet_label, facet_order))