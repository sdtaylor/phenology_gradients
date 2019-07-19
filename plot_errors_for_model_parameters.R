library(tidyverse)
library(ggrepel)
#source('flowering_gradient_generator.R')
#source('spatial_grid_estimator.R')
source('config.R')


simulation_metadata = read_csv('results/models/2019-07-05-0935/simulation_metadata.csv') %>%
  mutate(flowering_gradients = round(flowering_gradients, 1)) 
model_errors = read_csv('results/models/2019-07-05-0935/all_errors.csv')

model_errors = model_errors %>%
  left_join(simulation_metadata, by='simulation_id')
  


########################################################
# idw model
idw_errors = model_errors %>%
  filter(model == 'idw') %>%
  left_join(idw_model_parameters, by='model_id')

idw_errors %>%
  filter(sample_size %in% c(150,300), 
         flowering_lengths %in% c(15,30, 60)) %>%
ggplot(aes(x=interaction(clustering, sample_size), y=rmse, fill=as.factor(idw_power))) + 
  geom_boxplot() + 
  scale_fill_viridis_d(end=0.8) + 
  facet_wrap(flowering_gradients ~ flowering_lengths ~ spatial_gradient_types, labeller = label_both, scales = 'free', ncol=6)

#########################################################
# bi-squared model
bisq_errors = model_errors %>%
  filter(model == 'bisq') %>%
  left_join(bisq_model_parameters, by='model_id')

bisq_errors %>%
  filter(sample_size %in% c(150,300), 
         flowering_lengths %in% c(15,30, 60)) %>%
ggplot(aes(x=interaction(clustering, sample_size), y=rmse, fill=as.factor(b_distance))) + 
  geom_boxplot() + 
  scale_fill_viridis_d(end=0.8) + 
  facet_wrap(flowering_gradients ~ flowering_lengths ~ spatial_gradient_types, labeller = label_both, scales = 'free', ncol=6)

#########################################################
######################################
library(cowplot)
weibull_errors = model_errors %>%
  filter(model == 'weibull_grid') %>%
  left_join(weibull_model_parameters, by='model_id')

# Sweet gridded legend for the weibull grid model, which as 2 primary parameters
unique_box_sizes = unique(weibull_errors$box_size)
unique_buffer_sizes =  unique(weibull_errors$edge_buffer)
legend_data = expand.grid(box_size = as.factor(unique_box_sizes), edge_buffer = as.factor(unique_buffer_sizes)) 
legend_data$id = paste('id-',1:nrow(legend_data),sep = '')

gridded_legend = ggplot(legend_data, aes(x=box_size, y=fct_rev(edge_buffer), fill=id)) + 
  geom_raster() +
  #scale_fill_brewer(palette = 'Paired') +
  scale_fill_manual(values=c('grey60','grey10','#edca7b','#E69F00','#a6cee3','#1f78b4')) + 
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(size=10),
        axis.title = element_text(size=11),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color='NA'),
        plot.background = element_rect(color='black', size=1.5)) +
  labs(x='Box Size',y='Buffer Size')

weibull_error_plot = weibull_errors %>%
  filter(sample_size %in% c(150,300), 
         flowering_lengths %in% c(15,30, 60)) %>%
ggplot(aes(x=interaction(clustering, sample_size), y=rmse, fill=interaction(box_size, edge_buffer))) + 
  geom_boxplot() + 
  # Uncomment below to put a label on the bars and ensure the custom legend colors match
  #geom_label(aes(label =interaction(box_size, edge_buffer)), position = position_dodge(width=0.5)) + 
  scale_fill_manual(values=c('grey60','grey10','#edca7b','#E69F00','#a6cee3','#1f78b4')) + 
  facet_wrap(flowering_gradients ~ flowering_lengths ~ spatial_gradient_types, labeller = label_both, scales = 'free', ncol=6) +
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18),
        strip.text = element_text(size=18),
        legend.position = 'none')

ggdraw() + 
  draw_plot(weibull_error_plot, 0,0,1,1) + 
  draw_plot(gridded_legend, x=0.85, y=0.1, width=0.14, height = 0.18)

