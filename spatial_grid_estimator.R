source('weibull.R')
library(dplyr)
library(purrr)

PhenologyGridEstimator = function(doy_points,
                                  stratum_size_x=0.1,
                                  stratum_size_y=0.1,
                                  boxes_per_stratum=5,
                                  box_size=0.3,
                                  xlimits=c(0,1),
                                  ylimits=c(0,1),
                                  edge_buffer=0.1,
                                  not_enough_data_fallback='use_na'
                                  ){
  # From a collection of georeferenced  points indicating the Julian day (DOY) of
  # observing a flower, make a model estimating the onset,peak, and end of 
  # flowering across the spatial extent.
  # Uses bag of boxes strategy from Fink et al. 2010 combined with the 
  # Weibull method from Pearse et al. 2017 to estimate onset and end. 
  # Uses the mean doy for peak estimates. 
  
  if((xlimits[2]<xlimits[1]) | 
     (ylimits[2]<ylimits[1]) | 
     (length(xlimits) !=2  ) |
     (length(ylimits) !=2  )){stop(paste('xlimits and ylimits must each be of length 2, and have the lower',
                                         'and upper bounds in the 1st and 2nd nposition, respectively.'))}
  
  model_details = list()
  model_details$doy_points = doy_points
  model_details$stratum_size_x = stratum_size_x
  model_details$stratum_size_y = stratum_size_y
  model_details$boxes_per_stratum = boxes_per_stratum
  model_details$box_size = box_size
  model_details$xlimits = xlimits
  model_details$ylimits = ylimits
  model_details$edge_buffer = edge_buffer
  model_details$not_enough_data_fallback = not_enough_data_fallback
 
  generate_centers_within_stratum = function(min,max){
    runif(n = boxes_per_stratum, min, max)
  }
  
  # This data.frame describes the bottom left coordinate of the stratum cells.
  stratum_cells = expand.grid(x = seq(xlimits[1], xlimits[2] - stratum_size_x, by=stratum_size_x),
                              y = seq(ylimits[1], ylimits[2] - stratum_size_y, by=stratum_size_y))
  
  # uniform random centers, with each stratum cell getting an equal number.
  center_x = purrr::map2(stratum_cells$x, stratum_cells$x+stratum_size_x, generate_centers_within_stratum) %>%
    purrr::flatten_dbl()
  center_y = purrr::map2(stratum_cells$y, stratum_cells$y+stratum_size_y, generate_centers_within_stratum) %>%
    purrr::flatten_dbl()
  
  boxes = dplyr::tibble(center_x = center_x,
                        center_y = center_y)
  total_boxes = nrow(boxes)
  
  boxes$box_id = 1:total_boxes
  boxes$size = box_size
  boxes$half_size = boxes$size/2
  
  estimate_metrics = function(box_id, size, half_size, center_x, center_y){
    doy_points_subset = subset_points_to_box(doy_points, 
                                             box_center_x = center_x,
                                             box_center_y = center_y,
                                             box_half_size = half_size)
    estimates = list()
    estimates$box_id = box_id
    estimates$n_points = nrow(doy_points_subset)
    
    # If the weib.limit fails (due to low sample size), then revert to
    # using the first/last observed of the subset within a box
    if(not_enough_data_fallback=='first_observed'){
      fallback = function(doy){min(doy)}
    } else if(not_enough_data_fallback=='use_na'){
      fallback = function(doy){NA}
    } else {
      stop(paste('unknown option for not_enough_data_fallback: ',not_enough_data_fallback))
    }
    
    estimates$onset_estimate = tryCatch(as.numeric(weib.limit(doy_points_subset$doy)[1]),
                                        error = function(cond){fallback(doy_points_subset$doy)})
    estimates$end_estimate   = tryCatch(as.numeric(weib.limit(doy_points_subset$doy*-1)[1]) * -1,
                                        error = function(cond){fallback(doy_points_subset$doy * -1) * -1})
    estimates$peak_estimate  = mean(doy_points_subset$doy)
    
    return(estimates)
  }
  
  # Get estimates of all metrics for each box
  box_estimates = purrr::pmap_df(boxes, estimate_metrics)
  
  boxes = boxes %>%
    dplyr::left_join(box_estimates, by='box_id')
  
  model_details$boxes = boxes
  
  return(model_details)
  
}

predict.PhenologyGridEstimator = function(model,
                                          doy_points){
  
  outside_buffer = function(x,y){
    !within_bounds2(x,y,
                    x_low  = model$xlimits[1] + model$edge_buffer,
                    x_high = model$xlimits[2] - model$edge_buffer,
                    y_low  = model$ylimits[1] + model$edge_buffer,
                    y_high = model$ylimits[2] - model$edge_buffer)
  }
  
  estimate_metrics_from_model = function(x, y){
    box_subset = subset_boxes_to_point(x = x,
                                       y = y,
                                       boxes = model$boxes)
    
    estimates = list()
    estimates$x = x
    estimates$y = y
    if(outside_buffer(x,y)){
      estimates$onset_estimate = NA
      estimates$end_estimate = NA
      estimates$peak_estimate = NA
      estimates$outside_buffer = TRUE
    } else {
      estimates$onset_estimate = mean(box_subset$onset_estimate, na.rm=T)
      estimates$end_estimate = mean(box_subset$end_estimate, na.rm=T)
      estimates$peak_estimate = mean(box_subset$peak_estimate, na.rm=T)
      estimates$outside_buffer = FALSE
    }
    return(estimates)
  }
  
  point_estimates = purrr::pmap_df(doy_points[c('x','y')], estimate_metrics_from_model)
  
  outside_buffer_count = sum(point_estimates$outside_buffer)
  if(outside_buffer_count>0){
    warning(paste(outside_buffer_count,'points were outside the buffer and could not be estimated.'))
  }
  
  return(dplyr::select(point_estimates, -x, -y, -outside_buffer))
}


plot.PhenologyGridEstimator = function(model,
                                       plot_type='boxes'){
  # box_types: 
  #   boxes: draw all the boxes from the model on the x,y grid, 
  #          with a red dashed line marking the boundery
  #   density: plot the density (# of boxes) for all points
  # Plot the boxes generated from the model
  boundary_box = dplyr::tribble(
    ~x, ~y, ~xend, ~yend,
    model$xlimits[1], model$ylimits[1], model$xlimits[2], model$ylimits[1], #bottom
    model$xlimits[2], model$ylimits[1], model$xlimits[2], model$ylimits[2], #right
    model$xlimits[2], model$ylimits[2], model$xlimits[1], model$ylimits[2], #top
    model$xlimits[1], model$ylimits[2], model$xlimits[1], model$ylimits[1]  #left
  )
  
  if(plot_type=='boxes'){
    line_segments_left = with(model$boxes, data.frame(x = center_x - half_size,
                                                      y = center_y + half_size,
                                                      xend = center_x - half_size,
                                                      yend = center_y - half_size))
    line_segments_top = with(model$boxes, data.frame(x = center_x - half_size,
                                                      y = center_y + half_size,
                                                      xend = center_x + half_size,
                                                      yend = center_y + half_size))
    line_segments_right = with(model$boxes, data.frame(x = center_x + half_size,
                                                      y = center_y + half_size,
                                                      xend = center_x + half_size,
                                                      yend = center_y - half_size))
    line_segments_bottom = with(model$boxes, data.frame(x = center_x + half_size,
                                                      y = center_y - half_size,
                                                      xend = center_x - half_size,
                                                      yend = center_y - half_size))
    
    line_segments = line_segments_left %>%
      bind_rows(line_segments_top) %>%
      bind_rows(line_segments_right) %>%
      bind_rows(line_segments_bottom)
    
    ggplot() + 
      geom_segment(data=line_segments,aes(x=x, y=y, xend=xend, yend=yend)) +
      geom_segment(data=boundary_box, aes(x=x, y=y, xend=xend, yend=yend), color='red', size=2, linetype='dashed')
  
  } else if(plot_type=='density'){
    # An evenly spaced grid of 2500 points on the model's x,y plane
    equal_spaced_grid = expand.grid(x=seq(model$xlimits[1],model$xlimits[2],by=(model$xlimits[2]-model$xlimits[1])/50),
                                    y=seq(model$ylimits[1],model$ylimits[2],by=(model$xlimits[2]-model$xlimits[1])/50))
    
    # count the boxes behind each point
    box_count = function(x,y){
      return(nrow(subset_boxes_to_point(x=x,y=y,boxes=model$boxes)))
    }
    
    equal_spaced_grid = equal_spaced_grid %>%
      mutate(n_boxes = map2_int(x,y,box_count))
    
    ggplot(equal_spaced_grid, aes(x=x,y=y,fill=n_boxes)) + 
      geom_raster() +
      scale_fill_viridis_c(limits=c(min(equal_spaced_grid$n_boxes),max(equal_spaced_grid$n_boxes))) + 
      geom_segment(data=boundary_box, aes(x=x, y=y, xend=xend, yend=yend), color='red', size=2, linetype='dashed', inherit.aes = FALSE)
    
  } else {
    stop(paste0('unknown plot type: ',plot_type))
  }
  
}

#######################################
#######################################
# Helper Function
#######################################
#######################################
subset_points_to_box = function(points, box_center_x, box_center_y, box_half_size){
  #
  # Given a single box center (box_x, box_y) and it's half size, return points which 
  # are conatined in ti.
  #
  points %>%
    filter(x > box_center_x - box_half_size,
           x <= box_center_x + box_half_size,
           y > box_center_y - box_half_size,
           y <= box_center_y + box_half_size)
  
}

subset_boxes_to_point = function(x,y, boxes){
  #
  # Given a single point (x,y), return boxes which that point is in
  #
  boxes %>%
    filter(center_x - half_size < x,
           center_x + half_size >= x,
           center_y - half_size < y,
           center_y + half_size >= y)
}

within_bounds = function(x, low, high){
  x >= low & x <= high
}

within_bounds2 = function(x, y, 
                          x_low, x_high,
                          y_low, y_high){
  within_bounds(x, x_low, x_high) & within_bounds(y, y_low, y_high)
}
