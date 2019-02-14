source('weibull.R')

PhenologyGridEstimator = function(doy_points,
                                  n_boxes=100,
                                  max_box_size=0.3,
                                  min_box_size=0.2,
                                  grain = 0.05,
                                  xlimits=c(0,1),
                                  ylimits=c(0,1)){
  # From a collection of georeferenced indicating the Julian day (DOY) of
  # observing a flower, make a model estimating the onset,peak, and end of 
  # flowering across the spatial extent.
  # Uses bag of boxes strategy from Fink et al. 2010 combined with the 
  # Weibull method from Pearse et al. 2017 to estimate onset and end. 
  # Uses the mean doy for peak estimates. 
  
  model_details = list()
  model_details$n_boxes = n_boxes
  model_details$max_box_size = max_box_size
  model_details$min_box_size = min_box_size
  model_details$grain = grain
  model_details$xlimits = xlimits
  model_details$ylimits = ylimits
  
  boxes = data.frame(box_id=1:n_boxes)
  boxes$size = runif(n_boxes, min=min_box_size, max=max_box_size)
  boxes$half_size = boxes$size/2
  boxes$center_x = runif(n=n_boxes, min=xlimits[1], max=xlimits[2])
  boxes$center_y = runif(n=n_boxes, min=ylimits[1], max=ylimits[2])
  
  estimate_metrics = function(box_id, size, half_size, center_x, center_y){
    doy_points_subset = subset_points_to_box(doy_points, 
                                             box_center_x = center_x,
                                             box_center_y = center_y,
                                             box_half_size = half_size)
    estimates = list()
    estimates$box_id = box_id
    estimates$n_points = nrow(doy_points)
    if(estimates$n_points<3){
      estimates$onset_estimate = NA
      estimates$end_estimate = NA
      estimates$peak_estimate = NA
    } else{
      estimates$onset_estimate = as.numeric(weib.limit(doy_points_subset$doy)[1])
      estimates$end_estimate = as.numeric(weib.limit(doy_points_subset$doy*-1)[1]) * -1
      estimates$peak_estimate = mean(doy_points_subset$doy)
    }
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
  
  estimate_metrics_from_model = function(x, y){
    box_subset = subset_boxes_to_point(x = x,
                                       y = y,
                                       boxes = model$boxes)
    
    estimates = list()
    estimates$x = x
    estimates$y = y
    estimates$onset_estimate = mean(box_subset$onset_estimate, na.rm=T)
    estimates$end_estimate = mean(box_subset$end_estimate, na.rm=T)
    estimates$peak_estimate = mean(box_subset$peak_estimate, na.rm=T)
    return(estimates)
  }
  
  point_estimates = purrr::pmap_df(doy_points[c('x','y')], estimate_metrics_from_model)
  
  return(dplyr::select(point_estimates, -x, -y))
}


plot.PhenologyGridEstimator = function(model){
  # Plot the boxes generated from the model
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
  
  boundary_box = dplyr::tribble(
    ~x, ~y, ~xend, ~yend,
    model$xlimits[1], model$ylimits[1], model$xlimits[2], model$ylimits[1], #bottom
    model$xlimits[2], model$ylimits[1], model$xlimits[2], model$ylimits[2], #right
    model$xlimits[2], model$ylimits[2], model$xlimits[1], model$ylimits[2], #top
    model$xlimits[1], model$ylimits[2], model$xlimits[1], model$ylimits[1]  #left
  )
  
  line_segments = line_segments_left %>%
    bind_rows(line_segments_top) %>%
    bind_rows(line_segments_right) %>%
    bind_rows(line_segments_bottom)
  
  ggplot() + 
    geom_segment(data=line_segments,aes(x=x, y=y, xend=xend, yend=yend)) +
    geom_segment(data=boundary_box, aes(x=x, y=y, xend=xend, yend=yend), color='red', size=2, linetype='dashed')
  
  
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


