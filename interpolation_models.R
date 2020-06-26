




InterpolationEstimator = function(doy_points,
                                  xlimits=c(0,1),
                                  ylimits=c(0,1),
                                  resolution = 0.1,
                                  percentile = 0.99,
                                  method = 'idw',
                                  power = 2,
                                  bisq_distance = 0.5){
  model_details = list()
  model_details$doy_points = doy_points
  model_details$xlimits = xlimits
  model_details$ylimits = ylimits
  model_details$percentile = percentile
  model_details$resolution = resolution
  model_details$method = method
  model_details$power = power
  model_details$bisq_distance = bisq_distance
  
  # TODO: do some checks here.
  # all mathy stuff is done in the prediction function
  return(model_details)
}

.euclidian_distance = function(single_point, many_points){
  # single_point = c(x,y), 
  # many_points = matrix(x,y)
  # returns a 1d vector of distances to all of many_points from single_point
  sqrt( (many_points[,1] - single_point[1])**2 + (many_points[,2] - single_point[2])**2 )
}

predict.InterpolationEstimator = function(model,
                                          doy_points=NULL,
                                          type='onset'){
  # https://andrewpwheeler.wordpress.com/2016/11/02/some-inverse-distance-weighting-hacks-using-r-and-spatstat/
  # Predict on the fitting data
  if(is.null(doy_points)){
    doy_points = model$doy_points
  }
  
  if(model$method == 'idw'){
    weight_transform = function(x){1/(x**model$power)}
  } else if(model$method == 'bisq'){
    weight_transform = function(x){ifelse( x < model$bisq_distance, (1 - (x/model$bisq_distance)**2), 0)}
  }
  
  fitting_points = cbind(model$doy_points$x, model$doy_points$y)
  new_points = cbind(doy_points$x, doy_points$y)
  prediction.mean = rep(NA, nrow(new_points))
  prediction.var  = rep(NA, nrow(new_points))
  for(i in 1:nrow(new_points)){
    distances = .euclidian_distance(new_points[i,], fitting_points)
    distances = weight_transform(distances)
    if(all(distances==0)){
      # No nearby points with this bisq_distance parameter
      next
    }
    prediction.mean[i] = weighted.mean(model$doy_points$doy, w = distances)
    prediction.var[i] = spatstat::weighted.var(model$doy_points$doy, w = distances)
  }
  
  
  lower_p = (1 - model$percentile)/2
  upper_p = 1 - lower_p
  
  if(type == 'onset'){
    estimate = qnorm(lower_p, mean = prediction.mean, sd = sqrt(prediction.var))
  } else if(type == 'end'){
    estimate = qnorm(upper_p, mean = prediction.mean, sd = sqrt(prediction.var))
  } else if(type == 'peak'){
    estimate = prediction.mean
  } else {
    stop(paste('unknown prediction type: ',type))
  }
  
  return(estimate)
}