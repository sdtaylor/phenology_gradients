

#reqs 
# 3 different dist types, normal, beta, uniform  
# presence only, or presence and absence
# 2 different spatial sampling types, random or clustered

spatialFloweringSampler = function(n,
                                   x=NULL,
                                   y=NULL,
                                   sample_type='po',
                                   # po: presence only, pa: presence and absence
                                   fraction_precent=0.8,
                                   distribution_type='a',
                                   # a: normal-ish curve, b: beta-ish curve, c: uniform-ish curve
                                   xlimits = c(0,1),
                                   ylimits = c(0,1),
                                   start_doy = 180,
                                   flowering_length = 30,
                                   flowering_gradient = 10/0.1){
  
  # Generate n points from a spatial domain of flowering
  
  # generate n random points on the x,y grid
  if(is.null(x) && is.null(y)){
    x=runif(n=n, min=xlimits[1], max=xlimits[2])
    y=runif(n=n, min=ylimits[1], max=ylimits[2])
  } else {
    if(length(x)!=length(y)){stop('x and y must have the same length')}
    if(length(x)!=n){stop('x and y must be the same as n')}
  }

  # make weights for each julian day corresponding to the chance of observing an open flower
  # corresponds to flower abundance
  doy_probabilites = rep(0, 365)
  if(distribution_type=='a'){
    doy_probabilites[start_doy:(start_doy+flowering_length-1)] = dnorm(seq(from=-3,to=3, length.out = flowering_length))
  } else if(distribution_type=='b'){
    doy_probabilites[start_doy:(start_doy+flowering_length-1)] = dbeta(seq(from=0,to=1, length.out = flowering_length), shape1=2, shape2=5)
  } else if(distribution_type=='c'){
    stop('type c not implemented yet')
  }
  
  doy = sample(1:365)
  
  # presences and absences included
  if(sample_type=='pa'){
    num_presence = ceiling(n*fraction_precent)
    num_absence  = n - num_presence
    
    doy =  sample(1:365, size=num_presence, replace=T, prob = doy_probabilites)
    doy =  c(doy,sample(1:365, size=num_absence, replace=T, prob = max(doy_probabilites) - doy_probabilites))
    flower_present = c(rep.int(1,num_presence),rep.int(0,num_absence))
  # only presences
  } else if(sample_type=='po'){
    doy = sample(1:365, size=n, replace=T, prob = doy_probabilites)
    flower_present = rep.int(1,n)
  } else {
    stop(paste0('unknown sample type: ',sample_type))
  }
  
  # for each point, transform the probability according to the flowering gradient
  #y_center = ylimits[2] - ((ylimits[2] - ylimits[1])/2)
  #y_scaled = y - y_center
  
  doy = doy + round(y * flowering_gradient, 0)
  
  true_start_doy = rep.int(start_doy, n)
  true_start_doy = true_start_doy + round(y * flowering_gradient, 0)
  
  return(data_frame(x=x,y=y,doy=doy,flower_present=flower_present, true_start_doy=true_start_doy))
}
