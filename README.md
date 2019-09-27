## Opportunistically collected photographs can be used to estimate large-scale phenological trends

This is the code and primary data for the above study. 

The analysis develops and tests a method, the Weibull Grid, to estimate flowering onset across large-scales using only observations of flowering presence. It uses a Weibull estimator (from the [phest](https://github.com/willpearse/phest) package) integrated into a spatially explicit ensemble. The method is tested with a simulation study and then applied to actually phenology data derived from iNaturalist. 

Code for the Weibull Grid model and making flowering simulations is in the [flowergrids](https://github.com/sdtaylor/flowergrids) R package. 

Data sources:

[iNaturalist](https://www.inaturalist.org)  
[Phenocam](https://phenocam.sr.unh.edu/webcam/)  
[USA National Phenology Network*](https://www.usanpn.org)  
[iDigBio*](https://www.idigbio.org)  
*These data sources were only used to show historic trends in data collection, not in the main analysis.
