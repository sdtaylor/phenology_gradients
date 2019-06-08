context('Testing the flowering gradient simulator')

library(dplyr)
source('../../flowering_gradient_generator.R')

x = seq(0,1,0.2)
y = seq(0,1,0.2)

test_that('random gradient generators are reproducible given a seed', {
  expect_equal(round(generate_random_spatial_gradient(x,y,seed=5),7), c(0.2587186, 0.8261091, 0.9619575, 0.7931870, 0.6221046, 0.5546258))
  expect_equal(round(generate_random_spatial_gradient(x,y,seed=6),7), c(0.9789147, 0.6754984, 0.1602651, 0.2284302, 0.4881718, 0.4741535))
  expect_equal(round(generate_random_spatial_gradient(x,y,seed=6999),7), c(0.9482594, 0.8344935, 0.5773765, 0.2222792, 0.0051186, 0.1434193))
  
  true_onset_dates = spatialFloweringGrid(xlimits=c(0,1),
                                          ylimits = c(0,1),
                                          cell_size=0.4,
                                          start_doy=90,
                                          flowering_length = 30,
                                          flowering_gradient = 10/0.1,
                                          spatial_gradient_type = 'non-linear',
                                          seed=15)
  expect_equal(true_onset_dates$onset, c(90,  97, 116, 122, 119, 135, 104, 140, 167))
})

test_that('primary simulation function has correct output',{
  sample_size = 100
  simulated_data = spatialFloweringSampler(n=sample_size)
  
  # sample size specified is correct
  expect_equal(nrow(simulated_data), sample_size)
  # no na's
  expect_true(all(!is.na(simulated_data)))
  
  # should be flowering presence only
  simulated_data_po = spatialFloweringSampler(n=sample_size, sample_type = 'po')
  expect_true(all(simulated_data_po$flower_present==1))
  
  # should have flowering presence and absence, and in the proportion specified
  fraction_present = 0.5
  simulated_data_pa = spatialFloweringSampler(n=sample_size, 
                                              sample_type = 'pa',
                                              fraction_present = fraction_present)
  expect_equal(sort(unique(simulated_data_pa$flower_present)), c(0,1))
  expect_equal(mean(simulated_data_pa$flower_present), fraction_present)
  
  # DOY is always between 1-365
  expect_true(all(simulated_data$doy >= 1))
  expect_true(all(simulated_data$doy <= 365))
  expect_true(all(simulated_data_po$doy >= 1))
  expect_true(all(simulated_data_po$doy <= 365))
  # Need to fix this in pa generation, see https://github.com/sdtaylor/phenology_gradients/issues/3
  #expect_true(all(simulated_data_pa$doy >= 1))
  #expect_true(all(simulated_data_pa$doy <= 365))
})
