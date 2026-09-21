##################################################################
############ Prepare the environment for the analysis ############
##################################################################


######### Create output folders ----------
folder <- 'climate_human' # name of the mother folder

if (!dir.exists(paste0('output/', folder, '/'))) {
  dir.create(paste0('output/', folder, '/'))
  dir.create(paste0('output/', folder, '/models/'))
  dir.create(paste0('output/', folder, '/models/crossValidation/'))
  dir.create(paste0('output/', folder, '/models/evaluation/'))
  dir.create(paste0('output/', folder, '/models/rast/'))
  dir.create(paste0('output/', folder, '/modellingData/'))
  dir.create(paste0('output/', folder, '/modellingData/ModelRegion/'))
  dir.create(paste0('output/', folder, '/projections/'))
  dir.create(paste0('output/', folder, '/projections/rast'))
} else {
  print('output dirs already exist')
}

