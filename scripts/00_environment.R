##################################################################
############ Prepare the environment for the analysis ############
##################################################################


######### Create output folders ----------
folder <- 'baseline_models_yearly_runs' # name of the mother folder

if (!dir.exists(paste0('output/', folder, '/'))) {
  dir.create(paste0('output/', folder, '/'))
  dir.create(paste0('output/', folder, '/models/'))
  dir.create(paste0('output/', folder, '/models/crossValidation/'))
  dir.create(paste0('output/', folder, '/models/evaluation/'))
  dir.create(paste0('output/', folder, '/models/rast/'))
  dir.create(paste0('output/', folder, '/modellingData/'))
  dir.create(paste0('output/', folder, '/modellingData/ModelRegion/'))
  dir.create(paste0('output/', folder, '/mess/'))
  dir.create(paste0('output/', folder, '/niche/'))
  dir.create(paste0('output/', folder, '/ensemble/'))
} else {
  print('output dirs already exist')
}

