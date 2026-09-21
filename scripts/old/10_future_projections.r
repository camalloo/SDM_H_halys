library("terra")
library("maxnet")
library("mgcv")
library("embarcadero")
library("fuzzySim")


################# Load data_________________________________
folder <- 'world' # name of the mother folder

data <- read.csv(paste0('output/', folder, '/modellingData/modData.csv')) # train dataset
vars <- read.csv(paste0('output/', folder, '/modellingData/selectedVariables.csv'))[, 1] # variables after collinearity analysis
world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antarctica
region <- vect(paste0('output/', folder, '/modellingData/ModelRegion/modellingRegionPoly.shp')) # modelling region polygon
AUC <- unlist(read.csv(paste0('output/', folder, '/models/evaluation/AUCused.csv'))) # for AUC based weighted mean ensemble

## load the models
mxtMod <- readRDS(paste0('output/', folder, '/models/mxt_mod.rds'))
gamMod <- readRDS(paste0('output/', folder, '/models/gam_mod.rds'))
## BART MODEL cannot be saved so always re-fit it if not loaded on env (KEEP THE SEED)
bartMod <- dbarts::bart(x.train = data[, vars], y.train = data[, 'presence'],
                        ntree = 100, keeptrees = T, verbose = F, seed = 1312)

## set up the counters for the loops
outerPath <- 'data/CHELSAbioclim/Future'
periods <- c('2011_2040', '2041_2070')
scenarios <- c('ssp126', 'ssp370', 'ssp585')





################# Download future scenarios bioclims and process them_________________________________
######## download only bios used to train the model --------
# ## for each of the five climate models, in the three ssp scenarios, of the two time periods, download the selected variables
# source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/refs/heads/master/download_vars")
# ## loop over the two time periods
# for (period in periods) {
#   ## second loop over the three scenarios
#   for (s in scenarios) {
#     models <- list.files(paste0(outerPath, '/', period, '/', s), pattern = '\\.txt$', full.names = T)
#     ## third loop over the five models
#     for (m in models) {
#       print(paste('Working on file:', stringr::str_extract(m, '[^/]+$')))
#       ## select the correct folder in which to save the downloaded files (one folder for each climate model)
#       fld <-  file.path(outerPath, period, s, stringr::str_extract(m, paste0('(?<=', period, '_).*(?=_', s, ')')))
#       print(paste('Saving to folder', fld))
# 
#       cat('\n')
#       ## read the .txt containing the link of the downloadable files
#       downloadable <- read.table(m, header = F)[,1]
#       ## select only the bios that were used to train the model to save memory and time
#       toDownload <- grep(paste0('_', paste(vars, collapse = '_|_'), '_'), downloadable, value = T)
#       ## download the selected bios and save them to the correct folder
#       downloadVars(toDownload, fld)
#       cat('\n')
#       cat('\n')
#       cat('\n')
#     }
#   }
# }


######## crop the bios to include only land, save them, and plot them --------
## create a directory to store the processed bios for modelling
if (!dir.exists(paste0('output/', folder, '/modellingData/Future'))) {
  dir.create(paste0('output/', folder, '/modellingData/Future'))
} else {
  cat(paste0('Future dir already exists in:\n', '\toutput/', folder, '/modellingData/Future/'))
}

## loop over all future data, process them and save them for modelling
## loop over the two time periods
for (period in periods) {
  cat(paste0('Starting period ', period, '...\n\n'))
  ## create a folder to store the processed future bio period
  if (!dir.exists(paste0('output/', folder, '/modellingData/Future/', period))) {
    dir.create(paste0('output/', folder, '/modellingData/Future/', period))
  }
  ## second loop over the ssp scenarios
  for (s in scenarios) {
    cat(paste0('Starting scenario ', s, '...\n'))
    modelsFlds <- list.dirs(file.path(outerPath, period, s), recursive = F, full.names = T)
    ## third loop over the five models
    for (m in modelsFlds) {
      print(basename(m))
      bio <- rast(list.files(m, pattern = '\\.tif$', full.names = T))

      ## fix the names of the layers
      names(bio) <- gsub('.*(bio[0-9]+).*', '\\1', names(bio))
      ## check if raster names match with variables used for training models
      if (!all(names(bio) == vars)) {
        cat('Raster names don\'t match with training variables\n...aborting')
        break
      }

      ## create a folder to store the processed future bio ssp
      if (!dir.exists(paste0('output/', folder, '/modellingData/Future/', period, '/', s))) {
        dir.create(paste0('output/', folder, '/modellingData/Future/', period, '/', s))
      }

      ## aggregate to match the resolution of the training dataset
      bio <- aggregate(bio, fact = 10, fun = 'mean', na.rm = T, cores = 6)
      
      ## crop bioclim to include only land
      bio <- crop(bio, world, mask = T)

      ## save the bios in the folder
      writeRaster(bio, filename = paste0('output/', folder, '/modellingData/Future/', period, '/', s, '/', basename(m), '.tif'),
                  gdal = c('COMPRESS=DEFLATE'), filetype = 'GTiff', overwrite = T)
      ## plot them to check if everything correct
      plot(bio, nc = 2, nr = 3)
      title(main = paste(basename(m), s, period), outer = T, line = -1)
      cat('\n')
    }
    cat(paste('...Scenario', s, 'done\n'))
  }
  cat(paste('...Period', period, 'completed\n\n\n'))
}





################# Run the simulations_________________________________
## IT TAKES 6:30 HOURS!!
## create a directory to store the outputs of the model projections
if (!dir.exists(paste0('output/', folder, '/projections/Future'))) {
  dir.create(paste0('output/', folder, '/projections/Future'))
} else {
  cat(paste0('Future dir already exists in:\n', '\toutput/', folder, '/projections/Future/'))
}

## time tracking
(start <- Sys.time())

## loop over all future data and run the models
## loop over the two time periods
for (period in periods) {
  cat(paste0('Starting period ', period, '...\n\n'))
  ## create a folder to store the projection periods
  if (!dir.exists(paste0('output/', folder, '/projections/Future/', period))) {
    dir.create(paste0('output/', folder, '/projections/Future/', period))
  }
  
  ## second loop over the ssp scenarios
  for (s in scenarios) {
    cat(paste0('Starting scenario ', s, '...\n'))
    ## create a folder to store the projection ssps
    if (!dir.exists(paste0('output/', folder, '/projections/Future/', period, '/', s))) {
      dir.create(paste0('output/', folder, '/projections/Future/', period, '/', s))
    }
    climateModels <- list.files(paste0('output/', folder, '/modellingData/Future/', period, '/', s), pattern = '.tif')
    
    ## third loop over the five models
    for (m in climateModels) {
      ## get bioclim data for climate model m, for ssp scenario s, at time period
      cat(m, '\n')
      bio <- rast(paste0('output/', folder, '/modellingData/Future/', period, '/', s, '/', basename(m)))
      
      ## run projections
      ## maxent
      mxtFutProj <- predict(bio, mxtMod, type = 'cloglog', na.rm = T, cores = 4, cpkgs = 'maxnet')
      names(mxtFutProj) <- paste(sub('.tif', '', m), 'MXT.pred', sep = '_')
      cat('\tmaxent done\n')
      # ## plot
      # plot(mxtFutProj, main = paste('MXT', sub('.tif', '', m), s, period) , range = c(0, 1))

      ## gam
      gamFutProj <- predict(bio, gamMod, type = 'response')
      names(gamFutProj) <- paste(sub('.tif', '', m), 'GAM.pred', sep = '_')
      gamFutFav <- Fav(pred = gamFutProj, sample.preval = prevalence(model = gamMod))
      names(gamFutFav) <- paste(sub('.tif', '', m), 'GAM.fav', sep = '_')
      cat('\tgam done\n')
      # ## plot
      # par(mfrow = c(2, 1))
      # plot(gamFutProj, main = paste('GAM', sub('.tif', '', m), s, period), range = c(0, 1))
      # plot(gamFutFav, main = paste('GAM Fav', sub('.tif', '', m), s, period), range = c(0, 1))

      ## bart
      bartFutProj <- predict2.bart(bartMod, raster::stack(bio), quantiles = c(.05, .95), splitby = 5, quiet = T)
      names(bartFutProj) <- paste(sub('.tif', '', m), c('BART.pred', 'BART.lower.int', 'BART.upper.int'), sep = '_')
      bartFutProj <- rast(bartFutProj)
      crs(bartFutProj) <- 'epsg:4326'
      ## convert to prevalence-independent favourability
      bartFutFav <- Fav(pred = bartFutProj, sample.preval = prevalence(model = bartMod))
      names(bartFutFav) <- paste(sub('.tif', '', m), c('BART.fav', 'BART.fav.lower.int', 'BART.fav.upper.int'), sep = '_')
      ## add the uncertainty range
      bartFutProj$uncertainty.range <- bartFutProj[[3]] - bartFutProj[[2]] # to pred
      bartFutFav$uncertainty.range <- bartFutFav[[3]] - bartFutFav[[2]] # to fav
      cat('\tbart done\n')
      # ## plot
      # par(mfrow = c(1, 1))
      # plot(bartFutProj, range = c(0, 1))
      # title(main = paste('BART', sub('.tif', '', m), s, period), outer = T, line = -3)
      # plot(bartFutFav, range = c(0, 1))
      # title(main = paste('BART Fav', sub('.tif', '', m), s, period), outer = T, line = -3)

      ## combine projections in a single raster
      ## check crs
      if (!(same.crs(mxtFutProj, gamFutFav) && same.crs(gamFutFav, bartFutFav))) {
        cat('coordinates don\'t match\n...aborting')
        break
      }
      ## combine
      allFutProj <- c(mxtFutProj, gamFutFav, bartFutFav[[1]])
      allFutProj
      names(allFutProj) <- sub('.pred|.fav', '', names(allFutProj))
      plot(allFutProj, range = c(0, 1))
      title(main = paste('Future projections', sub('.tif', '', m), s, period), outer = T, line = -3)

      ## save to disk
      writeRaster(mxtFutProj, filename = paste0('output/', folder, '/projections/Future/', period, '/', s, '/',
                                                'mxt_fut_proj', '_', sub('.tif', '', m), '.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
      writeRaster(gamFutProj, filename = paste0('output/', folder, '/projections/Future/', period, '/', s, '/',
                                                'gam_fut_proj', '_', sub('.tif', '', m), '.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
      writeRaster(gamFutFav, filename = paste0('output/', folder, '/projections/Future/', period, '/', s, '/',
                                               'gam_fut_fav', '_', sub('.tif', '', m), '.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
      writeRaster(bartFutProj, filename = paste0('output/', folder, '/projections/Future/', period, '/', s, '/',
                                                 'bart_fut_proj', '_', sub('.tif', '', m), '.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
      writeRaster(bartFutFav, filename = paste0('output/', folder, '/projections/Future/', period, '/', s, '/',
                                                'bart_fut_fav', '_', sub('.tif', '', m), '.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
      ## all combined (for gam and bart favorability is used)
      writeRaster(allFutProj, filename = paste0('output/', folder, '/projections/Future/', period, '/', s, '/',
                                                'Fut_proj', '_', sub('.tif', '', m), '.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
    }
    cat(paste('...Scenario', s, 'done\n\n'))
  }
  cat(paste('...Period', period, 'completed\n\n\n'))
}

## time tracking
(end <- Sys.time())
(length <- difftime(end, start, units = 'mins'))
system('say Done')





################# Compute ensembles_________________________________
## create a directory to store the ensembles
if (!dir.exists(paste0('output/', folder, '/projections/Future/ensembles'))) {
  dir.create(paste0('output/', folder, '/projections/Future/ensembles'))
} else {
  cat(paste0('Future dir already exists in:\n', '\toutput/', folder, '/projections/Future/ensembles/'))
}


######## compute ensemble of the five climate model --------
## loop over all projections and compute the ensemble of the climate models for each SDM model
## loop over the two time periods
for (period in periods) {
  cat(paste0('Starting period ', period, '...\n\n'))
  
  ## second loop over the ssp scenarios
  for (s in scenarios) {
    cat(paste0('Starting scenario ', s, '...\n'))
    
    ## maxent projections of the climate models in a single layer
    mxtProjs <- rast(list.files(paste0('output/', folder, '/projections/Future/', period, '/', s), pattern = 'mxt_fut_proj', full.names = T))
    mxtEns <- mean(mxtProjs)
    mxtEns$variance <- app(mxtProjs, fun = var)
    plot(mxtEns, nc = 1, nr = 2)
    title(main = paste('Ensemble of climate models MXT', s, period), outer = T, line = -1)
    writeRaster(mxtEns, filename = paste0('output/', folder, '/projections/Future/ensembles/', period, '_', s, '_mxtFutEns.tif'),
                gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
    
    ## gam favourability of the climate models in a single layer
    gamProjs <- rast(list.files(paste0('output/', folder, '/projections/Future/', period, '/', s), pattern = 'gam_fut_fav', full.names = T))
    gamEns <- mean(gamProjs)
    gamEns$variance <- app(gamProjs, fun = var)
    plot(gamEns, nc = 1, nr = 2)
    title(main = paste('Ensemble of climate models GAM', s, period), outer = T, line = -1)
    writeRaster(gamEns, filename = paste0('output/', folder, '/projections/Future/ensembles/', period, '_', s, '_gamFutEns.tif'),
                gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
    
    ## bart favourability (exclude lower and upper intervals and uncertainty) of the climate models in a single raster
    bartProjs <- rast(list.files(paste0('output/', folder, '/projections/Future/', period, '/', s), pattern = 'bart_fut_fav', full.names = T), lyrs = 1)
    bartEns <- mean(bartProjs)
    bartEns$variance <- app(bartProjs, fun = var)
    plot(bartEns, nc = 1, nr = 2)
    title(main = paste('Ensemble of climate models BART', s, period), outer = T, line = -1)
    writeRaster(bartEns, filename = paste0('output/', folder, '/projections/Future/ensembles/', period, '_', s, '_bartFutEns.tif'),
                gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
    
    cat(paste('...Scenario', s, 'done\n\n'))
  }
  cat(paste('...Period', period, 'completed\n\n\n'))
}

## add the ensembles for each ssp in a single map for each model in the two time periods and plot them
for (period in periods) {
  for (model in c('mxt', 'gam', 'bart')) {
    maps <- list.files(path = paste0('output/', folder, '/projections/Future/ensembles'), 
                       pattern = paste0('^', period, '_ssp[0-9]{3}_', model, 'FutEns\\.tif'), full.names = T)
    ssps <- rast(maps)
    
    ## plot them
    par(mfrow = c(3,2), oma=c(1,2,3,1))
    at = 5/6
    for (lyr in 1:nlyr(ssps)) {
      plot(ssps[[lyr]], main = names(ssps[[lyr]]), range = c(0,1))
      if (lyr %% 2 == 0) {
        mtext(sub(".*(ssp[0-9]{3}).*", "\\1", basename(sources(ssps[[lyr]]))), side = 2, outer = T, at = at)
        at = at - 2/6
      }
    }
    title(main = paste(period, model), outer = T, cex.main = 2, line = 0)
    
    ## save to disk
    writeRaster(ssps, filename = paste0('output/', folder, '/projections/Future/ensembles/', period, '_all_ssps_', model, 'FutEns.tif'),
                gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  }
}


######## compute ensemble of the 3 sdm models --------
for (period in periods) {
  for (scenario in scenarios){
    maps <- list.files(path = paste0('output/', folder, '/projections/Future/ensembles'), 
                       pattern = paste0(period, '_', scenario, '_'), full.names = T)
    for (m in maps) {
      cat(m, '\n')
    }
    ens <- rast(maps)
    mods <- ens[[c(1,3,5)]]
    ens <- weighted.mean(mods, w = AUC)
    names(ens) <- 'w.mean'
    ens$variance <- app(mods, var)
    writeRaster(ens, filename = paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/', period, '_', scenario, 'Ens.tif'),
                gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
    cat('\n\n')
  }
}





################# Clean memory_________________________________
rm(list = ls())
gc()


