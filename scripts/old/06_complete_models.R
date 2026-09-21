library("terra")
library("predicts")
library("modEvA")
library("maxnet")
library("fuzzySim")
library("dbarts")
library("embarcadero")
library("randomForest")
library("mgcv")


################# Load data_________________________________
folder <- 'AllPredictors' # name of the mother folder

world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica
bioCropped <- rast(paste0('output/', folder, '/modellingData/ModelRegion/predsCropped.tif'))
data <- read.csv(paste0('output/', folder, '/modellingData/modData.csv')) # train dataset
selectedMods <- read.csv(paste0('output/', folder, '/models/crossValidation/selectedModels.csv'))
selectedMods
vars <- read.csv(paste0('output/', folder, '/modellingData/selectedVariables.csv'))[ , 1] # variables after collinearity analysis
region <- vect(paste0('output/', folder, '/modellingData/ModelRegion/modellingRegionPoly.shp')) # modelling region polygon





################# Run complete models_________________________________
## List of models used 'BC'   'MXT'   'GLM'   'GAM'   'BART'    'RF'
dataSPAT <- vect(data, geom = c('x', 'y'), crs = 'EPSG:4326')
## time tracking
(start <- Sys.time())

######## Bioclim - PRESENCE-ONLY (envelope) --------
if ('BC' %in% selectedMods$x) {
  bcMod <- predicts::envelope(x = subset(data, presence == 1, select = vars))
  plot(bcMod, a = 2, b = 7, p = 0.5)
  bcPred <- terra::predict(bioCropped, bcMod)
  ## plot
  plot(bcPred, main = 'Bioclim')
  plot(world, lwd = .2, add = T)
  plot(region, add = T, border = 'magenta', col = NA)
  # plot(dataSPAT[dataSPAT$presence == 1],
  #      col = 'magenta', pch = '.', add = T)
  ## reclassify based on quantiles
  bcReclass <- quantReclass(bcPred)
  ## plot
  plot(bcReclass, main = 'Bioclim reclassified')
  plot(world, lwd = .2, add = T)
  plot(region, add = T, border = 'magenta', col = NA)
  # plot(dataSPAT[dataSPAT$presence == 1],
  #      col = 'magenta', pch = '.', add = T)
  
  ## save to disk
  saveRDS(bcMod, paste0('output/', folder, '/models/bc_mod.rds'))
  writeRaster(bcPred, filename = paste0('output/', folder, '/models/bc_pred.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  writeRaster(bcReclass, filename = paste0('output/', folder, '/models/bc_reclass.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  
  ## clean unused memory
  gc()
} 

######## MaxEnt - PRESENCE/BACKGROUND (maximum enthropy) --------
if ('MXT' %in% selectedMods$x) {
  sum(is.na(data[, c('presence', vars)]))
  ## model
  mxtMod <- maxnet(p = data[, 'presence'], data = data[, vars],
                   f = maxnet.formula(p = data[, 'presence'], data = data[, vars]))
  ## prediction
  mxtPred <- predict(bioCropped, mxtMod, type = 'cloglog', na.rm = T,
                     cores = 2, cpkgs = 'maxnet')
  ## plot
  plot(mxtPred, main = 'MaxEnt')
  plot(world, lwd = .2, add = T)
  plot(region, add = T, border = 'magenta', col = NA)
  # plot(dataSPAT[dataSPAT$presence == 1],
  #      col = 'magenta', pch = '.', add = T)
  
  ## save to disk
  saveRDS(mxtMod, paste0('output/', folder, '/models/mxt_mod.rds'))
  writeRaster(mxtPred, filename = paste0('output/', folder, '/models/rast/mxt_predStudyArea.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  
  ## clean unused memory
  gc()
}

#### GLM - PRESENCE/ABSENCE (linear fitting) --------
if ('GLM' %in% selectedMods$x) {
  ## formula
  glmForm <- reformulate(termlabels = vars, response = 'presence')
  ## model
  glmMod <- glm(formula = glmForm, family = binomial, data = data)
  ## prediction
  glmPred <- terra::predict(bioCropped, glmMod, type = 'response')
  ## plot
  plot(glmPred, main = 'GLM')
  plot(world, lwd = .2, add = T)
  plot(region, add = T, border = 'magenta', col = NA)
  # plot(dataSPAT[dataSPAT$presence == 1],
  #      col = 'magenta', pch = '.', add = T)
  ## convert to prevalence-independent favourability
  glmFav <- Fav(pred = glmPred, sample.preval = prevalence(model = glmMod))
  ## plot
  plot(glmFav, main = 'GLM favourability', range = c(0, 1))
  plot(world, lwd = .2, add = T)
  plot(region, add = T, border = 'magenta', col = NA)
  # plot(dataSPAT[dataSPAT$presence == 1],
  #      col = 'magenta', pch = '.', add = T)
  
  ## save to disk
  saveRDS(glmMod, paste0('output/', folder, '/models/glm_mod.rds'))
  writeRaster(glmPred, filename = paste0('output/', folder, '/models/rast/glm_predStudyArea.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  writeRaster(glmFav, filename = paste0('output/', folder, '/models/rast/glm_favStudyArea.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  
  ## clean unused memory
  gc()
}

######## GAM - PRESENCE/ABSENCE (non-linear fitting) --------
if ('GAM' %in% selectedMods$x) {
  ## formula
  gamForm <- as.formula(paste0('presence ~ ', paste0('s(', vars, ')', collapse = '+')))
  ## model
  gamMod <- gam(gamForm, family = binomial, data = data)
  ## prediction
  gamPred <- predict(bioCropped, gamMod, type = 'response')
  ## plot
  plot(gamPred, main = 'GAM')
  plot(world, lwd = .2, add = T)
  plot(region, add = T, border = 'magenta', col = NA)
  # plot(dataSPAT[dataSPAT$presence == 1],
  #      col = 'magenta', pch = '.', add = T)
  ## convert to prevalence-independent favourability
  gamFav <- Fav(pred = gamPred, sample.preval = prevalence(model = gamMod))
  ## plot
  plot(gamFav, main = 'GAM favourability', range = c(0, 1))
  plot(world, lwd = .2, add = T)
  plot(region, add = T, border = 'magenta', col = NA)
  # plot(dataSPAT[dataSPAT$presence == 1],
  #      col = 'magenta', pch = '.', add = T)
  
  ## save to disk
  saveRDS(gamMod, paste0('output/', folder, '/models/gam_mod.rds'))
  writeRaster(gamPred, filename = paste0('output/', folder, '/models/rast/gam_predStudyArea.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  writeRaster(gamFav, filename = paste0('output/', folder, '/models/rast/gam_favStudyArea.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  
  ## clean unused memory
  gc()
}

######## BART - PRESENCE/ABSENCE (bayesian sum of trees) --------
if ('BART' %in% selectedMods$x) {
  ## model
  bartMod <- bart(x.train = data[, vars], y.train = data[, 'presence'],
                  ntree = 100, keeptrees = T, verbose = F, seed = 1312)
  ## prediction
  bartPred <- predict2.bart(bartMod, raster::stack(bioCropped), quantiles = c(.05, .95), 
                            splitby = 5, quiet = T)
  summary(bartPred)
  names(bartPred) <- c('BART prediction', 'BART lower.interval', 'BART upper.interval')
  bartPred <- rast(bartPred)
  ## convert to prevalence-independent favourability
  bartFav <- Fav(pred = bartPred, sample.preval = prevalence(model = bartMod))
  names(bartFav) <- c('BART favourability', 'BART favourability lower.interval', 'BART favourability upper.interval')
  ## add the uncertainty range
  bartPred$uncertainty.range <- bartPred[[3]] - bartPred[[2]] # to pred
  bartFav$uncertainty.range <- bartFav[[3]] - bartFav[[2]]    # to fav
  ## plot
  plot(bartPred)
  plot(bartFav)
  plot(bartFav[[1]], main = 'BART mean posterior favourability')
  plot(world, lwd = .2, add = T)
  plot(region, add = T, border = 'magenta', col = NA)
  # plot(dataSPAT[dataSPAT$presence == 1],
  #      col = 'magenta', pch = '.', add = T)
  
  ## BART MODEL CANNOT be saved to rds because it looses information
  writeRaster(bartPred, filename = paste0('output/', folder, '/models/rast/bart_predStudyArea.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  writeRaster(bartFav, filename = paste0('output/', folder, '/models/rast/bart_favStudyArea.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  
  ## clean unused memory
  gc()
}

######## Random Forest - PRESENCE/ABSENCE (average of trees) --------
if ('RF' %in% selectedMods$x) {
  data4rf <- data
  data4rf$presence <- as.factor(data4rf$presence)
  ## formula
  rfForm <- formula(paste0('presence ~ ', paste(vars, collapse = '+')))
  ## model
  rfMod <- randomForest(rfForm, data4rf, na.action = na.exclude)
  ## prediction
  rfPred <- predict(bioCropped, rfMod, type = 'prob', index = 2,
                    cores = 6, cpkgs = 'randomForest')
  ## plot
  plot(rfPred, main = 'RF')
  plot(world, lwd = .2, add = T)
  plot(region, add = T, border = 'magenta', col = NA)
  ## convert to prevalence-independent favourability
  rfFav <- Fav(pred = rfPred, sample.preval = prevalence(model = rfMod))
  ## plot
  plot(rfFav, main = 'RF favourability')
  plot(world, lwd = .2, add = T)
  plot(region, add = T, border = 'magenta', col = NA)
  # plot(dataSPAT[dataSPAT$presence == 1],
  #      col = 'magenta', pch = '.', add = T)
  
  ## save to disk
  saveRDS(rfMod, paste0('output/', folder, '/models/rf_mod.rds'))
  writeRaster(rfPred, filename = paste0('output/', folder, '/models/rast/rf_predStudyArea.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  writeRaster(rfFav, filename = paste0('output/', folder, '/models/rast/rf_favStudyArea.tif'), gdal = c('COMPRESS=DEFLATE'), overwrite = TRUE)
  
  ## clean unused memory
  gc()
}


## time tracking
(end <- Sys.time())
(length <- difftime(end, start, units = 'mins'))
system('say Done')


######## Predictions to df STUDY AREA --------
(start <- Sys.time())

## bioclim
if (exists('bcMod')) {
  data$bcPred <- predict(bcMod, data)
  data$bcReclass <- quantReclass(data$bcPred)
  ## clean unused memory
  gc()
}

## glm
if (exists('glmMod')) {
  data$glmPred <- predict(glmMod, data, type = 'response')
  data$glmFav <- Fav(pred = data$glmPred, sample.preval = prevalence(model = glmMod))
  ## clean unused memory
  gc()
}

## maxent
if (exists('mxtMod')) {
  data$mxtPred <- as.vector(predict(mxtMod, data, type = 'cloglog'))
  ## clean unused memory
  gc()
}

## gam
if (exists('gamMod')) {
  data$gamPred <- as.numeric(unlist(predict(gamMod, data, type = 'response')))
  data$gamFav <- Fav(pred = data$gamPred, sample.preval = prevalence(model = gamMod))
  ## clean unused memory
  gc()
}

## bart
if (exists('bartMod')) {
  source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df")
  bartDFPred <- predict_bart_df(bartMod, data, , quantiles = c(0.05, 0.95))
  names(bartDFPred) <- c('bartPred', 'bartPredLow', 'bartPredUp')
  bartDFFav <- bartDFPred
  for (i in 1:ncol(bartDFFav)) {
    bartDFFav[, i] <- Fav(pred = bartDFPred[, i], sample.preval = prevalence(model = bartMod))
  }
  names(bartDFFav) <- gsub('Pred', 'Fav', names(bartDFFav))
  data <- data.frame(data, bartDFPred, bartDFFav)
  head(data)
  ## clean unused memory
  gc()
}

## random forest
if (exists('rfMod')) {
  data$rfPred <- predict(rfMod, data4rf, type = 'prob', index = 2)[,2]
  data$rfFav <- Fav(pred = data$rfPred, sample.preval = prevalence(model = rfMod))
  ## clean unused memory
  gc()
}

## time tracking
(end <- Sys.time())
(length <- difftime(end, start, units = 'mins'))
system('say Done')


######## Save df to disk --------
data <- data[, -c(grep('PC', names(data)))] # remove bioclim columns
head(data)
write.csv(data, paste0('output/', folder, '/models/allModelsPredictions.csv'), row.names = F)





################# Map all predictions together_________________________________
dataSPAT <- terra::vect(data, geom = c('x', 'y'))
## map all prediction columns in one window:
names(dataSPAT)
predCols <- c('glmFav', 'gamFav', 'bartFav', 'rfFav')
clrs <- hcl.colors(100)
plot(dataSPAT, predCols, cex = 0.1, type = 'continuous', range = c(0, 1),
     axes = FALSE, nc = 2, mar = c(2, 2, 2, 5))
title('Models for H. halys in study area', outer = TRUE, line = -1.5)





################# Clean memory_________________________________
rm(list = setdiff(ls(), c('mxtMod', 'gamMod', 'bartMod', 'rfMod')))
gc()
