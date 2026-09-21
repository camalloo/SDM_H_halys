########################################################################
############ Run the models yearly observations iteratively ############
########################################################################

library(tidyverse)
library(terra)
library(fuzzySim)
library(maxnet)
library(mgcv)
library(dbarts)
library(embarcadero)
library(randomForest)
library(vip)
library(pdp)
source('https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df')

######### Load data ----------
folder <- 'baseline_models_yearly_runs' # name of the mother folder

world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica

preds <- rast(file.path('output', folder,'modellingData/predictors.tif'))
# predsDF <- as.data.frame(preds, na.rm = T, xy = T, cells = T)
bias <- geodata::footprint(year = 2009, path = 'data/geodata/') ## bias layer for absence selection
bias <- project(bias, preds)

native <- read.csv(file.path('output', folder, 'modellingData/native_occs.csv')) |> 
  mutate(year = year(eventDate),
         status = 'n')
invasive <- read.csv(file.path('output', folder, 'modellingData/invasive_occs.csv')) |> 
  mutate(year = year(eventDate),
         status = 'i') |> 
  ## filter obs before 1996 (first invasive discovery)
  filter(year >= 1996)

wd <- readRDS(file.path('output', folder, 'modellingData/ModelRegion/region_mult.rds'))

vars <- read.csv(file.path('output', folder, 'modellingData/selectedVariables.csv'))[, 2] # variables after collinearity analysis
######### 


######### Loop over the years and run the models ----------
## params
saveRootName <- 'yearly_no_weights'             # tag for model run for saved objects (just yearly runs)
absMult <- 10                                   # ratio of absences
weight <- F                                     # logical to use or not weights
dataStorage <- list()                           # save the data used for training in a list
allModels <- NULL                               # models results in df format
allFavAVG <- NULL                            # all means of yearly prediction models
ens <- NULL                                     # initialise ensemble object
regions <- NULL                                 # training regions of all models
allVarimp <- data.frame()                       # variable importance of all models
allRc <- data.frame()                           # response curves

years <- unique(invasive$year)
(years <- sort(years))

## plot the observations across all years
for (y in years) {
  plot(world, main = paste0('obs in year ', y))
  points(vect(native, geom = c('long', 'lat'), crs = crs(world)), col = 'blue', cex = 0.2)
  points(vect(invasive |> filter(year <= y),
              geom = c('long', 'lat'), crs = crs(world)), col = 'red')
}

## overall time tracking
(start <- Sys.time())
for (i in 1:length(years)) {
  y <- years[i]
  cat('Iteration', i, ' -  Year', y, '\n')
  ## yearly time tracking
  startTime <- proc.time()
  
  ## get training and testing dfs ----
  trainDF <- rbind(invasive[invasive$year <= y, ], native)
  #####

  ## modelling data ----
  occs <- vect(trainDF, geom = c('long', 'lat'), crs = crs(preds))
  # plot(occs)

  ## compute the buffer
  b <- buffer(occs, width = wd)
  ## get the region and assign an identifier
  region <- aggregate(b)
  region$year <- y
  # plot(region, col = 'gold', border = NA, main = 'Model training region')
  # plot(world, add = T)
  ## save the region to output stack
  if (is.null(regions)) {
    regions <- region
  } else {
    regions <- rbind(regions, region)
  }
  ## crop predictors
  predsCropped <- crop(preds, region, mask = T, touches = F)
  # plot(predsCropped[[1]])
  # plot(world, add = T)

  ## data gridding
  gridData <- gridRecords(rst = predsCropped, pres.coords = occs)
  ## check if NAs in gridded data
  nas <- gridData |>
    filter(if_any(everything(), is.na))
  # plet(vect(nas, geom = c('x', 'y'), crs = 'EPSG:4326'))
  gridData <- drop_na(gridData)
  ## crop bias layer for absence generation to study area
  biasCropped <- crop(bias, region, mask = T)
  ## absence generation
  data <- selectAbsences(gridData, sp.cols = "presence", coord.cols = c("x", "y"), CRS = crs(occs),
                         mult.p = absMult, seed = 1312, plot = F, df = T, bias = biasCropped, verbosity = 0)
  #####
  
  ## weights on presences and absences (set weight = T in the parameters to use them) ----
  # if (weight == T) {
  #   ## use ensemble of previous years predictions as weight for current year model
  #   if (y == min(unique(invasive$year))) { # don't use weights if iteration is fist year of invasion
  #     data$w <- 1
  #   } else {
  #     ## extract the suitability values
  #     vals <- data.frame(terra::extract(favAVG, vect(data, geom = c("x", "y"), crs = crs(occs))))[,2]
  #     vals <- ifelse(is.na(vals), 0, vals)
  #     ## assign weights only to presences
  #     data$w <- 1 ## all observations get a weight of 1
  #     ## filter for the extracted values of presences only
  #     wPres <- vals[data$presence == 1]
  #     ## assign high weight to low suitable presences (learn from mistakes in previous years)
  #     data$w[data$presence == 1] <- 1 + (1-wPres)
  #     ## filter for the extracted values of absences only
  #     wAbs <- vals[data$presence == 0]
  #     ## assign low weight to high suitable absences (learn from mistakes in previous years)
  #     data$w[data$presence == 0] <- 1 + (1-wAbs)
  #   }
  #   ## rescale weignts to have a mean of 1
  #   data$w <- data$w / mean(data$w)
  # }
  
  ## save the data used for training
  dataStorage[[as.character(y)]] <- data
  #####

  ## model runs ----
  # cat('\tMXT ')
  # ## MAXENT
  # mxtFormula <- maxnet.formula(p = data[, 'presence'], data = data[, vars], classes = 'lqh')
  # data$w_mxt <- ifelse(data$presence == 1, data$w, data$w*100) ## adjust weights for maxent
  # mxtMod <- maxnet(p = data[, 'presence'], data = data[, vars], f = mxtFormula, weights = data$w_mxt,
  #                  regmult = 2)
  # # mxtPred <- predict(preds, mxtMod, type = 'cloglog', na.rm = T, cores = 10, cpkgs = 'maxnet')
  # allModelsDF[[paste0('mxtPred_', y)]] <- as.vector(predict(mxtMod, predsDF[, vars], type = 'cloglog', na.rm = T, cores = 10, cpkgs = 'maxnet'))
  
  cat('\tGAM ')
  ## GAM
  gamForm <- as.formula(paste0('presence ~ ', paste0('s(', vars, ', bs = \'cr\')', collapse = '+')))
  gamMod <- bam(gamForm, family = binomial, data = data, weights = data$w, model = T,
                discrete = T, nthreads = parallel::detectCores() - 1)
  gamPred <- predict(preds, gamMod, type = 'response', na.rm = T)
  names(gamPred) <- paste0('gamPred_', y)
  gamFav <- Fav(pred = gamPred, sample.preval = prevalence(model = gamMod))
  names(gamFav) <- paste0('gamFav_', y)
  # allModelsDF[[paste0('gamPred_', y)]] <- predict(gamMod, predsDF[, vars], type = 'response', na.rm = T)
  # allModelsDF[[paste0('gamFav_', y)]] <- Fav(pred = allModelsDF[[paste0('gamPred_', y)]], 
  #                                            sample.preval = prevalence(model = gamMod))

  cat('- BART ')
  ## BART
  bartMod <- bart(x.train = data[, vars], y.train = data[, 'presence'], weights = data$w,
                  ntree = 200, keeptrees = T, verbose = F, seed = 1312)
  bartPred <- predict2.bart(bartMod, raster::stack(preds), quantiles = c(), splitby = 10, quiet = T)
  bartPred <- rast(bartPred)
  crs(bartPred) <- crs(preds)
  names(bartPred) <- paste0('bartPred_', y)
  bartFav <- Fav(pred = bartPred, sample.preval = prevalence(model = bartMod))
  names(bartFav) <- paste0('bartFav_', y)
  # allModelsDF[[paste0('bartPred_', y)]] <- predict_bart_df(bartMod, predsDF[, vars])
  # allModelsDF[[paste0('bartFav_', y)]] <- Fav(pred = allModelsDF[[paste0('bartPred_', y)]], 
  #                                             sample.preval = prevalence(model = bartMod))

  cat('- RF\n')
  ## RF
  rfForm <- formula(paste0('presence ~ ', paste(vars, collapse = '+')))
  rfMod <- randomForest(rfForm, data, weights = data$w, na.action = na.exclude)
  rfPred <- predict(preds, rfMod, na.rm = T, cores = parallel::detectCores() - 1, cpkgs = 'randomForest')
  names(rfPred) <- paste0('rfPred_', y)
  rfPred <- clamp(rfPred, lower = 0, upper = 1) # clamp to 0,1
  rfFav <- Fav(pred = rfPred, sample.preval = prevalence(model = rfMod))
  names(rfFav) <- paste0('rfFav_', y)
  # allModelsDF[[paste0('rfPred_', y)]] <- predict(rfMod, predsDF[, vars], na.rm = T, cores = 10, cpkgs = 'randomForest')
  # allModelsDF[[paste0('rfPred_', y)]] <- pmin(pmax(allModelsDF[[paste0('rfPred_', y)]], 0), 1)
  # allModelsDF[[paste0('rfFav_', y)]] <- Fav(pred = allModelsDF[[paste0('rfPred_', y)]], 
  #                                           sample.preval = prevalence(model = rfMod))
  
  ## save the models together
  allYear <- c(gamPred, gamFav,
               bartPred, bartFav,
               rfPred, rfFav)
  if (is.null(allModels)) {
    allModels <- allYear
  } else {
    allModels <- c(allModels, allYear)
  }
  #####
  
  
  ## subsample of data & wrapper functions ----
  ## subsample the data for variable importance and response curves
  set.seed(1312)
  subdata <- data |>
    dplyr::slice_sample(n = 5000)
  
  # mxt_pred_wrapper <- function(object, newdata) {
  #   p <- predict(object, newdata, type = 'cloglog')
  #   drop(p)
  # }
  gam_pred_wrapper <- function(object, newdata) {
    p <- predict(object, newdata, type = 'response')
  }
  rf_pred_wrapper <- function(object, newdata) {
    p <- predict(object, newdata)
  }
  bart_pred_wrapper <- function(object, newdata) {
    p <- predict_bart_df(object, newdata)
  }
  #####
  
  ## variable importance ----
  ## maxent
  # MXTvarimp <- vi_permute(mxtMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
  #                         pred_wrapper = mxt_pred_wrapper, nsim = 1,
  #                         parallel = T, .packages = 'maxnet')
  # vip(MXTvarimp) + ggtitle('MXT')
  # MXTvarimp <- cbind(mod = 'MXT', year = y, MXTvarimp)
  ## gam
  GAMvarimp <- vi_permute(gamMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
                          pred_wrapper = gam_pred_wrapper, nsim = 10,
                          parallel = T, .packages = 'mgcv')
  # vip(GAMvarimp) + ggtitle('GAM')
  GAMvarimp <- cbind(mod = 'GAM', year = y, GAMvarimp)
  ## rf
  RFvarimp <- vi_permute(rfMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
                         pred_wrapper = rf_pred_wrapper, nsim = 10,
                         parallel = T, .packages = 'randomForest')
  # vip(RFvarimp) + ggtitle('RF')
  RFvarimp <- cbind(mod = 'RF', year = y, RFvarimp)
  ## bart
  BARTvarimp <- vi_permute(bartMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
                           pred_wrapper = bart_pred_wrapper, nsim = 10,
                           parallel = T, .packages = 'dbarts')
  # vip(BARTvarimp) + ggtitle('BART')
  BARTvarimp <- cbind(mod = 'BART', year = y, BARTvarimp)
  ## append to results
  # allVarimp <- rbind(allVarimp, MXTvarimp, GAMvarimp, RFvarimp, BARTvarimp)
  allVarimp <- rbind(allVarimp, GAMvarimp, RFvarimp, BARTvarimp)
  cat('\tvariable importance done\n')
  #####
  
  ## response curves ----
  ## maxent
  # MXTrc <- lapply(vars, function(v) {
  #   p <- pdp::partial(mxtMod, pred.var = v, train = subdata[, vars],
  #                     grid.resolution = 50, pred.fun = mxt_pred_wrapper)
  #   p <- p |> 
  #     group_by(pick(1)) |> 
  #     summarise(yhat = mean(yhat),
  #               .groups = 'drop')
  #   ## out
  #   data.frame(mod = 'MXT',
  #              year = y,
  #              var = v,
  #              x = p[[v]],
  #              y = p$yhat)
  # }) |> bind_rows()
  ## gam
  GAMrc <- lapply(vars, function(v) {
    p <- pdp::partial(gamMod, pred.var = v, train = subdata[, vars],
                      grid.resolution = 50, pred.fun = gam_pred_wrapper)
    p <- p |> 
      group_by(pick(1)) |> 
      summarise(yhat = mean(yhat),
                .groups = 'drop')
    ## out
    data.frame(mod = 'GAM',
               year = y,
               var = v,
               x = p[[v]],
               y = p$yhat)
  }) |> bind_rows()
  ## rf
  RFrc <- lapply(vars, function(v) {
    p <- pdp::partial(rfMod, pred.var = v, train = subdata[, vars],
                      grid.resolution = 50, prob = T, pred.fun = rf_pred_wrapper)
    p <- p |> 
      group_by(pick(1)) |> 
      summarise(yhat = mean(yhat),
                .groups = 'drop')
    ## out
    data.frame(mod = 'RF',
               year = y,
               var = v,
               x = p[[v]],
               y = p$yhat)
  }) |> bind_rows()
  ## bart
  BARTrc <- lapply(vars, function(v) {
    p <- pdp::partial(bartMod, pred.var = v, train = subdata[, vars],
                      grid.resolution = 50, prob = T, pred.fun = bart_pred_wrapper)
    p <- p |> 
      group_by(pick(1)) |> 
      summarise(yhat = mean(yhat),
                .groups = 'drop')
    ## out
    data.frame(mod = 'BART',
               year = y,
               var = v,
               x = p[[v]],
               y = p$yhat)
  }) |> bind_rows()
  ## append to results df
  # allRc <- rbind(allRc, MXTrc, GAMrc, RFrc, BARTrc)
  allRc <- rbind(allRc, GAMrc, RFrc, BARTrc)
  cat('\tresponse curve done\n')
  #####
  
  ## ensemble ----
  favMean <- app(c(bartFav[[1]], gamFav, rfFav), mean)
  names(favMean) <- paste0('mods_mean_fav_', y)
  favVar <- app(c(bartFav[[1]], gamFav, rfFav), var)
  names(favVar) <- paste0('mods_var_fav_', y)
  favAVG <- c(favMean, favVar)
  if (is.null(allFavAVG)) {
    allFavAVG <- favAVG
  } else {
    allFavAVG <- c(allFavAVG, favAVG)
  }
  #####
  
  ## yearly time tracking
  endTime <- proc.time()
  cat('Elapsed year', y, 'time:', difftime(endTime[3], startTime[3], units = 'mins'), 'mins')
  cat('\n\n')
}
## time tracking
(end <- Sys.time())
(length <- difftime(end, start, units = 'auto'))
## check if all models are there (6 models x 23 years = 138 models)
names(allModels)
names(allFavAVG)
######### 


######### Save the results to file ----------
## save the training data
saveRDS(dataStorage, file = file.path('output', folder, 'models', paste0(saveRootName, '_dataStorage.Rdata')))

## save the rasters
writeRaster(allModels, filename = file.path('output', folder, 'models/rast', paste0(saveRootName, '_predictions.tif')), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)

## save the dataframe
# write.csv(allModelsDF, file.path('output', folder, 'models/rast/yearly_predictionsDF.csv'), row.names = F)

## save models averages
writeRaster(allFavAVG, filename = file.path('output', folder, 'models/rast', paste0(saveRootName, '_mods_avg_fav.tif')), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
## save training regions
writeVector(regions, file.path('output', folder, 'modellingData/ModelRegion/yearly_modellingRegions.gpkg'), overwrite = T)

## save the variable importance
write.csv(allVarimp, file.path('output', folder, 'models/evaluation', paste0(saveRootName, '_varImp.csv')), row.names = F)

## save response curves
write.csv(allRc, file.path('output', folder, 'models/evaluation', paste0(saveRootName, '_Rc.csv')), row.names = F)

######### 


######### Plot a model species for all years ----------
mod <-'gamPred_'
plot(allFavAVG[['ens_0']], main = 'ens native')
plot(allFavAVG[['ens_1']], main = 'ens invasive')
plot(allFavAVG[['ens_2025']], main = 'ens 2025')
for (y in years) {
  plot(allFavAVG[[paste0('ens_', y)]], main = paste0('ens ', y))
  points(vect(native, geom = c('long', 'lat'), crs = crs(ens)), col = 'red', cex = 0.1)
  points(vect(invasive |> filter(year <= y), geom = c('long', 'lat'), crs = crs(ens)), col = 'red', cex = 0.1)
}

plet(vect(invasive, geom = c('long', 'lat'), crs = crs(world)), col = 'red')
