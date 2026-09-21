#################################################################################################
############ Run the models with native observations and project on the entire world ############
#################################################################################################

library(terra)
library(fuzzySim)
library(maxnet)
library(mgcv)
library(dbarts)
library(embarcadero)
library(randomForest)
library(vip)
library(pdp)
library(future)

folder <- 'baseline_models_yearly_runs' # name of the mother folder

######### Load data ----------
world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica

native <- read.csv(file.path('output', folder, 'modellingData/native_occs.csv'))
occs <- vect(native, geom = c('long', 'lat'), crs = crs(world))
preds <- rast(file.path('output', folder,'modellingData/predictors.tif'))

selectedMods <- read.csv(file.path('output', folder, 'models/crossValidation/selectedModels.csv'))
vars <- read.csv(file.path('output', folder, 'modellingData/selectedVariables.csv'))[, 2] # variables after collinearity analysis
## modelling region multiplier
wd <- readRDS(file.path('output', folder, 'modellingData/ModelRegion/region_mult.rds'))
######### 


######### Compute the modelling data ----------  
b <- buffer(occs, width = wd)
region <- aggregate(b)
region <- crop(region, world)
plot(world, lwd = 0, main = 'Model training region - native', font.main = 1)
plot(region, col = 'gold', border = NA, add = T)
plot(world, add = T)
points(occs, cex = 0.1)
## save to disk
writeVector(region, file.path('output', folder, 'modellingData/ModelRegion/native_modellingRegionPoly.shp'), overwrite = T)

## crop predictors according to model region
predsCropped <- crop(preds, region, mask = T, touches = F)
plot(predsCropped)
## save to disk
writeRaster(predsCropped, filename = file.path('output', folder, 'modellingData/native_predsCropped.tif'),
            filetype = 'GTiff', datatype = 'FLT4S', gdal = c('COMPRESS=DEFLATE'), overwrite = T)
gc()

## data gridding
gridData <- gridRecords(rst = predsCropped, pres.coords = occs)
nrow(gridData)
## clean from NAs in predictors if any
nas <- gridData |> 
  filter(if_any(everything(), is.na))
if (nrow(nas) > 0) {
  nasSPAT <- vect(nas, geom = c('x', 'y'), crs = crs(world))
  gridData <- drop_na(gridData)
  
  plet(predsCropped[['dem']]) |>
    points(occs, col = 'red') |>
    points(nasSPAT)
}
table(gridData$presence)
# save the raw gridded data to csv so I have them ready to go
write.csv(data.frame(gridData), file.path('output', folder, 'modellingData/native_GriddedData.csv'), row.names = F)

## modelling data
nrow(gridData)
table(gridData$presence)  # number of presences and absences
bias <- geodata::footprint(year = 2009, path = 'data/geodata/')
bias <- project(bias, preds)
bias <- crop(bias, region, mask = T)
plot(bias)
data <- selectAbsences(gridData,
                       sp.cols = "presence", coord.cols = c("x", "y"), CRS = crs(occs),
                       mult.p = 10, seed = 1312, df = T, bias = bias)
plot(world, add = T, main = 'Modelling data - native')
nrow(data)
table(data$presence)
## save test and train dfs to disk
write.csv(data.frame(data), file.path('output', folder, 'modellingData/native_modData.csv'), row.names = F)
######### 


######### Run complete models on STUDY AREA ---------- 
## selected models
selectedMods$mod
## time tracking
(start <- Sys.time())

######## MAXENT - PRESENCE/BACKGROUND (maximum entropy)
# ## formula
# mxtFormula <- maxnet.formula(p = data[, 'presence'], data = data[, vars])
# ## model
# mxtMod <- maxnet(p = data[, 'presence'], data = data[, vars],
#                  f = mxtFormula)
# ## prediction
# mxtPred <- predict(predsCropped, mxtMod, type = 'cloglog', na.rm = T,
#                    cores = 10, cpkgs = 'maxnet')
# names(mxtPred) <- 'MXT.pred'
# ## save model to disk
# saveRDS(mxtMod, file.path('output', folder, 'models/native_mxt_mod.rds'))
# ## clean unused memory
# gc()

######## GAM - PRESENCE/ABSENCE (non-linear fitting)
## formula
gamForm <- as.formula(paste0('presence ~ ', paste0('s(', vars, ', bs = "cs", k = 20)', collapse = '+')))
## model
gamMod <- gam(gamForm, family = binomial, data = data, model = T, method = 'REML')
gam.check(gamMod)
## prediction
gamPred <- predict(predsCropped, gamMod, type = 'response', na.rm = T)
names(gamPred) <- 'GAM.pred'
## convert to prevalence-independent favourability
gamFav <- Fav(pred = gamPred, sample.preval = prevalence(model = gamMod))
names(gamFav) <- 'GAM.fav'
## save model to disk
saveRDS(gamMod, file.path('output', folder, 'models/native_gam_mod.rds'))
## clean unused memory
gc()

######## BART - PRESENCE/ABSENCE (bayesian sum of trees)
## model
bartMod <- bart(x.train = data[, vars], y.train = data[, 'presence'],
                ntree = 200, keeptrees = T, verbose = F, seed = 1312)
## prediction
bartPred <- predict2.bart(bartMod, raster::stack(predsCropped), quantiles = c(.05, .95), 
                          quiet = T)
summary(bartPred)
names(bartPred) <- c('BART.pred', 'BART.lower.interval', 'BART.upper.interval')
bartPred <- rast(bartPred)
## convert to prevalence-independent favourability
bartFav <- Fav(pred = bartPred, sample.preval = prevalence(model = bartMod))
names(bartFav) <- c('BART.fav', 'BART.fav.lower.interval', 'BART.fav.upper.interval')
## add the uncertainty range
bartPred$BART.pred.uncertainty.range <- bartPred[[3]] - bartPred[[2]] # to pred
bartFav$BART.fav.uncertainty.range <- bartFav[[3]] - bartFav[[2]]    # to fav
## BART MODEL CANNOT be saved to rds because it looses information
## clean unused memory
gc()

######## Random Forest - PRESENCE/ABSENCE (average of trees)
## formula
rfForm <- formula(paste0('presence ~ ', paste(vars, collapse = '+')))
## model
rfMod <- randomForest(rfForm, data, na.action = na.exclude)
## prediction
rfPred <- predict(predsCropped, rfMod, na.rm = T, cores = 10, cpkgs = 'randomForest')
rfPred <- clamp(rfPred, lower = 0, upper = 1) # clamp to 0,1
names(rfPred) <- 'RF.pred'
## convert to prevalence-independent favourability
rfFav <- Fav(pred = rfPred, sample.preval = prevalence(model = rfMod))
names(rfFav) <- 'RF.fav'
## save to disk
saveRDS(rfMod, file.path('output', folder, 'models/native_rf_mod.rds'))
## clean unused memory
gc()

## stack predictions in a single raster
studyArea <- c(gamPred, gamFav,
               bartPred[['BART.pred']], bartFav[['BART.fav']],
               rfPred, rfFav)
plot(studyArea)
## save to disk
writeRaster(studyArea, filename = file.path('output', folder, 'models/rast/native_predictionsStudyArea.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)

## predict to dataframe
predictionsDF <- data[, c(1:4)]
# predictionsDF$mxt.pred <- as.vector(predict(mxtMod, data, type = 'cloglog', na.rm = T,
#                                             cores = 10, cpkgs = 'maxnet'))
predictionsDF$gam.pred <- predict(gamMod, data, type = 'response', na.rm = T)
predictionsDF$gam.fav <- Fav(pred = predictionsDF$gam.pred, 
                             sample.preval = prevalence(model = gamMod))
predictionsDF$rf.pred <- predict(rfMod, data, na.rm = T,
                                 cores = 10, cpkgs = 'randomForest')
predictionsDF$rf.pred <- pmin(pmax(predictionsDF$rf.pred, 0), 1)
predictionsDF$rf.fav <- Fav(pred = predictionsDF$rf.pred, 
                            sample.preval = prevalence(model = rfMod))
source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df")
predictionsDF$bart.pred <- predict_bart_df(bartMod, data)
predictionsDF$bart.fav <- Fav(pred = predictionsDF$bart.pred, 
                              sample.preval = prevalence(model = bartMod))
## save to file
write.csv(predictionsDF, file.path('output', folder, 'models/native_predictionsStudyArea_DF.csv'), row.names = F)

## time tracking
(end <- Sys.time())
(length <- difftime(end, start, units = 'mins'))
#########


######### Evaluate models on study area ----------
# ## load models 
# mxtMod <- readRDS(file.path('output', folder, 'models/native_mxt_mod.rds'))
# gamMod <- readRDS(file.path('output', folder, 'models/native_gam_mod.rds'))
# rfMod <- readRDS(file.path('output', folder, 'models/native_rf_mod.rds'))
# ## reload bart because it cannot be saved
# bartMod <- bart(x.train = data[, vars], y.train = data[, 'presence'],
#                 ntree = 200, keeptrees = T, verbose = F, seed = 1312)

## subsample the data for variable importance and response curves
set.seed(1312)
subdata <- data |>
  dplyr::slice_sample(n = 5000)
table(subdata$presence)

## load predictions
studyArea <- rast(file.path('output', folder, 'models/rast/native_predictionsStudyArea.tif'))
plot(studyArea[[c('GAM.fav', 'BART.fav', 'RF.fav')]])

##### correlation between predictions ----
correlation <- cor(predictionsDF[, c('gam.fav', 'rf.fav', 'bart.fav')])
corrplot::corrplot(correlation, method = 'ellipse', type = 'upper',
                   addCoef.col = 'wheat3', addCoefasPercent = F,
                   tl.col = 'black', diag = T)
#####

##### variable importance ----
# ## maxent
# mxt_pred_wrapper <- function(object, newdata) {
#   p <- predict(object, newdata, type = 'cloglog')
#   drop(p)
# }
# MXTvarimp <- vi_permute(mxtMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
#                         pred_wrapper = mxt_pred_wrapper, nsim = 1,
#                         parallel = T,.packages = 'maxnet')
# vip(MXTvarimp) + ggtitle('MXT')
# MXTvarimp <- cbind(mod = 'MXT', MXTvarimp)
# cat('maxent variable importance done\n')

## gam
gam_pred_wrapper <- function(object, newdata) {
  p <- predict(object, newdata, type = 'response')
}
GAMvarimp <- vi_permute(gamMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
                        pred_wrapper = gam_pred_wrapper, nsim = 1,
                        parallel = T, .packages = 'mgcv')
vip(GAMvarimp) + ggtitle('GAM')
GAMvarimp <- cbind(mod = 'GAM', GAMvarimp)
cat('gam variable importance done\n')

## random forest
rf_pred_wrapper <- function(object, newdata) {
  p <- predict(object, newdata)
}
RFvarimp <- vi_permute(rfMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
                       pred_wrapper = rf_pred_wrapper, nsim = 1,
                       parallel = T, .packages = 'randomForest')
vip(RFvarimp) + ggtitle('RF')
RFvarimp <- cbind(mod = 'RF', RFvarimp)
cat('random forest variable importance done\n')

## bart
source('https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df')
bart_pred_wrapper <- function(object, newdata) {
  p <- predict_bart_df(object, newdata)
}
BARTvarimp <- vi_permute(bartMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
                         pred_wrapper = bart_pred_wrapper, nsim = 1,
                         parallel = T, .packages = 'dbarts')
vip(BARTvarimp) + ggtitle('BART')
BARTvarimp <- cbind(mod = 'BART', BARTvarimp)
cat('bart variable importance done\n')

## save all to df
ALLvarimp <- rbind(GAMvarimp, RFvarimp, BARTvarimp)
ggplot(ALLvarimp, aes(x = Importance, y = reorder(Variable, Importance))) +
  geom_col() +
  facet_wrap(~mod)
write.csv(ALLvarimp, file.path('output', folder, 'models/evaluation/native_varimp.csv'))
#####

##### response curves ----
# ## maxent
# MXTrc <- lapply(vars, function(v) {
#   p <- pdp::partial(mxtMod, pred.var = v, train = subdata[, vars],
#                     grid.resolution = 50, pred.fun = mxt_pred_wrapper)
#   p <- p |> 
#     group_by(pick(1)) |> 
#     summarise(yhat = mean(yhat),
#               .groups = 'drop')
#   ## out
#   data.frame(mod = 'MXT',
#              var = v,
#              x = p[[v]],
#              y = p$yhat)
# }) |> bind_rows()
# ggplot(MXTrc, aes(x = x, y = y)) +
#   geom_line() +
#   facet_wrap(~var, scales = 'free') +
#   ggtitle('MXT')

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
             var = v,
             x = p[[v]],
             y = p$yhat)
}) |> bind_rows()
ggplot(GAMrc, aes(x = x, y = y)) +
  geom_line() +
  facet_wrap(~var, scales = 'free') +
  ggtitle('GAM')

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
             var = v,
             x = p[[v]],
             y = p$yhat)
}) |> bind_rows()
ggplot(RFrc, aes(x = x, y = y)) +
  geom_line() +
  facet_wrap(~var, scales = 'free') +
  ggtitle('RF')

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
             var = v,
             x = p[[v]],
             y = p$yhat)
}) |> bind_rows()
ggplot(BARTrc, aes(x = x, y = y)) +
  geom_line() +
  facet_wrap(~var, scales = 'free') +
  ggtitle('BART')

## save all response curves to files
ALLrc <- rbind(GAMrc, RFrc, BARTrc)
write.csv(ALLrc, file.path('output', folder, 'models/evaluation/native_response_curves.csv'))

## plot all response curves together
ggplot(ALLrc, aes(x = x, y = y)) +
  geom_line(aes(colour = mod)) +
  facet_wrap(~var, scale = 'free') +
  ggtitle('Partial dependence functions')
#####

##### evaluate complete models on AUC & Miller ----
predictionsDF <- read.csv(file.path('output', folder, 'models/native_predictionsStudyArea_DF.csv'))
evalALL <- data.frame()
for (mod in c('gam.fav', 'rf.fav', 'bart.fav')) {
  cat('Model:', mod, '\n')
  
  ## auc
  auc <- pROC::roc(response = predictionsDF[, 'presence'], predictor = predictionsDF[, mod])
  pROC::plot.roc(auc)
  auc$auc
  # ## tss
  # tss <- modEvA::threshMeasures(obs = predictionsDF[, 'presence'], pred = predictionsDF[, mod],
  #                               thresh = 'maxTSS', measures = 'TSS',
  #                               simplif = TRUE, plot = F, standardize = FALSE)
  # tss
  ## miller
  miller <- modEvA::MillerCalib(obs = predictionsDF[, 'presence'], pred = predictionsDF[, mod],
                                main = paste(mod, 'Miller'))$slope
  miller
  
  ## save to df
  evalALL <- rbind(evalALL, data.frame(model = mod,
                                       auc = auc$auc,
                                       # tss = tss,
                                       miller = miller))
}
write.csv(evalALL, file.path('output', folder, 'models/evaluation/native_evaluationStudyArea.csv'), row.names = F)
#####
######### 


######### Run complete models on WORLD ----------
## selected models
selectedMods$mod
## time tracking
(start <- Sys.time())

######## MAXENT - PRESENCE/BACKGROUND (maximum entropy)
# ## prediction
# mxtPred <- predict(preds, mxtMod, type = 'cloglog', na.rm = T, cores = 10, cpkgs = 'maxnet',
#                    filename = file.path('output', folder, 'models/rast/native_mxt_predWorld.tif'), overwrite = T,
#                    wopt = list(gdal = c('COMPRESS=DEFLATE')))
# ## plot
# plot(mxtPred, main = 'MXT')
# plot(world, lwd = .2, add = T)
# plot(region, add = T, border = 'magenta', col = NA)
# ## clean unused memory
# gc()

######## GAM - PRESENCE/ABSENCE (non-linear fitting)
## prediction + save to raster
gamPred <- predict(preds, gamMod, type = 'response', na.rm = T,
                   filename = file.path('output', folder, 'models/rast/native_gam_predWorld.tif'), overwrite = T, 
                   wopt = list(gdal = c('COMPRESS=DEFLATE')))
## plot
plot(gamPred, main = 'GAM')
plot(world, lwd = .2, add = T)
plot(region, add = T, border = 'magenta', col = NA)
## convert to prevalence-independent favourability
gamFav <- Fav(pred = gamPred, sample.preval = prevalence(model = gamMod))
## plot
plot(gamFav, main = 'GAM favourability', range = c(0, 1))
plot(world, lwd = .2, add = T)
plot(region, add = T, border = 'magenta', col = NA)
## save to disk
writeRaster(gamFav, filename = file.path('output', folder, 'models/rast/native_gam_favWorld.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
## clean unused memory
gc()

######## BART - PRESENCE/ABSENCE (bayesian sum of trees)
## prediction
bartPred <- predict2.bart(bartMod, raster::stack(preds), quantiles = c(.05, .95), 
                          splitby = 10, quiet = T)
summary(bartPred)
names(bartPred) <- c('BART.pred', 'BART.lower.interval', 'BART.upper.interval')
bartPred <- rast(bartPred)
## convert to prevalence-independent favourability
bartFav <- Fav(pred = bartPred, sample.preval = prevalence(model = bartMod))
names(bartFav) <- c('BART.fav', 'BART.fav.lower.interval', 'BART.fav.upper.interval')
## add the uncertainty range
bartPred$BART.pred.uncertainty.range <- bartPred[[3]] - bartPred[[2]] # to pred
bartFav$BART.fav.uncertainty.range <- bartFav[[3]] - bartFav[[2]]    # to fav
## plot
plot(bartPred)
plot(bartFav)
plot(bartFav[[1]], main = 'BART mean posterior favourability')
plot(world, lwd = .2, add = T)
plot(region, add = T, border = 'magenta', col = NA)
## save to disk
writeRaster(bartPred, filename = file.path('output', folder, 'models/rast/native_bart_predWorld.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
writeRaster(bartFav, filename = file.path('output', folder, 'models/rast/native_bart_favWorld.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
## clean unused memory
gc()

######## Random Forest - PRESENCE/ABSENCE (average of trees)
## prediction + save to raster
rfPred <- predict(preds, rfMod, na.rm = T, cores = 10, cpkgs = 'randomForest',
                  filename = file.path('output', folder, 'models/rast/native_rf_predWorld.tif'), overwrite = T,
                  wopt = list(gdal = c('COMPRESS=DEFLATE')))
rfPred <- clamp(rfPred, lower = 0, upper = 1) # clamp to 0,1
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
## save to disk
writeRaster(rfFav, filename = file.path('output', folder, 'models/rast/native_rf_favWorld.tif'),
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
## clean unused memory
gc()

## time tracking
(end <- Sys.time())
(length <- difftime(end, start, units = 'mins'))
######### 


######### Ensemble ----------
## create a combined score 
# evalALL <- read.csv(file.path('output', folder, 'models/evaluation/evaluationStudyArea.csv'))
# scores <- evalALL |> 
#   mutate(miller = exp(-abs(miller - 1)),
#          combined = sqrt(auc * miller))
# whts <- setNames(scores$combined, scores$model)
# allproj <- c(bartFav[[1]], mxtPred, gamFav, rfFav)
# favEns <- weighted.mean(allproj, w = whts)
bartFav <- rast(file.path('output', folder, 'models/rast/native_bart_favWorld.tif'))
gamFav <- rast(file.path('output', folder, 'models/rast/native_gam_favWorld.tif'))
rfFav <- rast(file.path('output', folder, 'models/rast/native_rf_favWorld.tif'))

favMean <- app(c(bartFav[[1]], gamFav, rfFav), mean)
names(favMean) <- 'mods_mean_fav'
favVar <- app(c(bartFav[[1]], gamFav, rfFav), var)
names(favVar) <- 'mods_var_fav'
favAVG <- c(favMean, favVar)
writeRaster(favAVG, filename = file.path('output', folder, 'models/rast/native_mods_avg_fav.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
plot(favAVG)
######### 


######## Clean memory --------
rm(list = ls())
gc()
