#########################################################################
############ Perform cross validation and select best models ############
#########################################################################

library(terra)
library(blockCV)
library(fuzzySim)
library(modEvA)
library(predicts)
library(maxnet)
library(mgcv)
library(dbarts)
library(randomForest)


######### Load data ----------
folder <- 'climate_human' # name of the mother folder

predsCropped <- rast(file.path('output', folder, 'modellingData/predsCropped.tif'))
data <- read.csv(file.path('output', folder, 'modellingData/modData.csv')) # train dataset
dataSPAT <- vect(data, geom = c('x', 'y'), crs = 'EPSG:4326') # spatialize it

world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica

vars <- read.csv(file.path('output', folder, 'modellingData/selectedVariables.csv'))[, 2]

## check data
# plot(predsCropped[[1]])
# 
# nrow(data)
# table(data$presence)
# plot(world, main = 'Presence and background points')
# points(subset(data, data$presence == 0, select = c('x', 'y')), col = 'pink3', cex = 0.08)
# points(subset(data, data$presence == 1, select = c('x', 'y')), col = 'blue', cex = 0.08)
# plot(world, lw = 0.3, add = T)
# legend('bottomleft', horiz = T, inset = c(0, 0.1), xpd = T,
#        bty = 'n', bg = 'white', cex = 0.8,
#        legend = c('presence', 'background'), pch = 21, col = 'black', pt.bg = c('blue', 'pink3'))
# 
# vars
######### 


######## Blocks creation --------
## create blocks and define the number of folds for cross-validation
dataSPATsf <- sf::st_as_sf(data[1:3], coords = c('x', 'y'), crs = 'EPSG:4326') # transform to sf object
dataSPATsf <- sf::st_transform(dataSPATsf, crs = "+proj=moll") # use Mollweide projection for euclidean distances
## use spatial autocorrelation of points to select blocks size quantitatively
spatAutocorrPoints <- cv_spatial_autocor(x = dataSPATsf, column = 'presence', 
                                         plot = T, progress = T)
blockSize <- spatAutocorrPoints$range
rm(dataSPATsf)     # remove sf object because not needed
## blocks area in km
blockSize/1000

## clean unused memory
gc()

## define the k folds for cv and divide the region in blocks of k classes
k = 5
## create the blocks and divide them in five folds
blocks <- cv_spatial(dataSPAT,
                     column = 'presence', r = predsCropped, size = blockSize, 
                     k = k, selection = 'random', hexagon = T, seed = 1312, iteration = 1000)
## assign the belonging fold to each observation in the dataset (also the spatial one for plotting)
data$fold <- blocks$folds_ids
dataSPAT$fold <- blocks$folds_ids

## check the number of test data for each fold
(n <- table(data[, c('fold', 'presence')]))
(nt <- rowSums(n))  # n total rows per fold
(np <- n[ , 2])  # n presence rows per fold

## plot the folds and presence points
plot(dataSPAT, 'fold', cex = 0.4, col = hcl.colors(k, 'viridis'), legend = F)
plot(world, add = T, lw = .4)
# points(subset(dataSPAT, dataSPAT$presence == 1), pch = 4, cex = 0.2)
plot(vect(blocks$blocks), border = "blue", add = TRUE, lw = .7)
# legend("bottom", inset = c(0, -0.05), xpd = T, bty = 'n',
#        legend = paste("block size:", round(blockSize / 1000), "km"))
# legend("top", inset = c(0, 0), xpd = T, cex = 0.7, bg = 'white', horiz = T, 
#        legend = paste("fold", 1:k, ": N =", np, "/", nt), 
#        fill = hcl.colors(k, 'Peach'))

## save the blocked data to csv
write.csv(data.frame(data), file.path('output', folder, 'modellingData/FoldedData.csv'), row.names = F)

## clean unused memory
gc()
######### 


######## Block Cross Validation training --------
data <- read.csv(file.path('output', folder, 'modellingData/FoldedData.csv'))
folds <- sort(unique(data$fold))
## create a df to store results
predsCV <- data[ , c(1:4, grep("fold", names(data)))]
## create a df with as.factor(data$presence) for the RF model
data4rf <- data
data4rf$presence <- as.factor(data4rf$presence)

## iterate through the folds and train the models
source('https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df')
for (f in folds) {
  startTime <- proc.time() # start time of the fold
  ## train the model on all data but a fold f
  ## but do the prediction using all data
  foldTrainingData <- subset(data, data$fold != f)
  
  message(paste('fold', f, ' - Bioclim'))
  bcMod_fold <- envelope(x = subset(foldTrainingData, presence == 1, select = vars))
  bcPred_fold <- predict(bcMod_fold, data)
  predsCV[ , paste0('BC_fold_', f, '_pred')] <- as.vector(bcPred_fold)
  predsCV[ , paste0('BC_fold_', f, '_pred_reclass')] <- quantReclass(bcPred_fold)
  
  message(paste('fold', f, ' - Maxent'))
  mxtFormula <- maxnet.formula(foldTrainingData[, 'presence'], foldTrainingData[, vars])
  mxtMod_fold <- maxnet(p = foldTrainingData[, 'presence'], data = foldTrainingData[, vars],
                        f = mxtFormula)
  mxtPred_fold <- predict(mxtMod_fold, data, type = 'cloglog')
  predsCV[ , paste0('MXT_fold_', f, '_pred')] <- as.vector(mxtPred_fold)
  
  message(paste('fold', f, ' - GLM'))
  glmFormula <- reformulate(termlabels = vars, response = 'presence')
  glmMod_fold <- glm(formula = glmFormula, family = binomial, data = foldTrainingData)
  glmPred_fold <- predict(glmMod_fold, data, type = 'response')
  glmFav_fold <- Fav(pred = glmPred_fold, sample.preval = prevalence(model = glmMod_fold))
  predsCV[ , paste0('GLM_fold_', f, '_pred')] <- as.vector(glmPred_fold)
  predsCV[ , paste0('GLM_fold_', f, '_fav')] <- as.vector(glmFav_fold)
  
  message(paste('fold', f, ' - GAM'))
  gamFormula <- as.formula(paste0('presence ~ ', paste0('s(', vars, ', bs = \'cr\')', collapse = '+')))
  gamMod_fold <- gam(gamFormula, family = binomial, data = foldTrainingData)
  gamPred_fold <- predict(gamMod_fold, data, type = 'response')
  gamFav_fold <- Fav(pred = gamPred_fold, sample.preval = prevalence(model = gamMod_fold))
  predsCV[ , paste0('GAM_fold_', f, '_pred')] <- as.vector(gamPred_fold)
  predsCV[ , paste0('GAM_fold_', f, '_fav')] <- as.vector(gamFav_fold)
  
  message(paste('fold', f, ' - BART'))
  bartMod_fold <- bart(x.train = foldTrainingData[, vars], y.train = foldTrainingData[, 'presence'],
                       keeptrees = TRUE, verbose = FALSE, seed = 1312)
  bartPred_fold <- predict_bart_df(bartMod_fold, data)
  bartFav_fold <- Fav(pred = bartPred_fold, sample.preval = prevalence(model = bartMod_fold))
  predsCV[ , paste0('BART_fold_', f, '_pred')] <- as.vector(bartPred_fold)
  predsCV[ , paste0('BART_fold_', f, '_fav')] <- as.vector(bartFav_fold)
  
  message(paste('fold', f, ' - RF'))
  foldTrainingData$presence <- as.factor(foldTrainingData$presence)
  rfFormula <- formula(paste0('presence ~ ', paste(vars, collapse = '+')))
  rfMod_fold <- randomForest(rfFormula, foldTrainingData, na.action = na.exclude)
  rfPred_fold <- predict(rfMod_fold, data4rf, type = 'prob', index = 2)[,2]
  rfFav_fold <- Fav(pred = rfPred_fold, sample.preval = prevalence(model = rfMod_fold))
  predsCV[ , paste0('RF_fold_', f, '_pred')] <- as.vector(rfPred_fold)
  predsCV[ , paste0('RF_fold_', f, '_fav')] <- as.vector(rfFav_fold)
  
  endTime <- proc.time() # end time of the fold
  message('Elapsed fold ', f, ' time: ', (endTime-startTime)[3], ' sec') # elapsed fold time
  message('')
  gc()
}

## save outputs to disk
write.csv(predsCV, file.path('output', folder, 'models/crossValidation/crossPreds.csv'), row.names = FALSE)
######### 


######## Clean memory --------
rm(list = ls())
gc()
