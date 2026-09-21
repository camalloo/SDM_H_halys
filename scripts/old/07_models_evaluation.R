library("predicts")
library("maxnet")
library("mgcv")
library("randomForest")
library("pROC")


################# Load data_________________________________
folder <- 'world' # name of the mother folder

preds <- read.csv(paste0('output/', folder, '/models/allModelsPredictions.csv'))
data <- read.csv(paste0('output/', folder, '/modellingData/modData.csv'))
vars <- read.csv(paste0('output/', folder, '/modellingData/selectedVariables.csv'))[ , 1] # variables after collinearity analysis
## load the models
glmMod <- readRDS(paste0('output/', folder, '/models/glm_mod.rds'))
gamMod <- readRDS(paste0('output/', folder, '/models/gam_mod.rds'))
## BART MODEL cannot be saved so always re-fit it if not loaded on env (KEEP THE SEED)
bartMod <- dbarts::bart(x.train = data[, vars], y.train = data[, 'presence'],
                        ntree = 100, keeptrees = T, verbose = F, seed = 1312)
rfMod <- readRDS(paste0('output/', folder, '/models/rf_mod.rds'))





################# Correlation and PCA of predictions_________________________________
## select only comparable metrics (for glm and bart use favourability)
names(preds)
cols <- c('gamFav', 'bartFav', 'rfFav')


######## Correlation --------
correlation <- cor(preds[, cols])
## plot
corrplot::corrplot(correlation, method = 'ellipse', type = 'upper',
         addCoef.col = 'wheat3', addCoefasPercent = F,
         tl.col = 'black', diag = T)


######## PCA of models predictions --------
pca <- prcomp(preds[, cols])
summary(pca)
## plot
plot(pca)
biplot(pca, xlabs = rep('', nrow(pca$x)))





################# Variable importance and response curve_________________________________
## merge the predictions dataset with the input dataset
predsData <- merge(data, preds, all = T)
head(predsData)
## check the names of the bioclim used for modelling
vars == names(glmMod$levels) # same for all models


######## Variable importance --------
## through permutations (model-independent) for mxt, gam, rf
set.seed(1312)
## maxent
GLMvarimp <- varImportance(glmMod, y = predsData$presence, x = predsData[, vars], n = 50, stat = 'RMSE')
save(GLMvarimp, file = paste0('output/', folder, '/models/evaluation/glmVarImp.RData'))
# load(file = paste0('output/', folder, '/models/evaluation/glmVarImp.RData'))

## gam
GAMvarimp <- varImportance(gamMod, y = predsData$presence, x = predsData[, vars], n = 50, stat = 'RMSE')
save(GAMvarimp, file = paste0('output/', folder, '/models/evaluation/gamVarImp.RData'))
# load(file = paste0('output/', folder, '/models/evaluation/gamVarImp.RData'))

## rf
RFvarimp <- varImportance(rfMod, y = predsData$presence, x = predsData[, vars], n = 50, stat = 'RMSE', type = 'prob', index = 2)
save(RFvarimp, file = paste0('output/', folder, '/models/evaluation/rfVarImp.RData'))
# load(file = paste0('output/', folder, '/models/evaluation/rfVarImp.RData'))

## plots
barplot(sort(GLMvarimp), horiz = TRUE, las = 2, main = 'GLM')
barplot(sort(GAMvarimp), horiz = TRUE, las = 2, main = "GAM")
barplot(sort(RFvarimp), horiz = TRUE, las = 2, main = "RF")

## variable appearances in trees for bart
BARTvarimp <- embarcadero::varimp(bartMod, plots = T)
# save(bartVarImp, file = paste0('output/', folder, '/models/evaluation/bartVarImp.RData'))
# load(file = paste0('output/', folder, '/models/evaluation/bartVarImp.RData'))


######## Response curves --------
## maxent
plot(glmMod, pages = 1, shade = F)

## gam
plot(gamMod, pages = 1, shade = F)

## bart
bartRCs <- partial(bartMod, trace = T)
save(bartRCs, file = paste0('output/', folder, '/models/evaluation/bartResponsecurves.RData'))
load(file = paste0('output/', folder, '/models/evaluation/bartResponsecurves.RData'))
gridExtra::grid.arrange(grobs = bartRCs, nrow = 2, ncol = 3)

## rf
rfRCs <- list()
for (var in vars) {
  thisRC <- do.call(randomForest::partialPlot,
                    list(x = rfMod, pred.data = data, x.var = var, plot = F))
  rfRCs[[var]] <- thisRC
}
save(rfRCs, file = paste0('output/', folder, '/models/evaluation/rfResponsecurves.RData'))
load(file = paste0('output/', folder, '/models/evaluation/rfResponsecurves.RData'))
for (p in 1:length(rfRCs)) {
  thisRC <- rfRCs[[p]]
  plot(x = thisRC$x, y = thisRC$y, type = 'l',
       xlab = '', ylab = 'Response', main = names(rfRCs)[p])
}






################# Evaluation metrics_________________________________
par(mfrow = c(2, 2))
######## AUC --------
## maxent
with(preds, 
     roc(response = presence, predictor = mxtPred, plot = T, legacy.axes = T, 
         main = 'MXT',
         xlab = 'False Positive Rate (1-specificity)',
         ylab = 'True Positive Rate (sensitivity)',
         col = 'blue', lwd = 2, print.auc = T)
)

## gam
with(preds, 
     roc(response = presence, predictor = gamPred, plot = T, legacy.axes = T, 
         main = 'GAM',
         xlab = 'False Positive Rate (1-specificity)',
         ylab = 'True Positive Rate (sensitivity)',
         col = 'blue', lwd = 2, print.auc = T)
)

## bart
with(preds,
     roc(response = presence, predictor = bartPred, plot = T, legacy.axes = T, 
         main = 'BART',
         xlab = 'False Positive Rate (1-specificity)',
         ylab = 'True Positive Rate (sensitivity)',
         col = 'blue', lwd = 2, print.auc = T)
)

## rf
with(preds, 
     roc(response = presence, predictor = rfPred, plot = T, legacy.axes = T, 
         main = 'RF',
         xlab = 'False Positive Rate (1-specificity)',
         ylab = 'True Positive Rate (sensitivity)',
         col = 'blue', lwd = 2, print.auc = T)
)

## save AUC values to df
AUCs <- data.frame(mxtAUC = auc(response = preds$presence, predictor = preds$mxtPred),
                   gamAUC = auc(response = preds$presence, predictor = preds$gamPred),
                   bartAUC = auc(response = preds$presence, predictor = preds$bartPred),
                   rfAUC = auc(response = preds$presence, predictor = preds$rfPred))
write.csv(AUCs, paste0('output/', folder, '/models/evaluation/AUCall.csv'), row.names = F)





################# Clean memory_________________________________
rm(list = ls())
gc()
