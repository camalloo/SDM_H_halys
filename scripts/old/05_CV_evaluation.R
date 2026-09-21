library("modEvA")


################# Load Data_________________________________
folder <- 'AllPredictors' # name of the mother folder
predsCV <- read.csv(paste0('output/', folder, '/models/crossValidation/crossPreds.csv'))





################# Evaluate the cv models_________________________________
######## get the names of the models and define the evaluation measures --------
models <- unique(gsub(paste0('_fold.*'), '', grep('_fold', names(predsCV), value = T)))
## discrimination metrics
discrimination <- c('AUC', 'TSS', 'Kappa', 'Spearman')
## calibration metrics
calibration <- c('Miller', 'Brier', 'Boyce')
## select the threshold to distinguish predicted presence from absence
threshSelected = 'maxTSS'  # see modEvAmethods("getThreshold") for all thresholds ( can also be a value between 0-1 hard coded)


######## run the cross evaluation and save the results --------
## create df to store evals
folds <- sort(unique(predsCV$fold))
evalResults <- as.data.frame(matrix(nrow = length(folds), ncol = (length(discrimination)+length(calibration))*length(models)))
colnames(evalResults) <- c(outer(models, c(discrimination, calibration), FUN = paste, sep = '_'))

## evaluate each fold for each model and populate the df
for (f in folds) {
  par(mfrow = c(6, 3), mar = c(4, 2, 2, 2), xpd = FALSE)
  for (m in models) {
    message('Model: ', m, ' - fold: ', f)
    ## get the column to be evaluated for the model-fold combination according to the metric
    columnEval <- grep(paste0(m, '_fold_', f), names(predsCV), value = T)
    if (length(columnEval) > 1) {
      if (m == 'BC') {
        ## use the continuous reclassification for Bioclim
        columnEval <- columnEval[grep("_reclass", columnEval)]
      } else {
        ## separate favourability and probability for BART, GLM, GAM, and RF
        fav <- grep("_fav", columnEval)
        pred <- grep("_pred", columnEval)
      }
    }
    ## subset the data to the one you want to evaluate and paste the result in the right column
    toEval <- subset(predsCV, fold == f)
    if (m %in% c('BC', 'MXT')) {
      ######## discrimination --------
      ## AUC
      evalResults[f, paste(m, discrimination[1], sep = '_')] <- AUC(obs = toEval[, "presence"], pred = toEval[, columnEval],
                                                                 simplif = T, plot = T, main = paste(m, 'AUC'))
      
      ## TSS
      evalResults[f, paste(m, discrimination[2], sep = '_')] <- threshMeasures(obs = toEval[, "presence"], pred = toEval[, columnEval],
                                                                     thresh = threshSelected, measures = 'TSS',
                                                                     simplif = TRUE, plot = F, standardize = FALSE)
      
      ## Cohen's Kappa
      evalResults[f, paste(m, discrimination[3], sep = '_')] <- threshMeasures(obs = toEval[, "presence"], pred = toEval[, columnEval],
                                                                     thresh = threshSelected, measures = 'kappa',
                                                                     simplif = TRUE, plot = F, standardize = FALSE)[1]
      
      ## Spearman Rank
      evalResults[f, paste(m, discrimination[4], sep = '_')] <- cor.test(x = toEval[, "presence"], y = toEval[, columnEval], method = 'spearman')$estimate
      
      
      ######## calibration --------
      ## Miller
      evalResults[f, paste(m, calibration[1], sep = '_')] <- MillerCalib(obs = toEval[, "presence"], pred = toEval[, columnEval],
                                                                  main = paste(m, "Miller"))$slope
      
      ## Brier
      evalResults[f, paste(m, calibration[2], sep = '_')] <- errorMeasures(obs = toEval[, "presence"], pred = toEval[, columnEval])[[3]]
      
      ## Boyce
      evalResults[f, paste(m, calibration[3], sep = '_')] <- Boyce(obs = toEval[, "presence"], pred = toEval[, columnEval], 
                                                                   main = paste(m, 'Boyce'))$B
    } else {
      ######## discrimination --------
      ## AUC
      evalResults[f, paste(m, discrimination[1], sep = '_')] <- AUC(obs = toEval[, "presence"], pred = toEval[, columnEval[fav]],
                                                                 simplif = T, plot = T, main = paste(m, 'AUC'))
      
      ## TSS
      evalResults[f, paste(m, discrimination[2], sep = '_')] <- threshMeasures(obs = toEval[, "presence"], pred = toEval[, columnEval[fav]],
                                                                     thresh = threshSelected, measures = 'TSS',
                                                                     simplif = TRUE, plot = F, standardize = FALSE)
      
      ## Cohen's Kappa
      evalResults[f, paste(m, discrimination[3], sep = '_')] <- threshMeasures(obs = toEval[, "presence"], pred = toEval[, columnEval[fav]],
                                                                     thresh = threshSelected, measures = 'kappa',
                                                                     simplif = TRUE, plot = F, standardize = FALSE)[1]
      
      ## Spearman Rank
      evalResults[f, paste(m, discrimination[4], sep = '_')] <- cor.test(x = toEval[, "presence"], y = toEval[, columnEval[fav]], method = 'spearman')$estimate
      
      
      ######## calibration --------
      ## Miller
      evalResults[f, paste(m, calibration[1], sep = '_')] <- MillerCalib(obs = toEval[, "presence"], pred = toEval[, columnEval[fav]],
                                                                  main = paste(m, "Miller"))$slope
      
      ## Brier
      evalResults[f, paste(m, calibration[2], sep = '_')] <- errorMeasures(obs = toEval[, "presence"], pred = toEval[, columnEval[pred]])[[3]]
      
      ## Boyce
      evalResults[f, paste(m, calibration[3], sep = '_')] <- Boyce(obs = toEval[, "presence"], pred = toEval[, columnEval[fav]], 
                                                                   main = paste(m, 'Boyce'))$B
    }
  }
  par(mfrow = c(1, 1), xpd = T)
  mtext(paste('Fold:', f), side = 1, outer = T, xpd = T, line = -1) 
}


######## plots to better understand metrics and save --------
par(mfrow = c(1, 1), mar = c(7, 3, 2, 1), xpd = F)
## boxplot discrimination measures
boxplot(evalResults[, grepl(paste(discrimination, collapse = '|'), names(evalResults))], 
        col = rep(c('goldenrod1', 'seagreen', 'royalblue', 'tomato', 'cyan3', 'darkorange'), length(discrimination)), 
        las = 2, cex.axis = 0.8, 
        main = 'Discrimination metrics')
abline(h = 1, lty = 2)
abline(v = seq(0.5, ncol(evalResults[, grepl(paste(discrimination, collapse = '|'), names(evalResults))]) + 0.5, by = length(models)), 
       col = 'darkgrey', lwd = 1.5)
text(x = seq(from = 0.7, to = ncol(evalResults[, grepl(paste(discrimination, collapse = '|'), names(evalResults))]), 
             by = length(models)), y = max(evalResults[, grepl(paste(discrimination, collapse = '|'), names(evalResults))], na.rm = TRUE), 
     labels = discrimination, col = 'darkorchid', cex = 1, font = 2, adj = 0)
legend('bottomleft', horiz = F, 
       xpd = T, bty = 'y', bg = 'white',
       legend = c('bioclim', 'maxent', 'glm', 'gam', 'bart', 'rf'), 
       fill = c('goldenrod1', 'seagreen', 'royalblue', 'tomato', 'cyan3', 'darkorange'))

## boxplot calibration measures
boxplot(evalResults[, grepl(paste(calibration, collapse = '|'), names(evalResults))], 
        col = rep(c('goldenrod1', 'seagreen', 'royalblue', 'tomato', 'cyan3', 'darkorange'),length(calibration)), 
        las = 2, cex.axis = 0.8, 
        main = 'Calibration metrics')
abline(h = 1, lty = 2)
abline(v = seq(0.5, ncol(evalResults[, grepl(paste(calibration, collapse = '|'), names(evalResults))]) + 0.5, by = length(models)), 
       col = 'darkgrey', lwd = 1.5)
text(x = seq(0.7, ncol(evalResults[, grepl(paste(calibration, collapse = '|'), names(evalResults))]), 
             by = length(models)), y = max(evalResults, na.rm = TRUE), 
     labels = calibration, col = 'darkorchid', cex = 1, font = 2, adj = 0)
legend('bottomright', horiz = F, 
       xpd = T, bty = 'y', bg = 'white',
       legend = c('bioclim', 'maxent', 'glm', 'gam', 'bart', 'rf'), 
       fill = c('goldenrod1', 'seagreen', 'royalblue', 'tomato', 'cyan3', 'darkorange'))

write.csv(evalResults, paste0('output/', folder, '/models/crossValidation/crossModelsEval.csv'), row.names = FALSE)

## clean unused memory
gc()





################# Select the models_________________________________
## models selection is based on the average performance of the block cv runs
## selected models will be run using the complete dataset for training
evalResults <- read.csv(paste0('output/', folder, '/models/crossValidation/crossModelsEval.csv'))

#### compute the average performances
evalResultsAVG <- sapply(evalResults, mean, na.rm = TRUE)
evalResultsAVG

######## define thresholds for accepting a model and select the models --------
AUC = 0.75      # Swets, 1988
TSS = 0.6
Kappa = 0.4
Boyce = 0.6
Spearman = 0.4
Miller = 0.3    # between 0.3 and 1.3
Brier = 0.2

for (m in c(discrimination, calibration)) {
  if (m == 'AUC') {
    avgAUC <- evalResultsAVG[grep(m, names(evalResultsAVG), value = T)]
    print(avgAUC)
    avgAUC <- avgAUC[which(avgAUC >= AUC)]
    names(avgAUC) <- sub('_AUC', '', names(avgAUC))
    print(avgAUC)
    cat('\n')
  } else if (m == 'TSS') {
    avgTSS <- evalResultsAVG[grep(m, names(evalResultsAVG), value = T)]
    print(avgTSS)
    avgTSS <- avgTSS[which(avgTSS >= TSS)]
    names(avgTSS) <- sub('_TSS', '', names(avgTSS))
    print(avgTSS)
    cat('\n')
  } else if (m == 'Kappa') {
    avgKappa <- evalResultsAVG[grep(m, names(evalResultsAVG), value = T)]
    print(avgKappa)
    avgKappa <- avgKappa[which(avgKappa >= Kappa)]
    names(avgKappa) <- sub('_Kappa', '', names(avgKappa))
    print(avgKappa)
    cat('\n')
  } else if (m == 'Spearman') {
    avgSpearman <- evalResultsAVG[grep(m, names(evalResultsAVG), value = T)]
    print(avgSpearman)
    avgSpearman <- avgSpearman[which(avgSpearman >= Spearman)]
    names(avgSpearman) <- sub('_Spearman', '', names(avgSpearman))
    print(avgSpearman)
    cat('\n')
  } else if (m == 'Boyce') {
    avgBoyce <- evalResultsAVG[grep(m, names(evalResultsAVG), value = T)]
    print(avgBoyce)
    avgBoyce <- avgBoyce[which(avgBoyce >= Boyce)]
    names(avgBoyce) <- sub('_Boyce', '', names(avgBoyce))
    print(avgBoyce)
    cat('\n')
  } else if (m == 'Brier') {
    avgBrier <- evalResultsAVG[grep(m, names(evalResultsAVG), value = T)]
    print(avgBrier)
    avgBrier <- avgBrier[which(avgBrier <= Brier)]
    names(avgBrier) <- sub('_Brier', '', names(avgBrier))
    print(avgBrier)
    cat('\n')
  } else if (m == 'Miller') {
    avgMiller <- evalResultsAVG[grep(m, names(evalResultsAVG), value = T)]
    print(avgMiller)
    avgMiller <- avgMiller[which(avgMiller >= Miller & avgMiller <= 1+Miller)]
    names(avgMiller) <- sub('_Miller', '', names(avgMiller))
    print(avgMiller)
    cat('\n')
  } 
}

## keep only the models that passed all three thresholds
selMods <- Reduce(intersect, list(names(avgAUC), names(avgTSS), 
                                  names(avgKappa), names(avgSpearman), 
                                  names(avgBoyce), names(avgBrier), 
                                  names(avgMiller)))
selMods
write.csv(selMods, paste0('output/', folder, '/models/crossValidation/selectedModels.csv'), row.names = F)





################# Clean memory_________________________________
rm(list = ls())
gc()
