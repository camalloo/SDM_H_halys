######################################################################################
############ Evaluate the cross-validation results and select best models ############
######################################################################################

library(modEvA)
library(tidyverse)

######### Load data ----------
folder <- 'climate_human' # name of the mother folder

predsCV <- read.csv(file.path('output', folder, 'models/crossValidation/crossPreds.csv'))

## model names
# models <- c('BC', 'MXT', 'GLM', 'GAM', 'BART', 'RF')
models <- unique(gsub(paste0('_fold.*'), '', grep('_fold', names(predsCV), value = T)))
## discrimination metrics
discrimination <- c('AUC', 'TSS', 'Kappa', 'Spearman')
## calibration metrics
calibration <- c('Miller', 'Brier', 'Boyce')
######### 


######### Evaluate the cv models ----------
## select the threshold to distinguish predicted presence from absence
threshSelected = 'maxTSS'  # see modEvAmethods("getThreshold") for all thresholds (can also be a value between 0-1 hard coded)

## create df to store evals
folds <- sort(unique(predsCV$fold))
evalResults <- data.frame(fold = folds)
evalResults <- cbind(evalResults, as.data.frame(matrix(nrow = length(folds), 
                                                       ncol = (length(discrimination)+length(calibration))*length(models))))
colnames(evalResults)[-1] <- c(outer(models, c(discrimination, calibration), FUN = paste, sep = '_'))

## evaluate each fold for each model and populate the df
for (f in folds) {
  par(mfrow = c(6, 3), mar = c(4, 2, 2, 2), xpd = FALSE)
  for (m in models) {
    message('Model: ', m, ' - fold: ', f)
    ## get the column to be evaluated for the model-fold combination according to the metric
    columnEval <- grep(paste0(m, '_fold_', f, '_'), names(predsCV), value = T)
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
      ######## Discrimination --------
      ## AUC
      evalResults[f, paste(m, discrimination[1], sep = '_')] <- AUC(obs = toEval[, "presence"], pred = toEval[, columnEval],
                                                                    simplif = T, curve = 'PR', interval = 0.001, 
                                                                    plot = T, main = paste(m, 'AUC'))
      
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
      
      
      ######## Calibration --------
      ## Miller
      evalResults[f, paste(m, calibration[1], sep = '_')] <- MillerCalib(obs = toEval[, "presence"], pred = toEval[, columnEval],
                                                                         main = paste(m, "Miller"))$slope
      
      ## Brier
      evalResults[f, paste(m, calibration[2], sep = '_')] <- errorMeasures(obs = toEval[, "presence"], pred = toEval[, columnEval])[[3]]
      
      ## Boyce
      evalResults[f, paste(m, calibration[3], sep = '_')] <- Boyce(obs = toEval[, "presence"], pred = toEval[, columnEval], 
                                                                   main = paste(m, 'Boyce'))$B
    } else {
      ######## Discrimination --------
      ## AUC
      evalResults[f, paste(m, discrimination[1], sep = '_')] <- AUC(obs = toEval[, "presence"], pred = toEval[, columnEval[fav]],
                                                                    simplif = T, curve = 'PR', interval = 0.001, 
                                                                    plot = T, main = paste(m, 'AUC'))
      
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
      
      
      ######## Calibration --------
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
## save to disk
write.csv(evalResults, file.path('output', folder, 'models/crossValidation/foldsEval_ALL.csv'), row.names = FALSE)
######### 


######## define thresholds to accept models --------
thresh <- c(AUC = 0.8,      # Swets, 1988
            TSS = 0.75,
            Kappa = 0.8,
            Boyce = 0.8,
            Spearman = 0.75,
            Miller = 0.2,    # between 0.7 and 1.3
            Brier = 0.1)
######### 


######## Get average raw metric per model and plot --------
# evalResults <- read.csv(file.path('output', folder, 'models/crossValidation/foldsEval_ALL.csv'))
# ## compute the average performances metricXmodel pf raw performance data
# evalResultsAVG <- evalResults |> 
#   select(-fold) |> 
# evalResultsAVG <- enframe(sapply(evalResultsAVG, mean, na.rm = TRUE), name = 'model_metric', value = 'mean') |> 
#   separate(model_metric, into = c('model', 'metric'), sep = '_', extra = 'merge') |> 
#   mutate(sd = sapply(evalResultsAVG, sd, na.rm = TRUE, USE.NAMES = F)) |> 
#   mutate(mean_sd = paste0(round(mean, 4), ' ± ', round(sd, 2))) |> 
#   select(-mean, -sd) |> 
#   pivot_wider(names_from = model, values_from = mean_sd) |> 
#   column_to_rownames('metric')
# ## check result and save
# evalResultsAVG
# write.csv(evalResultsAVG, file.path('output', folder, 'models/crossValidation/foldsEval_AVG.csv'))
# 
# ## plot
# par(mfrow = c(1, 1), mar = c(7, 3, 2, 1), xpd = F)
# ## boxplot discrimination measures
# boxplot(evalResults[, grepl(paste(discrimination, collapse = '|'), names(evalResults))],
#         col = rep(c('goldenrod1', 'seagreen', 'royalblue', 'tomato', 'cyan3', 'darkorange'), length(discrimination)),
#         las = 2, cex.axis = 0.8, ylim = c(0, 1),
#         main = 'Discrimination metrics')
# vlines <- seq(0.5,
#               ncol(evalResults[, grepl(paste(discrimination, collapse = '|'), names(evalResults))]) + 0.5,
#               by = length(models))
# abline(v = vlines, col = 'darkgrey', lwd = 1.5)
# ## add thresholds to plot
# segments(x0= vlines[-length(vlines)],
#          x1 = vlines[-1],
#          y0 = c(thresh['AUC'], thresh['TSS'], thresh['Kappa'], thresh['Spearman']),
#          lty = 3)
# text(x = seq(from = 0.7, to = ncol(evalResults[, grepl(paste(discrimination, collapse = '|'), names(evalResults))]),
#              by = length(models)), y = 1,
#      labels = discrimination, col = 'darkorchid', cex = 1, font = 2, adj = 0)
# legend('bottomleft', horiz = F, cex = 0.8, x.intersp = 0.2,
#        xpd = T, bty = 'y', bg = 'white',
#        legend = c('thesh', 'bioclim', 'maxent', 'glm', 'gam', 'bart', 'rf'),
#        pch = c(NA, 15, 15, 15, 15, 15, 15),
#        lty = c(3, rep(NA, 6)),
#        col = c('black', 'goldenrod1', 'seagreen', 'royalblue', 'tomato', 'cyan3', 'darkorange'),
#        pt.cex = 1.5)
# 
# ## boxplot calibration measures
# boxplot(evalResults[, grepl(paste(calibration, collapse = '|'), names(evalResults))],
#         col = rep(c('goldenrod1', 'seagreen', 'royalblue', 'tomato', 'cyan3', 'darkorange'),length(calibration)),
#         las = 2, cex.axis = 0.8,
#         main = 'Calibration metrics')
# vlines <- seq(0.5, ncol(evalResults[, grepl(paste(calibration, collapse = '|'), names(evalResults))]) + 0.5, by = length(models))
# abline(v = vlines, col = 'darkgrey', lwd = 1.5)
# 
# ## add thresholds to plot
# segments(x0= vlines[-length(vlines)],
#          x1 = vlines[-1],
#          y0 = c(1, thresh['Brier'], thresh['Boyce']),
#          lty = 3)
# segments(x0 = vlines[1],
#          x1 = vlines[2],
#          y0 = c(1-thresh['Miller'], 1+thresh['Miller']),
#          lty = 3)
# text(x = seq(0.7, ncol(evalResults[, grepl(paste(calibration, collapse = '|'), names(evalResults))]),
#              by = length(models)), y = 2,
#      labels = calibration, col = 'darkorchid', cex = 1, font = 2, adj = 0)
# legend('topright', horiz = F, cex = 0.8, x.intersp = 0.2,
#        xpd = T, bty = 'y', bg = 'white',
#        legend = c('thesh', 'bioclim', 'maxent', 'glm', 'gam', 'bart', 'rf'),
#        pch = c(NA, 15, 15, 15, 15, 15, 15),
#        lty = c(3, rep(NA, 6)),
#        col = c('black', 'goldenrod1', 'seagreen', 'royalblue', 'tomato', 'cyan3', 'darkorange'),
#        pt.cex = 1.5)
# 
# ## clean memory
# gc()
######### 


######## Select the models --------
## 1. normalize the scores: -1 < score < 1
## if score is above metric threshold, normalize between thresh and max value:
##            (score - thresh) / (max - thresh)
## if score is below thresh, normalize using 0 as maximum:
##            (score - thresh) / thresh
evalResults <- read.csv(file.path('output', folder, 'models/crossValidation/foldsEval_ALL.csv'))
## prepare df
ev <- evalResults |> 
  pivot_longer(cols = -fold, names_to = 'mod', values_to = 'score') |> 
  arrange(mod) |> 
  separate(mod, into = c('mod', 'metric'), sep = '_') |> 
  mutate(thresh = thresh[metric],
         normScore = pmax( # force the lower limit to -1
           case_when(
             ## handle Miller who is an intorno of 1
             metric == 'Miller' ~ 1 - abs(score-1)/thresh,
             ## handle Brier where 0 is best result, so pivot it as 1 is best
             metric == 'Brier' & score <= thresh ~ ((1-score) - (1-thresh)) / (1 - (1-thresh)),
             metric == 'Brier' & score > thresh ~ ((1-score) - (1-thresh)) / (1-thresh),
             ## other cases
             score >= thresh ~ (score - thresh) / (1 - thresh),
             score < thresh ~ (score - thresh) / thresh),
           -1))
head(ev)

## plot evaluation metrics
## discrimination metrics
ggplot(data = ev |> 
         filter(metric %in% discrimination),
       aes(x = metric, y = normScore, col = mod)) +
  geom_point(position = position_dodge(width = 1), size = 2) +
  ylim(-1, 1.1) + 
  stat_summary(aes(group = mod), fun = mean, 
               geom = 'point', size = 2, shape = 4, color = 'black',
               position = position_dodge(width = 1)) +
  geom_vline(xintercept = seq(0.5, length(discrimination) + 0.5, by = 1),
             color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = 'Discrimination metrics ~ climate + human') +
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_blank())
## calibration metrics
ggplot(data = ev |> 
         filter(metric %in% calibration),
       aes(x = metric, y = normScore, col = mod)) +
  geom_point(position = position_dodge(width = 1), size = 2) +
  ylim(-1, 1.1) + 
  stat_summary(aes(group = mod), fun = mean, 
               geom = 'point', size = 2, shape = 4, color = 'black',
               position = position_dodge(width = 1)) +
  geom_vline(xintercept = seq(0.5, length(calibration) + 0.5, by = 1),
             color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = 'Calibration metrics ~ climate + human') +
  theme(axis.text.x = element_text(size = 10), axis.title.x = element_blank())

## 2. Investigate if differences in performance among models are significant
## investigate the data to check which test to perform:
## a. check normality of data
qqnorm(ev$normScore) # Q-Q plot
qqline(ev$normScore)
hist(ev$normScore) # histogram
shapiro.test(ev$normScore) # Shapiro–Wilk test (normal if p.value > 0.05)

## DATA ARE NOT NORMALLY DISTRIBUTED

## to get differences perform Friedman rank test on each metric (use folds as repeated observations)
## if significant difference perform Nemenyi test
testRes <- data.frame()
postHoc <- list()
## run the test
for (m in c(discrimination, calibration)) {
  ## filter for the metric
  df <- ev |> 
    filter(metric == m)
  ## run Friedman rank test
  t <- friedman.test(normScore ~ mod | fold, data = df)
  ## store result
  res <- data.frame(metric = m,
                    chi_sq = as.numeric(t$statistic),
                    p_value = as.numeric(t$p.value),
                    sig = ifelse(t$p.value < 0.05, 'H1', 'H0'))
  if (res$sig == 'H1') {
    cat('There is significant difference in', m, '\n', '...performing Nemenyi test\n')
    cat('\n')
    p <- PMCMRplus::frdAllPairsNemenyiTest(normScore ~ mod | fold, data = df)
    postHoc[[m]] <- data.frame(cbind(p$p.value))
  } else {
    cat('There is no significant difference in', m, '\n', '...no action')
    cat('\n')
  }
  ## populate output df
  testRes <- rbind(testRes, res)
}

print(postHoc)
names(postHoc)

## since there are significant differences in models performance
## get the average ranking of the models in each metric and overall rank
ranks <- ev |>
  # filter(metric %in% names(postHoc)) |> # only use metrics where difference is significant to create rank
  group_by(fold, metric) |> 
  mutate(rank = rank(-normScore)) |> 
  group_by(mod, metric) |> 
  summarise(rank = mean(rank), .groups = 'drop') |> 
  pivot_wider(names_from = metric, values_from = rank) |> 
  mutate(overall_mean = rowMeans(across(-mod)))

## keep only the models that passed all three thresholds
selMods <- ranks |> 
  arrange(overall_mean) |> 
  slice(1:3)
selMods
write.csv(selMods, file.path('output', folder, 'models/crossValidation/selectedModels.csv'), row.names = F)
######### 


######## Clean memory --------
rm(list = ls())
gc()
