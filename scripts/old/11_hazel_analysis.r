library("terra")


################# Load data_________________________________
folder <- 'world' # name of the mother folder

world <- geodata::world(path = 'data/geodata/', resolution = 1)
data <- read.csv(paste0('output/', folder, '/modellingData/modData.csv')) # train dataset
region <- vect(paste0('output/', folder, '/modellingData/ModelRegion/modellingRegionPoly.shp')) # modelling region polygon
vars <- read.csv(paste0('output/', folder, '/modellingData/selectedVariables.csv'))[, 1] # variables after collinearity analysis

## current climate ensemble
ensProj <- rast(paste0('output/', folder, '/projections/rast/ensemble_proj.tif'))
ensProj$variance <- rast(paste0('output/', folder, '/projections/rast/ensemble_variance.tif'))
## plot to check
par(oma=c(1,1,3,1))
plot(ensProj, nc = 1, nr = 2)
title(main = '1991-2010 ensemble of mxt, gam, and bart models', outer = T, cex.main = 1.5, line = 0)


######## each sdm separately (five future climate models ensemble) --------
ssps <- c('ssp126', 'ssp370', 'ssp585')
par(oma=c(1,2,3,2))
## maxent
mxt2011_2040 <- rast(paste0('output/', folder, '/projections/Future/ensembles/2011_2040_all_ssps_mxtFutEns.tif'))
mxt2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/2041_2070_all_ssps_mxtFutEns.tif'))
## plot 2011-2040
plot(mxt2011_2040, nc = 2, nr = 3)
at = 5/6
for (l in 1:length(ssps)) {
  print(ssps[l])
  mtext(ssps[l], side = 2, outer = T, at = at)
  at = at - 2/6
}
title(main = '2011_2040 mxt ensemble of five climate models', outer = T, cex.main = 1.5, line = 0)
## plot 2041-2070
plot(mxt2041_2070, nc = 2, nr = 3)
at = 5/6
for (l in 1:length(ssps)) {
  print(ssps[l])
  mtext(ssps[l], side = 2, outer = T, at = at)
  at = at - 2/6
}
title(main = '2041_2070 mxt ensemble of five climate models', outer = T, cex.main = 1.5, line = 0)

## gam
gam2011_2040 <- rast(paste0('output/', folder, '/projections/Future/ensembles/2011_2040_all_ssps_gamFutEns.tif'))
gam2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/2041_2070_all_ssps_gamFutEns.tif'))
## plot 2011-2040
plot(gam2011_2040, nc = 2, nr = 3)
at = 5/6
for (l in 1:length(ssps)) {
  print(ssps[l])
  mtext(ssps[l], side = 2, outer = T, at = at)
  at = at - 2/6
}
title(main = '2011_2040 gam ensemble of five climate models', outer = T, cex.main = 1.5, line = 0)
## plot 2041-2070
plot(gam2041_2070, nc = 2, nr = 3)
at = 5/6
for (l in 1:length(ssps)) {
  print(ssps[l])
  mtext(ssps[l], side = 2, outer = T, at = at)
  at = at - 2/6
}
title(main = '2041_2070 gam ensemble of five climate models', outer = T, cex.main = 1.5, line = 0)

## bart
bart2011_2040 <- rast(paste0('output/', folder, '/projections/Future/ensembles/2011_2040_all_ssps_bartFutEns.tif'))
bart2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/2041_2070_all_ssps_bartFutEns.tif'))
## plot 2011-2040
plot(bart2011_2040, nc = 2, nr = 3)
at = 5/6
for (l in 1:length(ssps)) {
  print(ssps[l])
  mtext(ssps[l], side = 2, outer = T, at = at)
  at = at - 2/6
}
title(main = '2011_2040 bart ensemble of five climate models', outer = T, cex.main = 1.5, line = 0)
## plot 2041-2070
plot(bart2041_2070, nc = 2, nr = 3)
at = 5/6
for (l in 1:length(ssps)) {
  print(ssps[l])
  mtext(ssps[l], side = 2, outer = T, at = at)
  at = at - 2/6
}
title(main = '2041_2070 bart ensemble of five climate models', outer = T, cex.main = 1.5, line = 0)


######## sdm ensembles --------
par(oma=c(1,1,1,1))
## 2011-2040
ssp126_2011_2040 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2011_2040_ssp126Ens.tif'))
plot(ssp126_2011_2040, nc = 1, nr = 2)
title(main = 'ssp126_2011_2040 ensemble of mxt, gam, and bart models', outer = T, cex.main = 2, line = -1)

ssp370_2011_2040 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2011_2040_ssp370Ens.tif'))
plot(ssp370_2011_2040, nc = 1, nr = 2)
title(main = 'ssp370_2011_2040 ensemble of mxt, gam, and bart models', outer = T, cex.main = 2, line = -1)

ssp585_2011_2040 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2011_2040_ssp585Ens.tif'))
plot(ssp585_2011_2040, nc = 1, nr = 2)
title(main = 'ssp585_2011_2040 ensemble of mxt, gam, and bart models', outer = T, cex.main = 2, line = -1)


# ## 2041-2070
# ssp126_2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2041_2070_ssp126Ens.tif'))
# plot(ssp126_2041_2070, nc = 1, nr = 2)
# title(main = 'ssp126_2041_2070 ensemble of mxt, gam, and bart models', outer = T, cex.main = 1.5, line = 0)
# 
# ssp370_2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2041_2070_ssp370Ens.tif'))
# plot(ssp370_2041_2070, nc = 1, nr = 2)
# title(main = 'ssp370_2041_2070 ensemble of mxt, gam, and bart models', outer = T, cex.main = 1.5, line = 0)
# 
# ssp585_2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2041_2070_ssp585Ens.tif'))
# plot(ssp585_2041_2070, nc = 1, nr = 2)
# title(main = 'ssp585_2041_2070 ensemble of mxt, gam, and bart models', outer = T, cex.main = 1.5, line = 0)





################# MESS analysis_________________________________
## method from Elith et al., 2010: compute the dissimilarity between the variables used in training and the ones used for projection

######## current climate --------
## load the bioclims
bioCurr <- rast('data/CHELSAbioclim/historicCropped/bio.tif')
## compute mess
messCurrent <- predicts::mess(x = bioCurr[[vars]], v = data[, vars])
## reproject to only display negative values
messCurrent <- ifel(messCurrent >= 0, NA, messCurrent)
plot(messCurrent, col = heat.colors(10), main = "MESS")
plot(world, add = T)

######## future climate --------
## make one mess for each climate model and then average
## set up the counters for the loops
periods <- c('2011_2040')#, '2041_2070')
scenarios <- c('ssp126')#, 'ssp370', 'ssp585')
for (period in periods) {
  cat(paste0('Starting period ', period, '...\n\n'))
  ## second loop over the ssp scenarios
  for (s in scenarios) {
    cat(paste0('Starting scenario ', s, '...\n'))
    ## list the five climate modles
    climateModels <- list.files(paste0('output/', folder, '/modellingData/Future/', 
                                       period, '/', s), pattern = '.tif')
    all <- NULL # placeholder to store the mess for each model
    ## third loop over the five models
    for (m in climateModels) {
      print(basename(m))
      ## load bioclims
      bioFut <- rast(paste0('output/', folder, '/modellingData/Future/',
                            period, '/', s, '/', basename(m)))
      ## compute mess
      mess <- predicts::mess(x = bioFut, v = data[, vars])
      ## add mess(m) to placeholder
      if (is.null(all)) {all <- mess} else {all <- c(all, mess)}
    }
    avg <- app(all, fun = base::mean)
    avg <- ifel(avg >= 0, NA, avg)
    plot(avg, col = heat.colors(10), main = paste('MESS', period, s))
    plot(world, add = T)
  }
}





################# Crop maps to hazelnut regions_________________________________
## usa
usa <- crop(world[world$GID_0 == 'USA'], ext(-130, -60, 0, 50))
## current climate
usa1991_2010 <- crop(ensProj, usa, mask = T)
plot(usa1991_2010, nc = 1, nr = 2)
title(main = '1991-2010 USA ensemble', outer = T, cex.main = 1.5, line = 0)
## future climate
usa2011_2040ssp126 <- crop(ssp126_2011_2040, usa, mask = T)
plot(usa2011_2040ssp126, nc = 1, nr = 2)
title(main = '2011-2040 ssp126 USA ensemble', outer = T, cex.main = 1.5, line = 0)

## italy
ita <- world[world$GID_0 == 'ITA']
## current climate
ita1991_2010 <- crop(ensProj, ita, mask = T)
plot(ita1991_2010, nc = 2, nr = 1)
title(main = '1991-2010 ITA ensemble', outer = T, cex.main = 1.5, line = 0)
## future climate
ita2011_2040ssp126 <- crop(ssp126_2011_2040, ita, mask = T)
plot(ita2011_2040ssp126, nc = 2, nr = 1)
title(main = '2011-2040 ssp126 ITA ensemble', outer = T, cex.main = 1.5, line = 0)

## france
fra <- world[world$GID_0 == 'FRA']
## current climate
fra1991_2010 <- crop(ensProj, fra, mask = T)
plot(fra1991_2010, nc = 2, nr = 1)
title(main = '1991-2010 FRA ensemble', outer = T, cex.main = 1.5, line = 0)
## future climate
fra2011_2040ssp126 <- crop(ssp126_2011_2040, fra, mask = T)
plot(fra2011_2040ssp126, nc = 2, nr = 1)
title(main = '2011-2040 ssp126 FRA ensemble', outer = T, cex.main = 1.5, line = 0)

## serbia
srb <- world[world$GID_0 == 'SRB']
## current climate
srb1991_2010 <- crop(ensProj, srb, mask = T)
plot(srb1991_2010, nc = 2, nr = 1)
title(main = '1991-2010 SRB ensemble', outer = T, cex.main = 1.5, line = 0)
## future climate
srb2011_2040ssp126 <- crop(ssp126_2011_2040, srb, mask = T)
plot(srb2011_2040ssp126, nc = 2, nr = 1)
title(main = '2011-2040 ssp126 SRB ensemble', outer = T, cex.main = 1.5, line = 0)

## turkey, georgia, and georgia
turazegeo <- world[world$GID_0 %in% c('GEO', 'AZE', 'TUR')]
## current climate
turazegeo1991_2010 <- crop(ensProj, turazegeo, mask = T)
plot(turazegeo1991_2010, nc = 1, nr = 2)
title(main = '1991-2010 GEO, AZE, TUR ensemble', outer = T, cex.main = 1.5, line = 0)
## future climate
turazegeo2011_2040ssp126 <- crop(ssp126_2011_2040, turazegeo, mask = T)
plot(turazegeo2011_2040ssp126, nc = 1, nr = 2)
title(main = '2011-2040 ssp126 GEO, AZE, TUR ensemble', outer = T, cex.main = 1.5, line = 0)

## chile and argentina
chiarg <- crop(world[world$GID_0 %in% c('ARG', 'CHL')], ext(-78, -40, -69, -10))
## current climate
chiarg1991_2010 <- crop(ensProj, chiarg, mask = T)
plot(chiarg1991_2010, nc = 2, nr = 1)
title(main = '1991-2010 CHI and ARG ensemble', outer = T, cex.main = 1.5, line = 0)
## future climate
chiarg2011_2040ssp126 <- crop(ssp126_2011_2040, chiarg, mask = T)
plot(chiarg2011_2040ssp126, nc = 2, nr = 1)
title(main = '2011-2040 ssp126 CHI and ARG ensemble', outer = T, cex.main = 1.5, line = 0)


