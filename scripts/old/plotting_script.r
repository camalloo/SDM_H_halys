library("terra")


################# Load data_________________________________
folder <- 'world' # name of the mother folder
world <- geodata::world(path = 'data/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica
## study area
mxtPred <- rast(paste0('output/', folder, '/models/rast/mxt_predStudyArea.tif'))
gamPred <- rast(paste0('output/', folder, '/models/rast/gam_predStudyArea.tif'))
gamFav <- rast(paste0('output/', folder, '/models/rast/gam_favStudyArea.tif'))
bartPred <- rast(paste0('output/', folder, '/models/rast/bart_predStudyArea.tif'))
bartFav <- rast(paste0('output/', folder, '/models/rast/bart_favStudyArea.tif'))
rfPred <- rast(paste0('output/', folder, '/models/rast/rf_predStudyArea.tif'))
rfFav <- rast(paste0('output/', folder, '/models/rast/rf_favStudyArea.tif'))

par(mfrow = c(2,2))
plot(bartPred, buffer = T)
plot(world, lwd = 0.5, add = T)
plot(bartFav, buffer = T, fun = function() {
  plot(world, lwd = 0.5, add = TRUE, border = "black")
})
title('bart favourability in the study area', outer = T, line = -1)
plot(bartFav, buffer = T)
plot(world, lwd = 0.5, add = T)

## world projections
mxtProj <- rast(paste0('output/', folder, '/projections/rast/mxt_proj.tif'))
gamProj <- rast(paste0('output/', folder, '/projections/rast/gam_proj.tif'))
gamProjFav <- rast(paste0('output/', folder, '/projections/rast/gam_projFav.tif'))
bartProj <- rast(paste0('output/', folder, '/projections/rast/bart_proj.tif'))
bartProjFav <- rast(paste0('output/', folder, '/projections/rast/bart_projFav.tif'))

allProj <- c(mxtProj, gamProjFav, bartProjFav[[1]])
allProj
names(allProj) <- sub('.pred|.fav', '', names(allProj))
plot(allProj, range = c(0, 1))
title(main = 'Model projections world', outer = T, line = -1)

## ensemble
ensProj <- rast(paste0('output/', folder, '/projections/rast/ensemble_proj.tif'))
ensVar <- rast(paste0('output/', folder, '/projections/rast/ensemble_variance.tif'))

par(mfrow = c(2, 1))
plot(ensProj, main = ' Ensemble mean weighted on AUC', range = c(0, 1))
plot(ensVar, main = ' Ensemble variance', range = c(0, 1))


## future climate
sspsMXT <- rast(paste0('output/', folder, '/projections/Future/ensembles/2011_2040_all_ssps_mxtFutEns.tif'))
sspsGAM <- rast(paste0('output/', folder, '/projections/Future/ensembles/2011_2040_all_ssps_gamFutEns.tif'))
sspsBART <- rast(paste0('output/', folder, '/projections/Future/ensembles/2011_2040_all_ssps_bartFutEns.tif'))
sspsMXT2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/2041_2070_all_ssps_mxtFutEns.tif'))
sspsGAM2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/2041_2070_all_ssps_gamFutEns.tif'))
sspsBART2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/2041_2070_all_ssps_bartFutEns.tif'))

ssps = sspsMXT
plot(ssps, nc = 2, nr = 3)
at = 5/6
for (lyr in 1:nlyr(ssps)) {
  if (lyr %% 2 == 0) {
    mtext(sub(".*(ssp[0-9]{3}).*", "\\1", basename(sources(ssps[[lyr]]))), side = 2, outer = T, at = at)
    at = at - 2/6
  }
}
title(main = paste(period, model), outer = T, cex.main = 2.5, line = -1)




## hazelnut 1
par(mfrow = c(3,2), oma = c(1,1,3,1))
plot(chiarg1991_2010[[1]], main="Chile and Argentina (1991-2010)", buffer = T)
plot(world, lwd = 0.5, add = T)
plot(chiarg2011_2040ssp126[[1]], main="Chile and Argentina (2011-2040 ssp126)", buffer = T)
plot(world, lwd = 0.5, add = T)

plot(ita1991_2010[[1]], main="Italy (1991-2010)", buffer = T)
plot(world, lwd = 0.5, add = T)
plot(ita2011_2040ssp126[[1]], main="Italy (2011-2040 ssp126)", buffer = T)
plot(world, lwd = 0.5, add = T)

plot(srb1991_2010[[1]], main="Serbia (1991-2010)", buffer = T)
plot(world, lwd = 0.5, add = T)
plot(srb2011_2040ssp126[[1]], main="Serbia (2011-2040 ssp126)", buffer = T)
plot(world, lwd = 0.5, add = T)
par(mfrow = c(1,1))
title('Hazelnut sourcing countries', outer = T, line = 1)

## hazelnut 2
par(mfrow = c(3,2), oma = c(1,1,3,1))
plot(turazegeo1991_2010[[1]], main="Turkey, Georgia, and Aerbaijan (1991-2010)", buffer = T)
plot(world, lwd = 0.5, add = T)
plot(turazegeo2011_2040ssp126[[1]], main="Turkey, Georgia, and Aerbaijan (2011-2040 ssp126)", buffer = T)
plot(world, lwd = 0.5, add = T)

plot(usa1991_2010[[1]], main="United States (1991-2010)", buffer = T)
plot(world, lwd = 0.5, add = T)
plot(usa2011_2040ssp126[[1]], main="United States (2011-2040 ssp126)", buffer = T)
plot(world, lwd = 0.5, add = T)

plot(fra1991_2010[[1]], main="France (1991-2010)", buffer = T)
plot(world, lwd = 0.5, add = T)
plot(fra2011_2040ssp126[[1]], main="France (2011-2040 ssp126)", buffer = T)
plot(world, lwd = 0.5, add = T)
par(mfrow = c(1,1))
title('Hazelnut sourcing countries', outer = T, line = 1)

c


## binary maps present and future
#present
binaryMap <- ifel(ensProj >= cutoff, 1, NA)
levels(binaryMap) <- data.frame(id = c(1, 0), value = 'presence')
plot(binaryMap[[1]], col = 'gold3', main = 'Predicted presence/absence')
plot(world, add = T, lwd = 0.4)
#2011-2040
ssp126_2011_2040 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2011_2040_ssp126Ens.tif'))
binaryMapFut2011_ssp126 <- ifel(ssp126_2011_2040 >= cutoff, 1, NA)
levels(binaryMapFut2011_ssp126) <- data.frame(id = c(1, 0), value = 'presence')
plot(binaryMapFut2011_ssp126[[1]], col = 'orange4', main = 'Predicted presence/absence')
plot(world, add = T, lwd = 0.4)
#2041-2070
ssp126_2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2041_2070_ssp126Ens.tif'))
binaryMapFut2041_ssp126 <- ifel(ssp126_2041_2070 >= cutoff, 1, NA)
levels(binaryMapFut2041_ssp126) <- data.frame(id = c(1, 0), value = 'presence')
plot(binaryMapFut2041_ssp126[[1]], col = 'red4', main = 'Predicted presence/absence')
plot(world, add = T, lwd = 0.4)

## plot toghether
par(mfrow = c(1,1), oma = c(1,1,1,1))
plot(binaryMapFut2041_ssp126[[1]], col = 'red4', main = 'Predicted presence', legend = F, buffer = T)
plot(binaryMapFut2011_ssp126[[1]], col = 'orange4', legend = F, add = T)
plot(binaryMap[[1]], col = 'gold3', legend = F, add = T)
plot(world, add = T, lwd = 0.4)
mtext(paste0('Cutoff threshold of ', round(cutoff, 3)), side = 1, line = -2)
legend('bottomleft', y.intersp = 1, inset = c(0.02, 0.07), legend = c('1991-2010', '2011-2040 ssp126', '2041-2070 ssp126'), pt.bg = c('gold3', 'orange4', 'red4'), 
       bty = 'n', xpd = F, cex = 1, pch = 22, pt.cex = 2)







################# Interactive plots_________________________________
library("terra")
## load data
folder <- 'world' # name of the mother folder

modData <- read.csv(paste0('output/', folder, '/modellingData/modData.csv')) # train dataset
modDataSPAT <- vect(modData, geom = c('x', 'y'), crs = 'EPSG:4326') # spatialize it
nativeOccs <- read.csv(paste0('output/', folder, '/modellingData/nativeOccs.csv'))
nativeOccsSPAT <- terra::vect(nativeOccs, geom = c('decimalLongitude', 'decimalLatitude'), crs = 'epsg:4326')
# hazelnut crop mask
hazel <- rast('data/CROPGRIDSv1.08_hazelnut.nc')
hazelPol <- as.polygons(ifel(hazel[[1]] > 0, 1, NA))
# current climate
ensProj <- rast(paste0('output/', folder, '/projections/rast/ensemble_proj.tif'))
ensVar <- rast(paste0('output/', folder, '/projections/rast/ensemble_variance.tif'))
ext() <- ext(ensProj)
# future climate
#2011-2040
ssp126Ens_2011_2040 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2011_2040_ssp126Ens.tif'))
ext(ssp126Ens_2011_2040) <- ext(ensProj)
ssp370Ens_2011_2040 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2011_2040_ssp370Ens.tif'))
ext(ssp370Ens_2011_2040) <- ext(ensProj)
ssp585Ens_2011_2040 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2011_2040_ssp585Ens.tif'))
ext(ssp585Ens_2011_2040) <- ext(ensProj)
#2041-2070
ssp126Ens_2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2041_2070_ssp126Ens.tif'))
ext(ssp126Ens_2041_2070) <- ext(ensProj)
ssp370Ens_2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2041_2070_ssp370Ens.tif'))
ext(ssp370Ens_2041_2070) <- ext(ensProj)
ssp585Ens_2041_2070 <- rast(paste0('output/', folder, '/projections/Future/ensembles/modsEnsembles/2041_2070_ssp585Ens.tif'))
ext(ssp585Ens_2041_2070) <- ext(ensProj)

## interactive plots
# ensProj, ssp126Ens_2011_2040[[1]], ssp370Ens_2011_2040[[1]], ssp585Ens_2011_2040[[1]],
#   ssp126Ens_2041_2070[[1]], ssp370Ens_2041_2070[[1]], ssp585Ens_2041_2070[[1]])
rasters <- ensProj
terra::plet(hazelPol)




