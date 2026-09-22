############################################
############ hazelnut locations ############
############################################

library(terra)

######### Load data ----------
folder <- 'baseline_models_yearly_runs' # name of the mother folder

world <- geodata::world(path = 'data/geodata/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica
ita <- world[world$GID_0 == 'ITA']
srb <- world[world$GID_0 == 'SRB']
tur <- world[world$GID_0 == 'TUR']
geo <- world[world$GID_0 == 'GEO']

ens <- rast(file.path('output', folder, 'ensemble/ens.tif'))

hazelnut <- rast('data/otherLayers/CROPGRIDSv1.08_hazelnut.nc')
######### 


######### Hazelnut regions ----------
## hazelnut regions
hazelIta <- as.polygons(ifel(crop(hazelnut[[1]], ita) > 5, 1, NA), aggregate = TRUE)
plot(hazelIta)
plot(world, add = T)
hazelSrb <- as.polygons(ifel(crop(hazelnut[[1]], srb) > 5, 1, NA), aggregate = TRUE)
plot(hazelSrb)
plot(world, add = T)
hazelTurGeo <- as.polygons(ifel(crop(hazelnut[[1]], turgeo) > 5, 1, NA), aggregate = TRUE)
plot(hazelTurGeo)
plot(world, add = T)

## suitability in hazelnut regions
ensIta <- crop(ens, ita, mask = T)
ensSrb <- crop(ens, srb, mask = T)
ensTur <- crop(ens, tur, mask = T)
ensGeo <- crop(ens, geo, mask = T)

defaultPlotOptions <- par(no.readonly = T) # store original plotting options
par(mfrow = c(2,2))
plot(ensGeo[[1]], main = 'georgia', font.main = 1, plg = list(size = c(1, 0.2)))
plot(world, add = T)
plot(ensIta[[1]], main = 'italy', font.main = 1, plg = list(size = c(1, 0.2)))
plot(world, add = T)
plot(ensSrb[[1]], main = 'serbia', font.main = 1, plg = list(size = c(1, 0.2)))
plot(world, add = T)
plot(ensTur[[1]], main = 'turkey', font.main = 1, plg = list(size = c(1, 0.2)))
plot(world, add = T)
par(defaultPlotOptions) # restore plotting window to original








