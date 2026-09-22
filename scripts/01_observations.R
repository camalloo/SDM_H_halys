######################################################
############ Process the raw observations ############
######################################################

library(tidyverse)
library(terra)

folder <- 'baseline_models_yearly_runs' # name of the mother folder

######### Prepare observation data ----------
## world map
world <- geodata::world(path = 'data/geodata', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antarctica

## raw data
# dirs <- list.dirs('data/observations', recursive = F)

GBIF <- 'data/observations/GBIF/occurrence.csv'
dataRaw <- read.csv(GBIF)

## clean the dataframe
dataClean <- dataRaw |> 
  select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, countryCode,
         decimalLatitude, decimalLongitude, verbatimEventDate, year, month, day, startDayOfYear, endDayOfYear) |> 
  mutate(eventDate = as.Date(eventDate, '%d/%m/%Y')) |> 
  ## make the date from day, month, year when available
  mutate(eventDate = coalesce(eventDate, as.Date(paste(day, month, year, sep = '/'), '%d/%m/%Y'))) |> 
  ## make the date from doy and year when available 
  mutate(eventDate = coalesce(eventDate, as.Date(startDayOfYear-1, origin = paste0(year, '/01/01')))) |>
  ## convert verbatim dates
  mutate(eventDate = coalesce(eventDate, as.Date(verbatimEventDate, format = '%d/%m/%Y %H:%M'))) |> 
  ## drop no coords or date
  drop_na(decimalLatitude, decimalLongitude, eventDate)

## remove duplicated records and check for other error in coordinates --> FINAL DF
dataClean <- dataClean |> 
  CoordinateCleaner::cc_dupl(lon = 'decimalLongitude', lat = 'decimalLatitude', addition = 'eventDate', value = 'clean') |> 
  CoordinateCleaner::clean_coordinates(lon = 'decimalLongitude', lat = 'decimalLatitude', countries = 'countryCode', 
                                       tests = c('centroids', 'equal', 'gbif', 'outliers', 'seas', 'zeros'),
                                       value = 'clean') |> 
  select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, year, month, day, countryCode,
         lat = decimalLatitude, long = decimalLongitude)
## additional cleaning with FuzzySim
dataClean <- fuzzySim::cleanCoords(dataClean, coord.cols = c('long', 'lat'), uncert.col = NULL)

## plot all
plot(world)
points(dataClean$long, dataClean$lat, pch = 20, cex = 0.3, col = 'red')

## occurrences from the native areas
native <- dataClean |> 
  filter(countryCode %in% c('CN', 'TW', 'JP', 'KP', 'KR', 'HK', 'MO'))
## plot
plot(world, main = expression(paste(italic('H. halys'), ' ditribution')))
points(native$long, native$lat, col = 'blue', pch = 16, cex = 0.3)

## occurrences from invaded areas
invasive <- dataClean |> 
  filter(!(countryCode %in% c('CN', 'TW', 'JP', 'KP', 'KR', 'HK', 'MO')))
## add to plot
points(invasive$long, invasive$lat, col = 'red', pch = 16, cex = 0.3)

## legend of the distribution plot
legend('bottomleft', horiz = F, legend = c('Native Range', 'Invaded Range'), 
       y.intersp = 0.2, xpd = T, inset = c(0, 0.15),
       col = c('blue', 'red'), pch = 16, pt.cex = 1, bty = 'n')

## save cleaned occs to .csv
write.csv(native |> select(eventDate, lat, long), file.path('output', folder, 'modellingData/native_occs.csv'), row.names = F)
write.csv(invasive |> select(eventDate, lat, long), file.path('output', folder, 'modellingData/invasive_occs.csv'), row.names = F)
write.csv(dataClean |> select(eventDate, lat, long), file.path('output', folder, 'modellingData/clean_occs.csv'), row.names = F)

## how many obs per year?
obsXyear <- data.frame(year = unique(dataClean$year))
obsXyear <- obsXyear |> 
  left_join(count(invasive, year, name = 'invasive'), by = 'year') |> 
  left_join(count(native, year, name = 'native'), by = 'year') |> 
  arrange(year)
## plot
ggplot(data = obsXyear, aes(x = year)) +
  geom_line(aes(y = invasive), linetype = 'dashed', na.rm = T) +
  geom_point(aes(y = invasive)) +
  geom_line(aes(y = native), linetype = 'dashed', col = 'red', na.rm = T) +
  geom_point(aes(y = native), col = 'red')


