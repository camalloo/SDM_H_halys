library("tidyverse")


################# Load and clean distribution data_________________________________
## raw data
dataPath <- "data/GBIF/raw/occurrence.csv"

dataRaw <- read.csv(dataPath)

## world map
world <- geodata::world(path = 'data/geodata/', resolution = 1)
# colnames(dataRaw)

## clean the dataframe
dataClean <- dataRaw |> 
  select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, countryCode,
         decimalLatitude, decimalLongitude, verbatimEventDate, year, month, day, startDayOfYear, endDayOfYear) |> 
  mutate(eventDate = as.Date(eventDate, "%d/%m/%Y")) |> 
  ## make the date from day, month, year when available
  mutate(eventDate = coalesce(eventDate, as.Date(paste(day, month, year, sep = "/"), "%d/%m/%Y"))) |> 
  ## make the date from doy and year when available 
  mutate(eventDate = coalesce(eventDate, as.Date(startDayOfYear-1, origin = paste0(year, "/01/01")))) |>
  ## convert verbatim dates
  mutate(eventDate = coalesce(eventDate, as.Date(verbatimEventDate, format = "%d/%m/%Y %H:%M"))) |> 
  ## drop no coords or date
  drop_na(decimalLatitude, decimalLongitude, eventDate) |> 
  ## remove native countries (China, Japan, Taiwan, Koreas, (Hong Kong and Macao included by me)) according to Lee et al., 2013
  filter(!(countryCode %in% c("CN", "TW", "JP", "KP", "KR", "HK", "MO")))

## remove duplicated records and check for other error in coordinates --> FINAL DF
dataClean <- dataClean |> 
  CoordinateCleaner::cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude", addition = "eventDate", value = "clean") |> 
  CoordinateCleaner::clean_coordinates(lon = "decimalLongitude", lat = "decimalLatitude", countries = "countryCode", 
                                       tests = c("centroids", "equal", "gbif", "outliers", "seas", "zeros"),
                                       value = "clean") |> 
  select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, countryCode,
         decimalLatitude, decimalLongitude)

terra::plot(world, main = 'GBIF observations (downloaded on 18/06/2025)')
terra::plot(terra::vect(dataClean, geom = c("decimalLongitude", "decimalLatitude"), crs = 'epsg:4326', keepgeom = T),
            col = 'red', cex = 0.3, add = T)
# points(nativeOccsSP, col = 'blue', cex = .5, pch = 4)
# legend('bottomleft', legend = c('Invasive observation', 'Native observation'), col = c('red', 'blue'),
#        pch = c(20, 4), xpd = F, inset = c(0.01, 0.12), bty = 'n')





################# First occurrences
# ## load first records extracted from literature
# firstOccLiterature <- readxl::read_xlsx("data/BMSB_First_Records.xlsx")
# firstOccLiterature <- firstOccLiterature |> 
#   subset(select = -c(state_region, district, municipality_city)) |> 
#   drop_na(decimalLatitude, decimalLongitude) |> 
#   mutate(occurrenceID = as.character(occurrenceID),
#          eventDate = as.Date(eventDate, format = "%d/%m/%Y"))
# ## merge it with gbif records
# all <- rbind(dataClean, firstOccLiterature)
# 
# ## get first occurrence per country
# firstOccurrences <- all |>
#   group_by(countryCode) |>
#   slice_min(order_by = eventDate, n = 1)
# 
# spatFstOcc <- terra::vect(firstOccurrences, geom = c("decimalLongitude", "decimalLatitude"), keepgeom = T) # spatialize the dataset
# plot(world)
# plot(spatFstOcc, add = T, col = year(firstOccurrences$eventDate))
# plot(firstOccurrences, col = year(firstOccurrences$eventDate), add = T)
# legend('bottom', legend = unique(year(firstOccurrences$eventDate)), fill = unique(year(firstOccurrences$eventDate)), ncol = 4, xpd = T)
