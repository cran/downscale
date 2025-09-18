## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(downscale)

## ----OriginalUpgrain, echo = FALSE, out.width = "100%", fig.align = "left", fig.cap = '"Upraining" of atlas data at a 10km cell width to 20km, 40km and 80 km. As the grain size is increased occupancy also increases.'----
knitr::include_graphics("figures/Original_upgrain.png")

## ----eval = FALSE-------------------------------------------------------------
# vignette("Upgraining", package = "downscale")

## ----Saturation, echo = FALSE, out.width = "100%", fig.align = "left", fig.cap = "Occupancy-area relationships (OAR) for two species showing a) the scale of saturation (the grain size at which all cells are occupied) and b) the scale of endemism (the scale at which only one cell is occupied)."----
knitr::include_graphics("figures/Saturation.png")

## ----Flow, echo = FALSE, out.width = "70%", fig.align = "left", fig.cap = "Structure of the `downscale` package showing the flow between the seven functions (yellow) and the three output object classes (orange). Black arrows represent the input of raw data of four types: a raster layer or a spatial points layer of presence-absence data, a data frame of cell coordinates and presence-absence for each cell, and a data frame of occupancies at coarse-grain sizes. Yellow arrows are where the output of one function may be used as the input for the next function."----
knitr::include_graphics("figures/Flow.png")

## ----eval = FALSE-------------------------------------------------------------
# install.packages("downscale")

## ----eval = FALSE-------------------------------------------------------------
# library("downscale")

## -----------------------------------------------------------------------------
library("sf")
library("terra")

## -----------------------------------------------------------------------------
occupancy <- data.frame(Cell.area = c(100, 400, 1600, 6400),
                        Occupancy = c(0.23, 0.56, 0.87, 1))

## -----------------------------------------------------------------------------
## fit logistic model to observed data using downscale
logisMod <- downscale(occupancies = occupancy,
                      model       = "Logis",
                      extent      = 384000)

## -----------------------------------------------------------------------------
logisMod

## -----------------------------------------------------------------------------
## new grain sizes to predict
areasPred <- c(1, 2, 5, 25, 100, 400, 1600, 6400)

## predict for the new grain sizes using the downscale object
logisPred <- predict(logisMod,
                     new.areas = areasPred,
                     plot = FALSE)
                      
## this creates an object of class 'predict.downscale'
## occupancy is given as a proportion (Occupancy) and area of occupancy (AOO)
logisPred$predicted

## ----fig.width = 6, fig.height = 4--------------------------------------------
plot(logisPred)

## -----------------------------------------------------------------------------
## if it is not already loaded, load in the package
library(downscale)

## -----------------------------------------------------------------------------
dataFile <- system.file("extdata", "atlas_data.txt", package = "downscale")
atlasData <- read.table(dataFile, header = TRUE)
head(atlasData)

## ----eval = FALSE-------------------------------------------------------------
# vignette("Upgraining", package = "downscale")

## ----fig.show = "hide"--------------------------------------------------------
## explore thresholds using upgrain.threshold
thresh <- upgrain.threshold(atlas.data = atlasData,
                            cell.width = 10,
                            scales     = 3)

## ----ThresholdPlots, echo = FALSE, out.width = "100%", fig.align = "left"-----
knitr::include_graphics("figures/Threshold_plots.png")

## ----ThresholdMaps, echo = FALSE, out.width = "100%", fig.align = "left"------
knitr::include_graphics("figures/Threshold_maps.png")

## -----------------------------------------------------------------------------
thresh$Thresholds

## ----fig.height = 5, fig.width = 7--------------------------------------------
## upgrain data (using All Occurrences threshold)
occupancy <- upgrain(atlas.data = atlasData,
                     cell.width = 10,
                     scales     = 3,
                     method     = "All_Occurrences",
                     plot       = TRUE)

## -----------------------------------------------------------------------------
## Improved Negative Binomial model
(inb <- downscale(occupancies = occupancy,
                  model = "INB"))

## -----------------------------------------------------------------------------
### Manually specifying the starting parameters
paramsNew <- list("C" = 0.1, "gamma" = 0.00001, "b" = 0.1)
inbNew <- downscale(occupancies = occupancy,
                     model = "INB",
                     starting_params = paramsNew)

## ----fig.width = 5, fig.height = 4--------------------------------------------
## plot the predictions of two FNB models using predict.downscale
inbPred <- predict(inb,
                   new.areas = c(1, 2, 5, 25, 100, 400, 1600, 6400),
                   plot = TRUE)
inbPredNew <- predict(inbNew,
                      new.areas = c(1, 2, 5, 25, 100, 400, 1600, 6400),
                      plot = TRUE)

## ----fig.width = 5, fig.height = 4--------------------------------------------
## Thomas model
thomas <- downscale(occupancies = occupancy,
                    model       = "Thomas",
                    tolerance   = 1e-3)
                    
## the tolerance can also be set for the predict function
thomas.pred <- predict(thomas,
                       new.areas = c(1, 2, 5, 25, 100, 400, 1600, 6400),
                       tolerance = 1e-6)

## ----fig.width = 5, fig.height = 4--------------------------------------------
plot(thomas.pred,
     col.pred = "green",  # change the colour of the prediction
     pch      = 16,       # change point character
     lwd.obs  = 3)        # change line width of the observed data

## ----fig.width = 5, fig.height = 4--------------------------------------------
## Hui model using a data frame as input
hui <- hui.downscale(atlas.data = atlasData,
                     cell.width = 10,
                     extent     = 228900,
                     new.areas  = c(1, 2, 5, 15, 50))

## the output is a normal 'predict.downscale' object	
plot(hui)

## ----fig.width = 5, fig.height = 4--------------------------------------------
huiStand <- hui.downscale(atlas.data = occupancy,
                          cell.width = 10,
                          new.areas  = c(1, 2, 5, 15, 50),
                          plot       = TRUE)

## -----------------------------------------------------------------------------
hui$predicted

## -----------------------------------------------------------------------------
huiStand$predicted

## -----------------------------------------------------------------------------
## hypothetical occupancy data
occupancy <- data.frame(Cell.area = c(100, 400, 1600, 6400),
                        Occupancy = c(0.23, 0.56, 0.87, 1))
                        
## grain sizes (cell areas) to predict
areasPred <- c(1, 2, 5, 25, 100, 400, 1600, 6400)

## ----fig.width = 7, fig.height = 6--------------------------------------------
ensemble <- ensemble.downscale(occupancies = occupancy,
                               new.areas   = areasPred,
                               extent      = 320000,
                               models      = c("PL",
                                               "Logis",
                                               "NB",
                                               "GNB",
                                               "INB"),
                               plot        = TRUE)

## -----------------------------------------------------------------------------
ensemble$Occupancy

## -----------------------------------------------------------------------------
ensemble$AOO

## ----fig.width = 7.5, fig.height = 6------------------------------------------
dataFile <- system.file("extdata", "atlas_data.txt", package = "downscale")
atlasData <- read.table(dataFile, header = TRUE)

## upgrain data (using "All Occurrences" threshold)
occupancy <- upgrain(atlas.data = atlasData,
                     cell.width = 10,
                     scales     = 3,
                     method     = "All_Occurrences",
                     plot       = FALSE)

## ensemble modelling
ensemble <- ensemble.downscale(occupancies = occupancy,
                               new.areas   = areasPred,
                               cell.width  = 10,
                               models      = c("Nachman",
                                               "PL",
                                               "Logis",
                                               "GNB",
                                               "FNB",
                                               "Hui"),
                               plot         = TRUE)

## ----fig.width = 7.5, fig.height = 6------------------------------------------
ensemble <- ensemble.downscale(occupancies   = occupancy,
                               new.areas     = areasPred,
                               cell.width    = 10,
                               models        = "all",
                               tolerance_mod = 1e-3,
                               plot          = TRUE)

## ----fig.width = 7.5, fig.height = 6------------------------------------------
## Specifying starting parameters for Nachman and GNB models
newParams <- list(Nachman = list("C" = 0.1, "z" = 0.01),
                  GNB = list("C" = 0.1, "z" = 1, "k" = 0.01))
newParams

ensemble <- ensemble.downscale(occupancies     = occupancy,
                               new.areas       = areasPred,
                               cell.width      = 10,
                               models          = "all",
                               tolerance_mod   = 1e-3,
                               starting_params = newParams,
                               plot            = TRUE)

## -----------------------------------------------------------------------------
## load in the necessary libraries
library(rgbif)

## ----eval = FALSE-------------------------------------------------------------
# ### NOT RUN in this tutorial
# ### Run this code to download the data yourself (results will slightly differ)
# records <- occ_search(scientificName = "Polyommatus coridon",
#                       country        = "GB",
#                       limit          = 10000,
#                       hasCoordinate  = TRUE,
#                       fields         = "minimal")
# records <- records$data

## -----------------------------------------------------------------------------
recordsFile <- system.file("extdata", "Polyommatus_coridon_gbif_records.txt",
                           package = "downscale")
records <- read.table(recordsFile, header = TRUE)

## -----------------------------------------------------------------------------
recordsCoords <- st_as_sf(records,
                          coords = c("decimalLongitude", "decimalLatitude"),
                          crs = "OGC:CRS84")

## -----------------------------------------------------------------------------
## reproject the coordinates to British National Grid
BNG <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs"
recordsCoords <- st_transform(recordsCoords,
                              crs = BNG)

## ----fig.width = 6, fig.height = 5--------------------------------------------
plot(st_geometry(recordsCoords), axes = TRUE)

## ----fig.width = 5, fig.height = 4--------------------------------------------
## set grain size as 20 km
cellWidth <- 20

# extract extent of coordinates
coordsExtent <- ext(recordsCoords)

## create a blank raster to fit the coordinates (note the addition of half a 
## cell width on all sides)
gbif_raster <- rast(xmin = coordsExtent$xmin - (cellWidth / 2),
                    xmax = coordsExtent$xmax + (cellWidth / 2),
                    ymin = coordsExtent$ymin - (cellWidth / 2),
                    ymax = coordsExtent$ymax + (cellWidth / 2),
                    res = cellWidth)
                      
## assign cells with presence records as 1
gbif_raster <- rasterize(recordsCoords, gbif_raster, field = 1, fun = "min")

## convert cells with NA (no records) to 0
gbif_raster[is.na(gbif_raster)] <- 0

plot(gbif_raster, legend = FALSE)

## ----fig.height = 5, fig.width = 7--------------------------------------------
occupancy <- upgrain(atlas.data = gbif_raster,
                     scales     = 2,
                     method     = "All_Sampled")

## -----------------------------------------------------------------------------
## The extent of the original atlas data
occupancy$occupancy.orig[1, 2]

## The extent of the standardised atlas data
occupancy$extent.stand

## ----fig.height = 6, fig.width = 7.5------------------------------------------
ensemble <- ensemble.downscale(occupancies   = occupancy,
                               models        = c("all"),
                               new.areas     = c(1, 10, 100, 400, 1600, 6400),
                               tolerance_mod = 1e-3)

## ----fig.height = 6, fig.width = 7.5------------------------------------------
ensemble <- ensemble.downscale(occupancy,
                               models = "all",
                               new.areas = c(1, 10, 100, 400, 1600, 6400),
                               tolerance_mod = 1e-3,
                               starting_params = list(INB = list(C = 10, 
                                                                 gamma = 0.01, 
                                                                 b = 0.1)))

## -----------------------------------------------------------------------------
ensemble$AOO[, c("Cell.area", "Means")]

## ----fig.height = 4, fig.width = 5--------------------------------------------
### read in the shapefile
uk <- system.file("extdata", "UK.shp", package = "downscale")
uk <- st_read(uk)

## plot our GBIF records on top of the UK polygon
plot(st_geometry(uk), axes = TRUE)
plot(recordsCoords[1], add = TRUE, col = "red")

## ----fig.height = 5, fig.width = 4--------------------------------------------
## create a blank raster with the same extent as the UK polygon
gbif_raster <- rast(ext = ext(uk),
                    res = cellWidth,
                    crs = crs(uk))

## assign cells with presence records as 1
gbif_raster <- rasterize(recordsCoords, gbif_raster, field = 1)

## convert cells with NA (no records) to 0
gbif_raster[is.na(gbif_raster)] <- 0

## mask the raster to the UK polygon, so cells outside the polygon are NA
gbif_raster <- mask(gbif_raster, uk)

## plot the masked atlas raster and overlay with the UK polygon
plot(gbif_raster, legend = FALSE)
plot(st_geometry(uk), add = TRUE)

## ----fig.height = 5, fig.width = 7--------------------------------------------
occupancy <- upgrain(gbif_raster,
                     scales = 2,
                     method = "All_Sampled")

## ----fig.height = 6, fig.width = 7.5------------------------------------------
ensembleUK <- ensemble.downscale(occupancy,
                                 models = "all",
                                 new.areas = c(1, 10, 100, 400, 1600, 6400),
                                 tolerance_mod = 1e-3,
                                 starting_params = list(INB = list(C = 10,
                                                                   gamma = 0.01,
                                                                   b = 0.1),
                                                        Thomas = list(rho = 1e-6,
                                                                      mu = 1,
                                                                      sigma = 1)))

## -----------------------------------------------------------------------------
ensembleUK$AOO[, c("Cell.area", "Means")]

