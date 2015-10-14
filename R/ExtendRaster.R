################################################################################
# 
# ExtendRaster.R
# Version 1.0
# 10/03/2015
#
# Takes the raster layer at grain size i and the raster layer at i^max and
# extends raster(i) to same extent as raster(i^max). All new cells added are
# given values of 0
#
# max_raster - the atlas data at the maximum grain size
# scaled_raster - the atlas data at grain size i
#
################################################################################

ExtendRaster <- function(max_raster, scaled_raster, scale) { 
  max_raster_scaled <- raster::disaggregate(max_raster, 2 ^ scale, fun = max)
  max_raster_scaled[max_raster_scaled == 1] <- 0
  scaled_raster <- raster::extend(scaled_raster, 
                                  raster::extent(max_raster_scaled))
  scaled_raster[is.na(scaled_raster)] <- 0
  scaled_raster <- raster::overlay(scaled_raster, max_raster_scaled, fun=sum)
  return(scaled_raster)
}