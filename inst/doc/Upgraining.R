## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(downscale)

## ----OriginalUpgrain, echo = FALSE, out.width = "100%", fig.align = "left", fig.cap = "Upgrained presence (red cells) and absence (white cells) maps for a UK species without standardising extent to the largest grain size. Unsampled cells are dark grey. As we upgrain the atlas data to larger grain sizes the total extent also increases."----
knitr::include_graphics("figures/Original_upgrain.png")

## ----AllInterior, echo = FALSE, out.width = "100%", fig.align = "left", fig.cap = "Upgrained presence (red cells) and absence (white cells) maps for a UK species after standardising extent to the largest grain size. Unsampled cells are dark grey. The extent of the atlas data is extended to that of the largest grain size by assigning absences to unsampled cells."----
knitr::include_graphics("figures/All_interior.png")

## ----InteriorOnly, echo = FALSE, out.width = "100%", fig.align = "left", fig.cap = "Upgrained presence (red cells) and absence (white cells) maps for a UK species after standardising extent to those cells at the largest grain size that solely contain sampled atlas data. Sampled cells outside the selected cells are assigned as No Data (dark grey)."----
knitr::include_graphics("figures/Interior_only.png")

## ----fig.show = "hide"--------------------------------------------------------
# The data may be a raster layer of presence (1) and absence (0) data or a
# data frame of cell center coordinates and presence-absence data in which case
# it must have these column names: "x", "y", "presence"
data.file <- system.file("extdata", "atlas_data.txt", package = "downscale")
atlas.data <- read.table(data.file, header = TRUE)

# run upgrain.threshold for three larger grain sizes
thresh <- upgrain.threshold(atlas.data = atlas.data,
                            cell.width = 10,
                            scales     = 3,
                            thresholds = seq(0, 1, 0.01))

## ----ThresholdPlots, echo = FALSE, out.width = "100%", fig.align = "left", fig.cap = 'Diagnostic plots produced by `upgrain.threshold` used to explore the trade-off  between assigning large areas of unsampled areas as absence, and discarding sampled areas and known presences. Two possible thresholds in the quantity of unsampled area allowed within cells at the largest grain size are identified: the "All_Occurrences" threshold (blue line) and the "Gain_Equals_Loss" threshold (red line).'----
knitr::include_graphics("figures/Threshold_plots.png")

## -----------------------------------------------------------------------------
thresh$Thresholds

## -----------------------------------------------------------------------------
head(thresh$Data)

## ----ThresholdMaps, echo = FALSE, out.width = "100%", fig.align = "left", fig.cap = "Maps of the atlas data (red = presence; light grey = absence; unsampled = dark grey) overlain with polygons showing the standardised extent after applying each of four possible thresholds."----
knitr::include_graphics("figures/Threshold_maps.png")

## ----eval = FALSE-------------------------------------------------------------
#  # run upgrain.threshold for two larger grain sizes
#  thresh <- upgrain.threshold(atlas.data = atlas.data,
#                              cell.width = 10,
#                              scales     = 2,
#                              thresholds = seq(0, 1, 0.01))

## ----ThresholdMaps2, echo = FALSE, out.width = "100%", fig.align = "left", fig.cap = "Maps for the four possible thresholds after upgraining across two scales."----
knitr::include_graphics("figures/Threshold_maps_2_scales.png")

## ----eval = FALSE-------------------------------------------------------------
#  # run upgrain.threshold for four larger grain sizes
#  thresh <- upgrain.threshold(atlas.data = atlas.data,
#                              cell.width = 10,
#                              scales     = 4,
#                              thresholds = seq(0, 1, 0.01))

## ----ThresholdMaps4, echo = FALSE, out.width = "100%", fig.align = "left", fig.cap = "Maps for the four possible thresholds after upgraining across four scales."----
knitr::include_graphics("figures/Threshold_maps_4_scales.png")

