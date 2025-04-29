### Code for Thompson et al 2025
### "Mixed effect of protected areas on avian food webs"


## Run sensitivity analysis for testing the effect of  % of overlap threshold 
## between PA polygons and grid cells

# We changed overlap rules for selecting grid cells for contrast analysis
# 90% overlap = at least 90% of grid cells must overlap with a PA for them to be considered protected
# 50% overlap = at least 50% of grid cells must overlap with a PA for them to be considered protected
# in both cases non-protected grid cells must have 0% overlap



library(stringr)
library(ggplot2)
library(plyr)
library(rnaturalearth)
library(sf)
library(raster)
library(foreign)
library(data.table)
# devtools::install_github("ropensci/rnaturalearthhires")

setwd("E:/TheseSwansea/WDPA")

europeRaster <- raster::raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")
cells_info <- read.dbf('data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img.vat.dbf')
cells_info = cbind(cells_info, coordinates(europeRaster)[which(!is.na(raster::values(europeRaster))),])
europe = ne_countries(returnclass = "sf", scale = 'large')
europe = st_transform(europe, crs(europeRaster))
europe = st_crop(europe, extent(europeRaster))



###Load datasets -----

## dataset with food web metrics per grid cell
res = readRDS("outputs/1.Metrics/NetworkMetrics.Rarefied.nb15.effort50.Over5yo.AfterDesig.NewEffort.rds")
res = merge(res, readRDS("outputs/1.Metrics/NetworkMetrics.Rarefied.nb15.effort50.Over5yo.AfterDesig.NewEffort.BodyMass.rds"))


## MAP FO RAREFIED SURVEY EFFORT
ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = effort.min), shape = 15, size = 0.5) + scale_color_gradientn(colours = rainbow(5))


## from '3.EnvironmentalVariables.R'
## set of environmental variables formatted to the grid cell level
## to combine with FW metrics
load("outputs/ContextVariables.Over5yo.RData")
## script with functions to combine environmental variables and 
## run GAMMs
source('Script/FunctionsAnalysis.R')

## Dataframe of WDPA shp turned into grid cell form and subset for europe TETRA EU extent
## it contains all PA info on protected grid cell 
PA.grid.nodups = readRDS("outputs/WDPA.Grid/distance/Grid.NoDups.PA_10km_Europe.POLYGONS.PointDistanceToEdge.Over5yo.RDS")

## add info about protection level and type of protection (directed at birds or not)
# Turning IUCN protection categories into a quantitative measure of protetcion level
levels = data.frame(IUCN_CAT = levels(factor(PA.grid.nodups$IUCN_CAT)), score = c(1,2,3,4,5,10,10,10,6,7))
PA.grid.nodups = merge(PA.grid.nodups, levels, by = 'IUCN_CAT')
PA.grid.nodups$Directed_at_birds = 
  ifelse(PA.grid.nodups$DESIG_ENG %in% c("Ramsar Site, Wetland of International Importance", "Special Protection Area (Birds Directive)"),
         1, 0)
PA.inf = unique(PA.grid.nodups[,c("ID","area")])
PA.grid.nodups = PA.grid.nodups[-which(!PA.grid.nodups$Value %in% res$Value),] # remove cells that weren't surveyed


# merged.PAs = table with all grid cells together (non protected and protected) and with 
# "NeighbouringPAID" whcih associates non protected grid cells to its PA network
merged.PAs = readRDS("outputs/WDPA.Grid/GridCells.NeighbouringPAID.new.Over5yo.RDS")
cells.to.seperate = c('555534657.1','555534704.1','555592567.1')
merged.PAs$NeighbouringPAID.new[which(merged.PAs$NeighbouringPAID %in% cells.to.seperate)] = '5.1'


## ADD CONTEXT VARIABLES
res = add_context_var(res = res)

setDF(res)
res = res[,-which(str_detect(names(res), 'NA'))]
res = subset(res, STATUS_YR != 0 | is.na(STATUS_YR))


## original way of selecting protected grid cells -- cells whose middle overlap 
## a PA polygon ("PA == 1")
## Number of protected and non protected grid cells under these rules:
ddply(res, c("Bioregion","NeighbouringPAID.new"), summarise, NeighbouringPAID.new = unique(NeighbouringPAID.new), 
      nb.inside = length(which(PA == 1)), nb.outside = length(which(PA == 0)), nb = unique(nb.inside + nb.outside))



### Extract overlap between PA polygons and grid cells -------------


# ## get cell IDs for grid cells that overlap PA boundary (see l 112 if already saved)
# PA.grid.nodups.shp = sf::st_read("outputs/WDPA.Shp/NewID.MergedPAs.1km.20012023.Over5yo.shp")
# ## we need to union the shapefile to remove any polygon overlap
# PA.grid.nodups.shp.combined = st_combine(unique(PA.grid.nodups.shp[,c('ID','ID_new','WDPAID')]))
# PA.grid.nodups.shp.combined = st_cast(PA.grid.nodups.shp.combined, 'POLYGON') 

## file has one row per grid cell, but we don't want that resolution here, we only want to keep one 
## polygon per unique original WDPA polygons
## union the polygon to eliminate overlap:

## if need to re-run union: (was run on laptop 21/03/2024)

# PA.grid.nodups.shp.union = st_union(PA.grid.nodups.shp.combined)
# st_write(PA.grid.nodups.shp.union, 'outputs/WDPA.Shp/Union.MergedPAs.1km.shp')

## or load if already save:
PA.grid.nodups.shp.union = st_read('outputs/WDPA.Shp/Union.MergedPAs.1km.shp')

## europeRaster needs to be a raster object 
europeRaster <- raster::raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")

# we use exact_extract function to get % coverage of PA polygon with each grid cell in europeRaster
# if the % coverage per grid cell is more than 1 if there are overlapping polygons. 
overlap = exactextractr::exact_extract(europeRaster, PA.grid.nodups.shp.union, include_cell = T)
overlap = overlap[[1]]
overlap = subset(overlap, !is.na(value))

summary(overlap) # coverage fraction should be max 1

## add coverage fraction of food web metrics
res = merge(res, overlap[,c('value','coverage_fraction')], by.x = 'Value', by.y = 'value', all.x = T)

summary(res$coverage_fraction) # the NAs are non-protected grid cells

dim(res[res$coverage_fraction >= 0.9,])

## check on a map
## 
## library(sf)

# ## checking that everything went as expected
# PA.grid.nodups.shp.eu = st_transform(PA.grid.nodups.shp, crs = crs(europeRaster))
# PA.grid.nodups.shp.eu.sub = merge(PA.grid.nodups.shp.eu, res[,c('Value','ID.new.Bioregion','coverage_fraction')], by = 'Value')
# 
# europeRaster$Coverage_fraction = res$coverage_fraction[match(values(europeRaster),res$Value)]
# raster_to_pol_sub = st_as_sf(raster::rasterToPolygons(europeRaster$Coverage_fraction))

# map = ggplot() + 
#   geom_sf(data = PA.grid.nodups.shp.eu.sub['ID_new' == 7,], fill = NA) + 
#   geom_sf(data = raster_to_pol_sub, aes(fill = Coverage_fraction), alpha = 0.5) +
#   xlim(4705960, 5005960) + ylim(4871218,5171218)
# 
# ggsave(plot = map, paste0('Plots/Figures/test_map.png'), dpi = 200,scale = 1.5)
 
## checking that everything looks fine on qgis
# st_write(raster_to_pol_sub,'outputs/WDPA.Shp/raster_coverage_fraction.shp' )


## tested 50 and 90% thresholds
## choose desired threshold
thresh = 0.9 # threshold for considering that a gridcell is protected - if it has > thresh % coverage by a PA
M.merged = ddply(res, "ID.new.Bioregion", summarise, NeighbouringPAID.new = unique(NeighbouringPAID.new), 
      nb.inside = length(which(coverage_fraction >= thresh & PA == 1)), 
      # there were a few cases where the PA didn't cover the middle of the grid cell (PA == 0) but 
      # still had > 50% coverage by PA, these grid cells are removed here (by specifying PA == 1)
      nb.outside = length(which(is.na(coverage_fraction))), 
      nb = unique(nb.inside + nb.outside), protection_coverage = thresh)

M.merged = subset(M.merged, nb.inside >= 15 & nb.outside >= 15) ## we only keep 33 if we choose 50% coverage fraction


## original way of selecting protected grid cells -- cells whose middle overlap 
## a PA polygon ("PA == 1")
## Number of protected and non protected grid cells under these rules:
M.merged_og = ddply(res, "ID.new.Bioregion", summarise, NeighbouringPAID.new = unique(NeighbouringPAID.new), 
                 nb.inside = length(which(PA == 1)), nb.outside = length(which(PA == 0)), nb = unique(nb.inside + nb.outside))
M.merged_og = subset(M.merged_og, nb.inside >= 15 & nb.outside >= 15) ## we only keep 33 if we choose 50% coverage fraction

# graph of evolution of data set size depending on chosen threshold
dt_size_conservative = data.frame()
dt_size = data.frame()
count = 0
for (prop in seq(0,1,0.05)){
  
  count = count + 1
  M.merged.temp_conservative = ddply(res, "ID.new.Bioregion", summarise, NeighbouringPAID.new = unique(NeighbouringPAID.new), 
                                     nb.inside = length(which(coverage_fraction >= prop)), nb.outside = length(which(is.na(coverage_fraction))), nb = unique(nb.inside + nb.outside))
  
  M.merged.temp_conservative = subset(M.merged.temp, nb.inside >= 15 & nb.outside >= 15) 
  
  M.merged.temp = ddply(res, "ID.new.Bioregion", summarise, NeighbouringPAID.new = unique(NeighbouringPAID.new), 
                        nb.inside = length(which(coverage_fraction > prop)), nb.outside = length(which(coverage_fraction <= prop | is.na(coverage_fraction))), nb = unique(nb.inside + nb.outside))
  
  M.merged.temp = subset(M.merged.temp, nb.inside >= 15 & nb.outside >= 15) 
  
  dt_size[count,c("prop","nb_sites")] = c(prop, nrow(M.merged.temp))
  dt_size_conservative[count,c("prop","nb_sites")] = c(prop, nrow(M.merged.temp_conservative))
  
}

# using the coverage as a threshold means that we have a hard line between what is protected and what is not
# so that we never reach that 46 site number that we do when using the 'middle overlap' rule 

par(mfrow = c(1,1))
plot(dt_size)
points(dt_size_conservative, col = 'blue')


### subset grid cells according to new coverage threshold
res = subset(res, ID.new.Bioregion %in% M.merged$ID.new.Bioregion & ((coverage_fraction > thresh & PA == 1) | is.na(coverage_fraction)))


# checks
print(paste("Should be TRUE: ", nrow(res) == sum(M.merged$nb))) 
# chekc that the nb of grid cells (= rows) match what we calculated previously

ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = m.gHP), size = 0.8) + scale_colour_gradientn(colors = rainbow(5))
#ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = m.THPI))

summary(res)
summary(res$coverage_fraction) # should be > thresh (protected) or NAs (unprotected grid cells)
all(res[which(!is.na(res['coverage_fraction'])),'PA'] == 1) # should be TRUE

length(unique(res$ID.new.Bioregion))
anyNA(res$LC1) # should be False








### RUN GAMMS for sensitivity analysis ------
### to compute difference in FW metrics between 
### protected and unprotected grid cells


library(mgcv)
# indices for which we want to compute the contrast in network metrics between
# protected and unprotected grid cells
Indices.list = c("S","C","L.S", "L", "modularity.2", "mfcl", "omnivory.1", "Gen", 
                 "Vul", "TL.mean", "S_top", "S_basal","S_intermediate","Ftop",
                 "Fbasal","FInt","SDGen","SDVul","mass_basal","mass_top","mass_int")

source('Script/FunctionsAnalysis_MixedModels_RandomSlopeControl.R')

set.seed(1997)
res$LC2 = as.factor(res$LC2)
res$LC = as.factor(res$LC)
setDF(res)

res$mass_top[which(is.na(res$mass_top))] = 0
res$mass_int[which(is.na(res$mass_int))] = 0


### change name of folder depending on the threshold chosen: (/90coverage or /50coverage)
outputs = run_GAMMs(res_resampled = res, Indices.list = Indices.list)
ranef = outputs[[1]]
fixef = outputs[[2]]

write.csv(ranef, 'outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/90coverage/GAMM-randomEffect-FullDataset-RandomSlope-FixedSmooths.csv')
write.csv(fixef, 'outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/90coverage/GAMM-fixedEffect-FullDataset-RandomSlope-FixedSmooths.csv')




