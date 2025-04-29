### Code for "Mixed effect of protected areas on avian food webs"


## Running GAMMs on food web metrics to obtain 
## difference in FW metrics between protected and 
## non-protected grid cells across PA networks



# Load libraries and maps ######################################################


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



# Load datasets and add context variables ######################################

## dataset of food web metrics per grid cell
res = readRDS("outputs/1.Metrics/NetworkMetrics.Rarefied.nb15.effort50.Over5yo.AfterDesig.NewEffort.rds")
## add body mass metrics to 'res' 
res = merge(res, readRDS("outputs/1.Metrics/NetworkMetrics.Rarefied.nb15.effort50.Over5yo.AfterDesig.NewEffort.BodyMass.rds"))

res$mass_top[which(is.na(res$mass_top))] = 0
res$mass_int[which(is.na(res$mass_int))] = 0


## from '3.EnvironmentalVariables.R'
## set of environmental variables formatted to the grid cell level
## to combine with FW metrics
load("outputs/ContextVariables.Over5yo.RData")
## script with functions to combine environmental variables and 
## run GAMMs
source('Script/FunctionsAnalysis_MixedModels_RandomSlopeControl.R')

## Dataframe of WDPA shp turned into grid cell form and subset for europe TETRA EU extent
## it contains all PA info on protected grid cell 
PA.grid.nodups = readRDS("outputs/WDPA.Grid/distance/Grid.NoDups.PA_10km_Europe.POLYGONS.PointDistanceToEdge.Over5yo.RDS")

# Turning IUCN protection categories into a quantitative measure of protetcion level
levels = data.frame(IUCN_CAT = levels(factor(PA.grid.nodups$IUCN_CAT)), score = c(1,2,3,4,5,10,10,10,6,7))
PA.grid.nodups = merge(PA.grid.nodups, levels, by = 'IUCN_CAT')
PA.grid.nodups$Directed_at_birds = 
  ifelse(PA.grid.nodups$DESIG_ENG %in% c("Ramsar Site, Wetland of International Importance", "Special Protection Area (Birds Directive)"),
         1, 0)
PA.inf = unique(PA.grid.nodups[,c("ID","area")])
PA.grid.nodups = PA.grid.nodups[-which(!PA.grid.nodups$Value %in% res$Value),] # remove cells that weren't surveyed



# CREATE A TABLE WITH ALL THE CARACTERISTICS OF EACH GRID PA ##################

# merged.PAs = table with all grid cells together (non protected and protected) and with 
# "NeighbouringPAID" whcih associates non protected grid cells to its PA network
merged.PAs = readRDS("outputs/WDPA.Grid/GridCells.NeighbouringPAID.new.Over5yo.RDS")
cells.to.seperate = c('555534657.1','555534704.1','555592567.1')
merged.PAs$NeighbouringPAID.new[which(merged.PAs$NeighbouringPAID %in% cells.to.seperate)] = '5.1'


## Add PAs context variables to food web metric data ###########################
res = add_context_var(res = res)

setDF(res)
res = res[,-which(str_detect(names(res), 'NA'))]
res = subset(res, STATUS_YR != 0 | is.na(STATUS_YR))
res = res[-which(res$S == 0),] # remove one community which had non interacting birds

# res = res[-which(is.na(res$temp)),]


removed = subset(res, ID.new.Bioregion == "2_Continental Bio-geographical Region")
res = subset(res, ID.new.Bioregion != "2_Continental Bio-geographical Region")


## MAP FO RAREFIED SURVEY EFFORT
ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = effort.min), shape = 15, size = 0.5) + scale_color_gradientn(colours = rainbow(5))


ddply(res, c("Bioregion","NeighbouringPAID.new"), summarise, NeighbouringPAID.new = unique(NeighbouringPAID.new), 
      nb.inside = length(which(PA == 1)), nb.outside = length(which(PA == 0)), nb = unique(nb.inside + nb.outside))

## get number of protected and non-protected grid cells 
## here a cell is considered protected if its middle is covered by a PA polygon
M.merged = ddply(res, "ID.new.Bioregion", summarise, NeighbouringPAID.new = unique(NeighbouringPAID.new), 
                 nb.inside = length(which(PA == 1)), nb.outside = length(which(PA == 0)), nb = unique(nb.inside + nb.outside))
M.merged = subset(M.merged, nb.inside >= 15 & nb.outside >= 15) 
# M.merged = M.merged[-which(M.merged$ID.new.Bioregion == "2_Continental Bio-geographical Region"),]
## this site is in the east and CLC doesn't cover that far east

res = subset(res, ID.new.Bioregion %in% M.merged$ID.new.Bioregion)



### Summary stats on food web metrics for manuscript #############################
summary(res)
anyNA(res$LC1)

## mfcl
summary(res$mfcl)

## number of top species
summary(res$S_top)

## mean trophic level
summary(res$TL.mean)

## number of sites
length(unique(res$ID.new.Bioregion))


# RUN GAMMS #####################################################


library(mgcv)
Indices.list = c("S","C","L.S", "L", "modularity.2", "mfcl", "omnivory.1", "Gen", 
                 "Vul", "TL.mean", "S_top", "S_basal","S_intermediate","Ftop","Fbasal",
                 "FInt","SDGen","SDVul","mass_basal","mass_top","mass_int")

source('Script/FunctionsAnalysis_MixedModels_RandomSlopeControl.R')
set.seed(1997)
## 22/08/2024 - running without random intercept
# res = add_context_var(res)
res$LC2 = as.factor(res$LC2)
res$LC = as.factor(res$LC)
res$ID.new.Bioregion = factor(res$ID.new.Bioregion)
setDF(res)


outputs = run_GAMMs(res_resampled = res, Indices.list = Indices.list)
ranef = outputs[[1]]
fixef = outputs[[2]]

write.csv(ranef, 'outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/randomEffect-FullDataset-RandomSlope-FixedSmooths.csv')
write.csv(fixef, 'outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/fixedEffect-FullDataset-RandomSlope-FixedSmooths.csv')

