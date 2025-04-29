### Code for Protected areas with clear management goals enhance avian food webs
################################################################################
## PREPARING ALL CONTEXT VARIABLES
## (1) extract env variable for all grid cells
## (2) save into csv
## (3) merge into an .RData
################################################################################

## Bioregion, land cover, remoteness, elevation, THPI, HP, Human density

setwd('F:/TheseSwansea/WDPA')

library(data.table)
library(sf)
library(raster)

europeRaster <- terra::rast(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")
cells_info <- foreign::read.dbf('data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img.vat.dbf')
europeRaster1 <- raster::raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")
cells_info = cbind(cells_info, raster::coordinates(europeRaster1)[which(!is.na(raster::values(europeRaster1))),])

################################################################################
#### adding bioregions

bioregions <- st_read("data/Galiana2021_network-area-europe-master/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp")
bioregions = st_transform(bioregions, crs = europeRaster@crs)
# bioregions_buff = st_buffer(bioregions, dist = 2000)

sf_master = st_as_sf(cells_info, coords = c('x','y'), crs = europeRaster@crs)

## extracting the bioregion for each point
bioregions_grids = st_intersects(sf_master,bioregions)
names(bioregions_grids) = sf_master$Value
bioregions_grid.dt = data.frame(PK_UID = do.call(rbind, bioregions_grids))
bioregions_grid.dt$Value = rownames(bioregions_grid.dt)
bioregions_grid.dt = merge(bioregions_grid.dt, st_drop_geometry(bioregions[,c('PK_UID', 'name','short_name')]))

setdiff(bioregions_grid.dt$Value, cells_info$Value)

write.csv(bioregions_grid.dt, 'outputs/Env.Grid/Grid.Bioregions.csv')


library(ggplot2)
ggplot() + geom_point(data = merge(bioregions_grid.dt, cells_info), aes(x = x, y = y, colour = short_name))

## region nb 10 is 'outside' so should be NA
# bioregions_grids[which(bioregions_grids$PK_UID == 10),] = NA


################################################################################
################################################################################
#### adding landcover

CLC = raster('E:/TheseSwansea/Bioclim/CorineLandCover2018/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif')
rasterOptions(progress = 'text')


plot(europeRaster)
plot(CLC, add = T)

ext = c(extent(9e5,4.4e6,9e5,2.4e6), extent(9e5,4.4e6,2.4e6,5.5e6), extent(4.4e6,6e6,9e5,2.4e6),extent(4.4e6,6e6,2.4e6,5.5e6),extent(6e6,7.4e6,9e5,2.4e6),extent(6e6,7.4e6,2.4e6,5.5e6))

for (i in 1:length(ext)){
  CLC1 = crop(CLC, ext[[i]])
  europeRaster1 = crop(europeRaster, ext[[i]])
  europePL1 = rasterToPolygons(europeRaster1)
  
  plot(CLC1)
  plot(europePL1, add = T)
  
  library(plyr)
  df_grid = raster::extract(CLC1, europePL1, exact = T, weights = T) 
  grid = lapply(df_grid, function(x) if (!is.null(x)) {ddply(data.frame(x), 'value', summarise, weight = sum(weight))} else NA)
  names(grid) = europePL1$reference_grid_10km
  grid = rbindlist(grid, idcol = T)
  
  write.csv(grid, paste0('outputs/Env.Grid/CLC/Grid.CLC.',i,'.csv'))
}

# transform the grid into a polygon
europePL = rasterToPolygons(europeRaster)
df_grid = raster::extract(CLC, europePL, exact = T, weights = T) 
save.image(file='outputs/Env.Grid/CLC_env_Grid.RData')

grid = lapply(df_grid, function(x) if (!is.null(x)) {ddply(data.frame(x), 'value', summarise, weight = sum(weight))} else NA)
names(grid) = europePL$reference_grid_10km
grid = rbindlist(grid, idcol = T)


################################################################################
################################################################################
#### adding remoteness

Remoteness = raster('E:/TheseSwansea/Bioclim/Remoteness/2015_accessibility_to_cities_v1.0/2015_accessibility_to_cities_v1.0.tif')
rasterOptions(progress = 'text')


plot(europeRaster)
Remoteness = crop(Remoteness, extent(-25,70,30,100))
plot(Remoteness)

europePL = rasterToPolygons(europeRaster)

ext = c(extent(-25,10,30,50), extent(10,50,30,50), extent(-25,10,50,100),extent(10,50,50,100),extent(50,70,30,50),extent(50,70,50,100))

for (i in c(1:length(ext))){
  
  Remoteness1 = crop(Remoteness, ext[[i]])
  Remoteness1 = projectRaster(Remoteness1, crs = europeRaster@crs)
  
  europePL1 = crop(europePL, extent(Remoteness1))
  plot(Remoteness1)
  plot(europePL1, add = T)
  
  library(plyr)
  df_grid = raster::extract(Remoteness1, europePL1, exact = T) 
  grid = lapply(df_grid, function(x) if (!is.null(x)) {Remoteness = mean(x, na.rm = T)} else NA)
  names(grid) = europePL1$reference_grid_10km
  grid1 = grid[which(!is.na(grid))]
  grid1 = rbindlist(lapply(grid1, function(x) data.frame(x)), idcol = T)
  names(grid1) = c('Value','Remoteness')
  
  write.csv(grid1, paste0('outputs/Env.Grid/Remoteness/Grid.Remoteness.',i,'.csv'))
}

grid1$Value = as.numeric(grid1$Value)
grid1 = merge(grid1, cells_info, by = 'Value', all.x = T)
library(ggplot2)
ggplot(data = grid1, aes(x = x, y = y, colour = Remoteness)) + geom_point() + scale_colour_gradientn(colours = rainbow(5))


################################################################################
################################################################################
#### adding elevation

El = raster('E:/Bioclim/SRTM_Elevation/wc2.1_30s_elev.tif')
Slope = terrain(El, opt = 'slope')
plot(El)
El@crs

# transform the grid into a polygon
europePL = rasterToPolygons(europeRaster)


ext = c(extent(-25,10,30,50), extent(10,50,30,50), extent(-25,10,50,100),extent(10,50,50,100),extent(50,70,30,50),extent(50,70,50,100))

rasterOptions(progress = 'text')
for (i in 1:length(ext)){
  
  El1 = crop(El, ext[[i]])
  El1 = projectRaster(El1, crs = europeRaster@crs)
  
  europePL1 = crop(europePL, extent(El1))
  
  plot(El1)
  plot(europePL1, add = T)
  df_grid1 = raster::extract(El1, europePL1, exact = T)
  
  library(plyr)
  fgrid = lapply(df_grid1, FUN = function(x) {m.elev = mean(x, na.rm = T); sd.elev = sd(x, na.rm = T);return(data.frame(m.elev,sd.elev))})
  names(fgrid) = europePL1$reference_grid_10km
  grid = rbindlist(fgrid, idcol = T)
  grid = grid[-which(is.na(grid$m.elev)),]
  
  write.csv(grid, paste0('outputs/Env.Grid/Elevation/Grid.Elevation.',i,'.csv'))
}

ext = c(extent(-25,10,30,50), extent(10,50,30,50), extent(-25,10,50,100),extent(10,50,50,100),extent(50,70,30,50),extent(50,70,50,100))

rasterOptions(progress = 'text')
for (i in 1:length(ext)){
  
  Slope1 = crop(Slope, ext[[i]])
  Slope1 = projectRaster(Slope1, crs = europeRaster@crs)
  
  europePL1 = crop(europePL, extent(Slope1))
  
  plot(Slope1)
  plot(europePL1, add = T)
  df_grid1 = raster::extract(Slope1, europePL1, exact = T)
  
  library(plyr)
  fgrid = lapply(df_grid1, FUN = function(x) {m.slope = mean(x, na.rm = T);return(data.frame(m.slope))})
  names(fgrid) = europePL1$reference_grid_10km
  grid = rbindlist(fgrid, idcol = T)
  grid = grid[-which(is.na(grid$m.slope)),]
  
  write.csv(grid, paste0('outputs/Env.Grid/Slope/Grid.Slope.',i,'.csv'))
  
}





################################################################################
################################################################################
#### adding human pressure

gHP = raster('E:/Bioclim/GlobalHumanPressure/gHM/gHM.tif')
gHP = projectRaster(gHP, crs = europeRaster@crs)
plot(gHP)

# transform the grid into a polygon
europePL = rasterToPolygons(europeRaster)


ext = c(extent(9e5,4.4e6,9e5,2.4e6), extent(9e5,4.4e6,2.4e6,5.5e6), extent(4.4e6,6e6,9e5,2.4e6),extent(4.4e6,6e6,2.4e6,5.5e6),extent(6e6,7.4e6,9e5,2.4e6),extent(6e6,7.4e6,2.4e6,5.5e6))

rasterOptions(progress = 'text')
for (i in 1:length(ext)){
  
  gHP1 = crop(gHP, ext[[i]])
  
  plot(europePL)
  plot(gHP1, add = T)
  df_grid1 = raster::extract(gHP1, europePL, exact = T)
  
  library(plyr)
  fgrid = lapply(df_grid1, FUN = function(x) {gHP = mean(x, na.rm = T); sd.gHP = sd(x, na.rm = T);return(data.frame(gHP,sd.gHP))})
  names(fgrid) = europePL$reference_grid_10km
  grid = rbindlist(fgrid, idcol = T)
  grid = grid[-which(is.na(grid$gHP)),]
  
  write.csv(grid, paste0('outputs/Env.Grid/gHP/Grid.gHP.',i,'.csv'))
  
}


################################################################################
################################################################################
#### adding human density
## nb of persons per km2 m

HD = raster('E:/Bioclim/PopulationDensity/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec_ClippedEurope.tif')
HD = projectRaster(HD, crs = europeRaster@crs)
plot(HD)

# transform the grid into a polygon
europePL = rasterToPolygons(europeRaster)


ext = c(extent(9e5,4.4e6,9e5,2.4e6), extent(9e5,4.4e6,2.4e6,5.5e6), extent(4.4e6,6e6,9e5,2.4e6),extent(4.4e6,6e6,2.4e6,5.5e6),extent(6e6,7.4e6,9e5,2.4e6),extent(6e6,7.4e6,2.4e6,5.5e6))

rasterOptions(progress = 'text')
for (i in 1:length(ext)){
  
  HD1 = crop(HD, ext[[i]])
  
  plot(europePL)
  plot(HD1, add = T)
  df_grid1 = raster::extract(HD1, europePL, exact = T)
  
  library(plyr)
  fgrid = lapply(df_grid1, FUN = function(x) {HD = mean(x, na.rm = T); sd.HD = sd(x, na.rm = T);return(data.frame(HD,sd.HD))})
  names(fgrid) = europePL$reference_grid_10km
  grid = rbindlist(fgrid, idcol = T)
  grid = grid[-which(is.na(grid$HD)),]
  
  write.csv(grid, paste0('outputs/Env.Grid/HD/Grid.HD.',i,'.csv'))
  
}


################################################################################
################################################################################
#### adding climate data

# BIO1 = Annual Mean Temperature
# 
# BIO5 = Max Temperature of Warmest Month
# 
# BIO6 = Min Temperature of Coldest Month
# 
# BIO12 = Annual Precipitation
# 
# BIO13 = Precipitation of Wettest Month
# 
# BIO14 = Precipitation of Driest Month
# 
# BIO15 = Precipitation Seasonality (Coefficient of Variation)

terra::terraOptions(progress = 0)
BIO1 = terra::rast('F:/Bioclim/WorldClim/wc2.1_30s_bio_1.tif')
BIO1 = terra::project(BIO1,europeRaster)

# transform the grid into a polygon
europePL = terra::as.polygons(europeRaster)

plot(europePL)
plot(BIO1, add = T)
df_grid1 = terra::extract(BIO1, europePL)
  
names(df_grid1) = c("Value","temp")
grid = df_grid1[-which(is.na(df_grid1$temp)),]
  
write.csv(grid, paste0('outputs/Env.Grid/Temperature/Grid.temp.csv'))
  

#### precipitation
BIO12 = terra::rast('F:/Bioclim/WorldClim/wc2.1_30s_bio_12.tif')
BIO12 = terra::project(BIO12,europeRaster)

plot(europePL)
plot(BIO12, add = T)

df_grid1 = terra::extract(BIO12, europePL)

names(df_grid1) = c("Value","precipitation")
grid = df_grid1[-which(is.na(df_grid1$precipitation)),]

write.csv(grid, paste0('outputs/Env.Grid/Precipitation/Grid.precipitation.csv'))



###########################################################################################################

###########################################################################################################
### Create .RData with all env data sets



setwd('F:/TheseSwansea/WDPA')
library(data.table)

Bioregions = read.csv('outputs/Env.Grid/Grid.Bioregions.csv'); names = names(Bioregions)
names(Bioregions)[2:5] = c('Bioregion.ID',"Value",'Bioregion','Bioregions_shortname')

### THPI : changes in human pressure from 1990 to 2010 (based on human density,
### land use change, nightlight) - 10km2 at the equator
### -100 to 100, where positive values mean increased human pressure and negative values mean decreased human pressure.

## gHP
gHP = rbindlist(lapply(paste0('outputs/Env.Grid/gHP/', list.files('outputs/Env.Grid/gHP/')), read.csv))  ; gHP = gHP[,-1]
names(gHP)[1] = c('Value')
gHP = ddply(gHP,'Value', summarise, m.gHP = mean(gHP), sd.gHP = mean(sd.gHP))

## Elevation
Elevation = rbindlist(lapply(paste0('outputs/Env.Grid/Elevation/', list.files('outputs/Env.Grid/Elevation/')), read.csv)) ; Elevation = Elevation[,-1]
names(Elevation)[1] = c('Value')
Elevation = ddply(Elevation,'Value', summarise, m.elev = mean(m.elev), sd.elev = mean(sd.elev))

## Slope
Slope = rbindlist(lapply(paste0('outputs/Env.Grid/Slope/', list.files('outputs/Env.Grid/Slope/')), read.csv)) ; Slope = Slope[,-1]
names(Slope)[1] = c('Value')
Slope = ddply(Slope,'Value', summarise, m.slope = mean(m.slope))

## Remoteness
Remoteness = rbindlist(lapply(paste0('outputs/Env.Grid/Remoteness/', list.files('outputs/Env.Grid/Remoteness/')), read.csv)) ; Remoteness = Remoteness[,-1]
Remoteness = ddply(Remoteness,'Value', summarise, Remoteness = mean(Remoteness))

## Land cover - need to add label names
CLC = rbindlist(lapply(paste0('outputs/Env.Grid/CLC/', list.files('outputs/Env.Grid/CLC/')), read.csv)); CLC = CLC[,-1]
names(CLC) = c('Value','GRID_CODE','lc.weight')
CLC.legend = read.csv('E:/Bioclim/CorineLandCover2018/u2018_clc2018_v2020_20u1_raster100m/Legend/CLC_legend.csv', sep=',')
CLC = merge(CLC, CLC.legend, all.x = T)
CLC1 = pivot_wider(CLC[,c('Value','LABEL3','lc.weight')], names_from = LABEL3,values_from = lc.weight, values_fill = 0, names_prefix = 'lc.')
CLC2 = pivot_wider(ddply(CLC[,c('Value','LABEL2','lc.weight')], c("Value","LABEL2"), summarise, lc.weight = sum(lc.weight, na.rm = T)), names_from = LABEL2,values_from = lc.weight, values_fill = 0, names_prefix = 'lc2.')
CLC3 = pivot_wider(ddply(CLC[,c('Value','LABEL1','lc.weight')], c("Value","LABEL1"), summarise, lc.weight = sum(lc.weight, na.rm = T)), names_from = LABEL1,values_from = lc.weight, values_fill = 0, names_prefix = 'lc1.')

anyDuplicated(CLC$Value)

## temperature
# temp = rbindlist(lapply(paste0('outputs/Env.Grid/Temperature/', list.files('outputs/Env.Grid/Temperature/')), read.csv), fill = T)
temp = read.csv('outputs/Env.Grid/Temperature/Grid.temp.csv')

# precipitation
precipitation = read.csv('outputs/Env.Grid/Precipitation/Grid.precipitation.csv')


save.image("outputs/ContextVariables.Over5yo.RData")




