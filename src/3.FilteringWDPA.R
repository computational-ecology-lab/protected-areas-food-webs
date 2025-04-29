################################################################################
### Code for Protected areas with clear management goals enhance avian food webs
################################################################################
### Load and filter WDPA dataset
### Calculate distance between protected areas and creating sites
### Calculate distance of grid cells to PA boundaries
################################################################################

setwd('E:/TheseSwansea/WDPA')


library(sf)
library(data.table)
library(stringr)
library(raster)
library(terra)
library(ggplot2)
library(rgdal)
library(plyr)
library(exactextractr)

europeRaster <- terra::rast(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")
cells_info <- foreign::read.dbf('data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img.vat.dbf')
cells_info = cbind(cells_info, raster::coordinates(raster::raster(europeRaster))[which(!is.na(raster::values(raster::raster(europeRaster)))),])

# europeRaster <- raster::raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")
# cells_info <- foreign::read.dbf('data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img.vat.dbf')
# cells_info = cbind(cells_info, coordinates(europeRaster)[which(!is.na(raster::values(europeRaster))),])

# load("E:/TheseSwansea/WDPA/Analysis/1.FilteringWDPA.Env.RData")


## keep only sites "designated", "inscribed" or "established" and polygons
## remove marine PA and man and biosphere reserves (include buffers and transition zones not actually protected)
## remove PAs smaller than 10 km2


## the WDPA file was divided into 3 sub files:
files = paste0("data/",str_subset(list.files("data/",recursive = T),'polygons.shp'))

# function that reads in raw WDPA files and applies first filters, and crops 
# WDPA to Europe 
readPA = function(file){
  
  f = st_read(file)
  f = subset(f, MARINE != 2) #0 (predominantly or entirely terrestrial), 1 (Coastal: marine and terrestrial), and 2 (predominantly or entirely marine). 
  f = subset(f, DESIG_ENG != "UNESCO-MAB Biosphere Reserve")
  f = subset(f, STATUS %in% c("Designated","Established","Inscribed"))
  f = subset(f, GIS_AREA >= 10) # larger than 10km2
  f$file = file
  
  # st_write(f, paste0('outputs/WDPA.Shp/2.FilteringWDPA/First_filter_',str_split(file, '/')[[1]][3], '.shp'), append=TRUE)
  
  print(paste('Initial filtering',nrow(f)))
  
  # ## Extract european grid values
  # #  A cell is covered if its center is inside the PA polygon
  # ext = raster::extract(terra::project(europeRaster, crs(f)), f)
  # ## match the ID from the extract function
  # f$ID = 1:nrow(f)
  # f_grid = merge(st_drop_geometry(f), ext, by = 'ID')
  # 
  # # saveRDS(f_grid, paste0('outputs/WDPA.Grid/2.FilteringWDPA/PA_grid_values_',str_split(file, '/')[[1]][3],'.RDS'))
  # 
  # ext = subset(ext, !is.na(PageName))
  # grid_list = plyr::ddply(ext, 'ID', summarise, PageName_list = list(PageName))
  # 
  # f = subset(f, ID %in% grid_list$ID)
  # 
  # st_write(f, paste0('outputs/WDPA.Shp/2.FilteringWDPA/First_filter_',str_split(file, '/')[[1]][3], '.shp'), append=TRUE)
  # print(paste('Europe filtering - number of rows:',nrow(f))) 
  
  return(f)
}

# around 5 mins
PA.crop = rbindlist(lapply(files, readPA), fill = T) 
PA.crop = st_as_sf(PA.crop)
# save.image("outputs/WDPA.Grid/PA.RData")

## Simple feature collection with 85795 features and 30 fields >> "Initial filtering 16447"
## "Initial filtering 15461"
## "Initial filtering 10279"

################################################################################################"
## CREATE SIMPLE POLYGONS FROM MULTIPOLYGONS

# # if need to reload the dataset: merge together the 3 files (Very fast)
# PA.crop = rbindlist(lapply(paste0("outputs/WDPA.Shp/2.FilteringWDPA/",
#                                   str_subset(str_subset(list.files("outputs/WDPA.Shp/2.FilteringWDPA/"),'\\.shp'),'_10km2')), st_read))


### transform into individual polygons instead of multipolygons (10 mins)
PA.crop.pol = data.frame()
for(i in 1:nrow(PA.crop)){
  
  if(i%%50 == 0){ cat(i) } # print progress
  
  P = PA.crop[i,]
  if (st_geometry_type(P) == 'GEOMETRYCOLLECTION'){
    P = st_collection_extract(P)
  }
  # cast to 'polygon' type
  PA.crop.pol = rbind(PA.crop.pol, st_cast(P, to = 'POLYGON'))
}
# compute area of individual polygon (<= GIS_AREA above)
PA.crop.pol$area = st_area(PA.crop.pol)
st_write(PA.crop.pol[,-which(names(PA.crop.pol) == 'PA.grid')], 'outputs/WDPA.Shp/2.FilteringWDPA/PAcrop.IndividualPolygons.wArea.shp')

print(paste('Individual poly area filtering - number of rows (polygons):',nrow(PA.crop.pol)))
print(paste('Individual poly area filtering - number of WDPAID:',length(unique(PA.crop.pol$WDPAID))))

################################################################################################"
## EXTRACTING GRID CELL IDs FOR EACH POLYGON

### extract grid cell IDs from europeRaster for the PA.crop with no multipolygons
PA.crop.values = PA.crop.pol
PA.crop.values$PA.grid = raster::extract(europeRaster, PA.crop.pol, exact = F, normalizeWeights=FALSE) # normalize weights makes that sum of weights is equal to 1 for each polygon
# exact = F take into account polygons that are only partly overlayed w/ europe Raster
save.image("outputs/WDPA.Shp/2.FilteringWDPA/PAcrop.IndividualPolygons.wArea.1.Merged.10km.PA.grid.RData")

# load("outputs/WDPA.Shp/2.FilteringWDPA/PAcrop.IndividualPolygons.wArea.1.Merged.10km.PA.grid.RData")
PA_grid = rbindlist(apply(PA.crop.values, 1, function(z) if (length(unlist(z$PA.grid))==0|all(is.na(unlist(z$PA.grid)))) data.frame(ID = z$ID, Value = NA_integer_) else data.frame(ID = rep(z$ID,length(unlist(z$PA.grid))),Value = unlist(z$PA.grid))))
saveRDS(PA_grid, 'outputs/WDPA.Grid/2.FilteringWDPA/PA_10km2_individualPolygons_grid_values.RDS')

### create a new ID for the newly created indiv polygons and save shapefile
anyDuplicated(PA.crop[,c('WDPAID','area')])
# PA.crop = PA.crop[-which(duplicated(PA.crop[,c('WDPAID','area')])),]
PA.crop$ID = paste0(PA.crop$WDPAID,'.',with(PA.crop, ave(as.numeric(WDPAID),WDPAID, FUN = seq_along)))
### save shapefile
st_write(PA.crop, 'outputs/WDPA.Shp/2.FilteringWDPA/Merged_PA_IndividualPolygons_10km2.shp')


# PA_grid = readRDS('outputs/WDPA.Grid/2.FilteringWDPA/PA_10km2_individualPolygons_grid_values.RDS')
# PA.crop = st_read('outputs/WDPA.Shp/2.FilteringWDPA/Merged_PA_IndividualPolygons_10km2.shp')
names(PA.crop)[which(names(PA.crop) == 'ID')] = 'WDPAID.new'

summary(PA_grid$PageName) # NAs are cells that fall inside a protected area but do not have a PageName in europeRaster 
length(unique(PA_grid$WDPAID.new)) 

## the missing values are cells that are not in the europeRaster
PA_missing = PA_grid[which(is.na(PA_grid$PageName)),] # for plotting below

# remove those PAs that don't overlap grid cells
PA_grid = merge(PA_grid, cells_info, by = 'PageName')

## looking at what the map looks like:
g = ggplot() + geom_sf(data = PA.crop, colour = ifelse(PA.crop$WDPAID.new %in% PA_missing$WDPAID.new, 'red', 'green')) + 
  geom_sf(data = st_transform(st_as_sf(PA_grid, coords = c('x','y'), crs = crs(europeRaster)), crs(PA.crop)), size = 0.1) 
ggsave('Plots/temp.Figures/PA.Map.initial.png', plot = g)

length(unique(PA_grid$WDPAID.new)) 


################################################################################
# remove PAs that were designated after 2015 
PA_grid = PA_grid[-which(PA_grid$STATUS_YR > 2015),] 
length(unique(PA_grid$WDPAID.new))  

################################################################################
## THERE ARE GRID CELLS WHICH HAVE SEVERAL 'LAYERS' OF PROTECTION. 


# To only keep one we keep the highest level of protection for each grid cell 
# which is overlapped by several PAs

length(which(duplicated(PA_grid))) 
PA_grid = PA_grid[-which(duplicated(PA_grid)),] # remove any rows that are completely identical

length(unique(PA_grid$WDPAID.new)) 
# this should decrease (we will be removing overlapping PAs)
length(unique(PA_grid$Value)) 
# this should be the same at the end of the process 

## first, we save the initial number of layers of protection per grid cell
PA_grid = merge(PA_grid, ddply(PA_grid, 'Value', summarise, layers_protection = length(Value), 
                Directed_at_birds_BeforeFiltering = ifelse(any(DESIG_ENG %in% c("Ramsar Site, Wetland of International Importance", "Special Protection Area (Birds Directive)")),
                                                           1, 0)), by = 'Value')
summary(PA_grid$layers_protection) # from 1 to 8 layers of protection
summary(factor(PA_grid$Directed_at_birds_BeforeFiltering)) 

dups.list = unique(PA_grid$Value[duplicated(PA_grid$Value)]) 
dups = PA_grid[which(PA_grid$Value %in% dups.list),]
PA.nodups = PA_grid[-which(PA_grid$Value %in% dups.list),]

length(unique(dups$Value)) + length(unique(PA.nodups$Value)) == length(unique(PA_grid$Value))
length(unique(c(PA.nodups$WDPAID.new, dups$WDPAID.new))) == length(unique(PA_grid$WDPAID.new))

######## three layers of filtering:
######## 1. IUCN PROTECTION, 2. AGE, 3. SIZE

# scoring system where less is better
levels = data.frame(IUCN_CAT = levels(factor(PA_grid$IUCN_CAT)), score = c(1,2,3,4,5,10,10,10,6,7))
dups = merge(dups, levels, by = 'IUCN_CAT')
dups = merge(dups, ddply(dups, c('Value'), summarise, to_keep_score = min(score), to_keep_ID = ID[which.min(score)]), by = 'Value')
length(unique(dups$Value)) + length(unique(PA.nodups$Value)) == length(unique(PA_grid$Value))

# filter by IUCN Protection level
dups.filtered = dups[which(dups$score == dups$to_keep_score),]
dups.list = dups.filtered$Value[which(duplicated(dups.filtered$Value))]
to_keep = dups.filtered[which(!dups.filtered$Value %in% dups.list),]
dups = dups.filtered[which(dups.filtered$Value %in% dups.list),]

length(unique(to_keep$Value)) + length(unique(dups$Value)) + length(unique(PA.nodups$Value)) == length(unique(PA_grid$Value))

setDF(to_keep)
PA.nodups = rbind(PA.nodups, to_keep[,names(PA.nodups)])

# filter by year of designation
dups.list = unique(PA_grid$Value[which(!PA_grid$Value %in% PA.nodups$Value)]) # 1893 grid cells
dups = PA_grid[which(PA_grid$Value %in% dups.list),]
length(unique(dups$Value)) + length(unique(PA.nodups$Value)) == length(unique(PA_grid$Value))

dups.Noyear = unique(dups$Value[which(dups$STATUS_YR == 0)]) # put aside cells which have no info on year of implementation
dups.Noyear.dt = dups[which(dups$Value %in% dups.Noyear),] ; length(unique(dups.Noyear.dt$Value))
dups = dups[-which(dups$Value %in% dups.Noyear),] ; length(unique(dups$Value))

length(unique(dups$Value)) + length(unique(dups.Noyear.dt$Value)) + length(unique(PA.nodups$Value)) == length(unique(PA_grid$Value))

dups = merge(dups, ddply(dups, 'Value', summarise, to_keep_yr = min(STATUS_YR)), by = 'Value') ; length(unique(dups$Value))
dups.filtered = dups[which(dups$STATUS_YR == dups$to_keep_yr),] ; length(unique(dups.filtered$Value))
dups.list = unique(dups.filtered$Value[which(duplicated(dups.filtered$Value))])
dups = dups.filtered[which(dups.filtered$Value %in% dups.list),]
to_keep = dups.filtered[which(!dups.filtered$Value %in% dups.list),]

length(unique(dups$Value)) + length(unique(to_keep$Value)) + length(unique(dups.Noyear.dt$Value)) + length(unique(PA.nodups$Value)) == length(unique(PA_grid$Value))

setDF(dups)
dups = rbind(dups[,names(dups.Noyear.dt)], dups.Noyear.dt)
length(unique(dups$Value))

length(unique(to_keep$Value)) + length(unique(dups$Value)) + length(unique(PA.nodups$Value)) 
anyDuplicated(to_keep$Value) ## should be false

setDF(to_keep)
PA.nodups = rbind(PA.nodups, to_keep[,names(PA.nodups)])

# filter by area
dups.list = unique(PA_grid$Value[which(!PA_grid$Value %in% PA.nodups$Value)]) 
dups = PA_grid[which(PA_grid$Value %in% dups.list),] ; length(unique(dups$Value))

dups = merge(dups, ddply(dups, 'Value', summarise, to_keep_area = max(area)), by = 'Value'); length(unique(dups$Value))
dups.filtered = dups[which(dups$area == dups$to_keep_area),]; length(unique(dups.filtered$Value))
dups.list = dups.filtered$Value[which(duplicated(dups.filtered$Value))]
dups = dups.filtered[which(dups.filtered$Value %in% dups.list),]
to_keep = dups.filtered[which(!dups.filtered$Value %in% dups.list),]

length(unique(dups$Value)) + length(unique(to_keep$Value)) + length(unique(PA.nodups$Value)) == length(unique(PA_grid$Value))
anyDuplicated(to_keep$Value)

setDF(to_keep)
PA.nodups = rbind(PA.nodups, to_keep[,names(PA.nodups)])

# keep one at random among the rest
dups.list = unique(PA_grid$Value[which(!PA_grid$Value %in% PA.nodups$Value)])
dups = PA_grid[which(PA_grid$Value %in% dups.list),] ; length(unique(dups$Value))

to_keep = dups[!duplicated(dups$Value),]
setDF(to_keep)
PA.nodups = rbind(PA.nodups, to_keep[,names(PA.nodups)])
nrow(PA.nodups) == length(unique(PA_grid$Value))
length(unique(PA.nodups$WDPAID.new)) # 3299 - lost about half of the initial PAs

# check
anyDuplicated(PA.nodups$Value)
setdiff(PA.grid$Value,PA.nodups$Value)

saveRDS(PA.nodups, 'outputs/WDPA.Grid/Grid.NoDups.PA_10km_Europe.POLYGONS.Over5yo_2023.RDS')
sf::st_write(subset(PA.crop, WDPAID.new %in% PA.nodups$WDPAID.new), 
             'outputs/WDPA.Shp/2.FilteringWDPA/Grid.NoDups.PA_10km_Europe.POLYGONS.Over5yo_2023.shp')


################################################################################
## DISTANCE OF ALL GRID CELLS WITHIN 200 KMs TO PA 

PA.grid.nodups = readRDS('outputs/WDPA.Grid/Grid.NoDups.PA_10km_Europe.POLYGONS.Over5yo_2023.RDS')
PA.crop = st_read('outputs/WDPA.Shp/2.FilteringWDPA/Grid.NoDups.PA_10km_Europe.POLYGONS.Over5yo_2023.shp')

seq = seq(0, nrow(PA.crop), 500) # separate PA shapefile into manageable sub units
library(nngeo)

## yields the list and distance of the 10 000 neighbours within 200 km of each protected area 
for (i in 1:(length(seq)-1)){
  start = seq[i] + 1
  stop = seq[i + 1]
  PA.sub = PA.crop[start:stop,]
  PA.sub = st_transform(PA.sub, crs = crs(europeRaster))
  print(start)
  dist = st_nn(PA.sub,st_as_sf(cells_info, coords = c('x','y'), crs = crs(europeRaster)), k = 10000, maxdist = 200e3, returnDist = T, parallel = 4)
  
  names(dist$nn) = names(dist$dist) = PA.sub$ID
  dt = do.call(rbind, Map(data.frame, Value=dist$nn, dist.PA=dist$dist))
  dt$ID = unlist(lapply(rownames(dt), function(x) paste0(str_split(x,'\\.')[[1]][1:2], collapse = '.')))
  
  saveRDS(dt, file = paste0('outputs/WDPA.Grid/distance/NearestNeighbours/NN_PA_POLYGONS_k10e3k_200km.',i,'.rds'))
  
  closeAllConnections()
  gc()
}



## function that take as argument the output of st_nn above and keeps the smallest distance
## between the PA and its neighbouring PAs
fun_distPA = function(start, stop, d.PA.mat = dist.PA.mat, d.PA = dist.PA, P.grid = PA.grid){
  
  library(foreach)
  library(doParallel)
  library(doSNOW) # for progression bar
  cores <- 2
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  
  iterations <- stop - start
  # iterations = 20
  
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # start_time <- Sys.time()
  dist <- foreach(P = start:stop, .combine=rbind, .options.snow = opts) %dopar%{
    library(plyr)
    #for (P in 1:iterations){
    #cat(P)
    IDP = unique(d.PA$ID)[P]
    # get all the distances from the grid cells that constitute the PA (ID)
    sub.dist.PA = subset(d.PA, Value %in% P.grid$Value[which(P.grid$ID == IDP)])
    
    # only keep the smallest distance to each neighbouring PA
    sub.dist.PA.nb = ddply(sub.dist.PA, "ID", summarise, min.dist = min(dist.PA), min.Value = Value[which.min(dist.PA)])
    
    if(nrow(sub.dist.PA.nb)>0){  
      sub.dist.PA.nb$ID.source = IDP
      print(sub.dist.PA.nb)
    }
  }
  # end_time <- Sys.time()
  # end_time - start_time
  return(dist)
  
  close(pb)
  stopCluster(cl) 
  
} 
### end function


dt = rbindlist(lapply(paste0("outputs/WDPA.Grid/distance/NearestNeighbours",
                             str_subset(list.files("outputs/WDPA.Grid/distance/NearestNeighbours"), 'rds')), readRDS))
saveRDS(dt, file = paste0('outputs/WDPA.Grid/distance/NearestNeighbours/NN_PA_POLYGONS_k10e3k_200km.rds'))

## distance to each neighbouring PA
dist.dt = fun_distPA(start = 1, stop = length(unique(dt$ID)), d.PA = dt)
saveRDS(dist.dt, "outputs/WDPA.Grid/distance/distancePA.dist.400km.10e3k.POLYGONS.RDS")


dist.dt = readRDS("outputs/WDPA.Grid/distance/PA.dist.200km.10e3k.POLYGONS.RDS")
#missing = setdiff(PA.crop$ID, unique(dist.dt$ID))


plot(PA.crop[which(PA.crop$ID %in% missing),1])
## all the missing PAs are outside our study zone

# create a distance matrix
# input the minimum dist into the distance matrix
dist.PA.mat = matrix(nrow = length(unique(dt$ID)), ncol = length(unique(dt$ID)))
rownames(dist.PA.mat) = colnames(dist.PA.mat) = unique(dt$ID)
diag(dist.PA.mat) = 0

for (P in unique(dt$ID)){
  cat(P)
  s.dist = subset(dist.dt, ID.source == P)
  dist.PA.mat[match(s.dist$ID,rownames(dist.PA.mat)),P] = s.dist$min.dist
}

dist.PA.mat = dist.PA.mat[-which(rowSums(dist.PA.mat, na.rm = T) + colSums(dist.PA.mat, na.rm = T) == 0),-which(rowSums(dist.PA.mat, na.rm = T) + colSums(dist.PA.mat, na.rm = T) == 0)]
saveRDS(dist.PA.mat, paste0("outputs/WDPA.Grid/distance/PA.distMatrix.Asym.400km.10e3kNN.RDS"))

dist.PA.mat = readRDS(paste0("outputs/WDPA.Grid/distance/PA.distMatrix.Asym.400km.10e3kNN.RDS"))


# replace with the correct (smallest) distance inside the matrix so that the matrix is symmetric
N.PA = length(colnames(dist.PA.mat))
for(i in 1:N.PA){
  for(j in 1:N.PA){
    
    if(!all(is.na(c(dist.PA.mat[i,j], dist.PA.mat[j,i])))){
      dist.PA.mat[i,j] = min(dist.PA.mat[i,j], dist.PA.mat[j,i], na.rm = T)
    }
  }}

# save.image("E:/TheseSwansea/WDPA/Analysis/Grids/PA.dist.400km.10e3k.RData")
saveRDS(dist.PA.mat, "outputs/WDPA.Grid/distance/PA.distMatrix.400km.10e3kNN.RDS")


ggplot() + geom_sf(data = st_transform(PA[which(PA$WDPAID == ID),], crs = crs(europeRaster)), aes(fill = factor(WDPAID)), alpha = 0.5) + 
  geom_point(data = test, aes(x = x, y = y))

save.image("E:/TheseSwansea/WDPA/Analysis/1.FilteringWDPA.Env.distance.RData")



################################################################################
## Merging together close (<1km) PAs 

# from 1.FilteringWDPA.GridCreation.LAB2.R
dist.PA.mat = readRDS("outputs/WDPA.Grid/distance/PA.distMatrix.400km.10e3kNN.RDS")
dim(dist.PA.mat)

# keep only non overlapping polygons within the study area
PA.grid.nodups = readRDS('outputs/WDPA.Grid/Grid.NoDups.PA_10km_Europe.POLYGONS.Over5yo.RDS')
dist.PA.mat = dist.PA.mat[-which(!rownames(dist.PA.mat) %in% PA.grid.nodups$ID), 
                          -which(!colnames(dist.PA.mat) %in% PA.grid.nodups$ID)]

library(data.table)
# choose a threshold (10km is too large)
thresh = 1e3 # km
grid.dist1 = apply(dist.PA.mat, 1, FUN = function(x) data.frame(rownames(dist.PA.mat)[which(x<thresh)]))
grid.dist = rbindlist(grid.dist1, idcol = T)
names(grid.dist) = c('ID','ID.nn')
## taking the most connected PAs and removing their nieghbouring PAs from the pool of PAs

# empty dataframe to store merged PAs - new PA ID
merged.PAs = data.frame(WDPAID = unique(rownames(dist.PA.mat)), WDPAID.new = NA)

library(plyr)
nn.neighbours = ddply(grid.dist, 'ID', summarise, nb = length(ID.nn))
x = nn.neighbours$nb; names(x) = nn.neighbours$ID
x = sort(x, decreasing = T)

# Merge PAs close together (1km threshold - see above)
ind = 1
while (length(x)>0) {
  
  ## take the most connected PA
  P = names(x)[1]
  cat(ind)
  
  ## look at its close neighbours in terms of distance 
  nn.10 = grid.dist$ID.nn[which(grid.dist$ID == P)]
  while(length(setdiff(grid.dist$ID.nn[which(grid.dist$ID %in% nn.10)], nn.10))>0){
    nn.10 = c(nn.10,grid.dist$ID.nn[which(grid.dist$ID %in% nn.10)])
    nn.10 = unique(nn.10)
  }
  
  # assign the new WDPAID 
  merged.PAs[which(merged.PAs$WDPAID %in% nn.10), 'WDPAID.new'] = ind
  
  # remove the PAs from the pool of potential PAs to draw from:
  x = x[-which(names(x) %in% nn.10)]
  x = sort(x, decreasing = T)
  
  grid.dist = grid.dist[-which(grid.dist$ID.nn %in% nn.10|grid.dist$ID %in% nn.10),]
  
  ind = ind + 1
  
}

min(dist.PA.mat[merged.PAs[which(is.na(merged.PAs$WDPAID.new)),'WDPAID'],], na.rm = T)
## PAs with NA value have no neighbours closer than 1km therefore their WDPAID.new is themselves

merged.PAs[which(is.na(merged.PAs$WDPAID.new)),'WDPAID.new'] = merged.PAs[which(is.na(merged.PAs$WDPAID.new)),'WDPAID']
anyNA(merged.PAs$WDPAID.new)
# average size of new PAs
names(merged.PAs) = c("ID","ID.new")
nb = ddply(merged.PAs, 'ID.new', summarise, x = length(ID))
summary(nb)

# save just the table
saveRDS(merged.PAs, "outputs/WDPA.Grid/NewID.MergedPAs.1km.20012023.Over5yo.RDS")

# save the shapefile
PA.grid.nodups = merge(PA.grid.nodups, merged.PAs, by = 'ID')
PA.grid.nodups = st_as_sf(PA.grid.nodups)
st_write(PA.grid.nodups, "outputs/WDPA.Shp/NewID.MergedPAs.1km.20012023.Over5yo.shp")



