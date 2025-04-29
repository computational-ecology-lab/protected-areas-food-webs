### Code for Protected areas with clear management goals enhance avian food webs



################################################################################
### Filter GBIF raw data
### Merging ebird and GBIF filtered data
### Filter on survey effort
### Calculating network metrics
################################################################################




## Content:

#####################################
## 1. Read and filter GBIF raw split data
#####################################
## 2 - assign site ID (and protection status) to each grid cell
#####################################
## 3 - create list of events and filter by date
####################################
## 4 - Compute number of visits per grid cell (= survey effort)
#####################################
## 5 - Filter on survey effort (>80) and at least 15 protected AND 15 non protected cells
#####################################
## 6 - compute network metrics 


## Load packages and functions -----
setwd("E:/TheseSwansea/WDPA")

library(data.table)
library(dplyr)
library(sf)
library(stringr)
library(rnaturalearth)
library(plyr)
library(tidyr)
library(ggplot2)

## grid of europe at a 10x10km resolution
europeRaster <- raster::raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")
cells_info <- foreign::read.dbf('data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img.vat.dbf')
cells_info = cbind(cells_info, raster::coordinates(europeRaster)[which(!is.na(raster::values(europeRaster))),])
source("data/Galiana2021_network-area-europe-master/whois-function.R")

europe = rnaturalearth::ne_countries(continent = "Europe", returnclass = "sf", scale = 'large')
europe = sf::st_transform(europe, crs = raster::crs(europeRaster))


## A. Read split files for each download key (= one polygon in Europe) + filter data, add grid cell info 

## 1 - assign PA ID + other contextual variables (distance to PA, bioregion etc) to each grid cell
## 2 - create list of events with number of visits (=events) per grid cell (= survey effort)
## 3 - subset grid cells with larger survey effort than threshold (=80)
## 4 - compute network metrics for each grid cell


################################################################################
## FILTERING OF GBIF DATA - 11/04/2022 -----------------------------------------
################################################################################

## A. Read split files for each download key (= one polygon in Europe) + filter data, add grid cell info 
## generate eventList

library(rgbif)

m_info = read.csv("data/Galiana2021_network-area-europe-master/species_codes_and_taxonomy.csv", sep = ";")
gbif_names = name_backbone_checklist(m_info$Species) 
gbif_taxon_keys <- gbif_names$usageKey

### Ebird data was downloaded straight from the ebird website, country by country
### all records since 2005, all bird sp

### list of downloaded GBIF files using first script (1.Download.GBIF.R):

## 0222825-210914110416597 # NM poly2 41,283,168  ## north east
## 0222868-210914110416597 ## poly1   9,792,753   ## north                
## 0222874-210914110416597 ## poly4   9,321,985   ## south west              
## 0222877-210914110416597 ## poly5   23,486,246  ## south middle
## 0222905-210914110416597 ## poly3.1 44,833,015  ## north middle west
## 0222907-210914110416597 ## poly3.2 55,818,583  ## north
## 0222933-210914110416597 ## poly6   629,279     ## south east
## 0282245-210914110416597 ## poly8               ## south east 1

## these files were split to load them into R without exhausting memory

################################################################################
### FUNCTION FOT READING AND FILTERING RAW SPLIT GBIF FILES and generating eventList
read.file = function(file, write = F, eventList = eventList){
  
  if(dir.exists(paste0("outputs/dta.Occurrence/GBIF/",cardinal,'/',key,'/Split'))){
    setwd(paste0("outputs/dta.Occurrence/GBIF/",cardinal,'/',key,'/Split'))
  } else{
    setwd(paste0("outputs/dta.Occurrence/GBIF/",cardinal,'/',key))
  }

  if(write){
    print(file)
    df = fread(file, nThread = 10)
    
    if(key == '0222907-210914110416597'){
      names(df) = c('V1',"gbifID","datasetKey","occurrenceID","class","order",'species','taxonRank','scientificName','countryCode', names(df)[c(10:17,19:25)])
    }
    
    ### https://data-blog.gbif.org/post/gbif-filtering-guide/
    
    df = df %>%
      filter(taxonRank == "SPECIES") %>%
      filter(occurrenceID != "") %>%
      filter(species != "") %>%
      filter(!str_detect(issue, "NEGATED|SWAPPED")) %>%
      filter(coordinateUncertaintyInMeters <= 10e3 | !is.na(coordinateUncertaintyInMeters))
    
    sf_df = st_as_sf(df, coords = c("decimalLongitude","decimalLatitude"), crs = 4326) ## spatial points for the occurrence data GBIF
    sf_df = st_transform(sf_df, crs = crs(europeRaster))
    
    # extract grid cell ID
    df_grid = raster::extract(europeRaster, sf_df)
    sf_df = cbind(sf_df,df_grid)
    
    eventList = unique(st_drop_geometry(sf_df[,c('eventDate','datasetKey','df_grid', 'class')]))
    print(dim(eventList))
    
    write.csv(eventList, paste0("outputs/dta.Occurrence/GBIF/EventList/EventList.", file))
    write.csv(st_drop_geometry(sf_df), paste0("outputs/dta.Occurrence/GBIF/Filtered.effort/Filtered.effort.masterValues.", file))
    # save species names for taxonomy
    write.csv(unique(st_drop_geometry(sf_df[,c('species','order','class','taxonKey')])), paste0("outputs/dta.Occurrence/GBIF/SpeciesList/SpeciesList.", file))
    
  } else { 
    
    sf_df = fread(paste0("outputs/dta.Occurrence/GBIF/Filtered.effort/Filtered.effort.masterValues.", file))
    
    if(cardinal == 'North'){ ## some keys still have the geometry column
      if(key == "0282225-210914110416597"){
        sf_df = sf_df[,-c(1)]
        names(sf_df) = c(names(sf_df)[2:51], "geometry")
      }else if(key == "0222868-210914110416597"|key == "0222905-210914110416597"){
        sf_df = sf_df
      }else{
        sf_df = sf_df[,-c(1,2)]
        names(sf_df) = c(names(sf_df)[2:25], "geometry")
      }
    }
    
    eventList = unique(sf_df[,c('eventDate','datasetKey','df_grid','class')])
    write.csv(eventList, paste0("outputs/dta.Occurrence/GBIF/EventList/EventList.", file))
    
  }
} ### end function



### RUN read.file function
all.south.keys = c('0222933-210914110416597','0282245-210914110416597','0222877-210914110416597','0222874-210914110416597')
all.north.keys = c('0222868-210914110416597','0222905-210914110416597','0222825-210914110416597','0222907-210914110416597','0282225-210914110416597')

for(key in all.south.keys){
  cardinal = "South"
  
  if(dir.exists(paste0("outputs/dta.Occurrence/GBIF/",cardinal,'/',key,'/Split'))){
    setwd(paste0("outputs/dta.Occurrence/GBIF/",cardinal,'/',key,'/Split'))
  } else{
    setwd(paste0("outputs/dta.Occurrence/GBIF/",cardinal,'/',key))
  }
  lapply(str_subset(list.files(),"csv"),read_file)
}


for(key in all.north.keys){
  cardinal = "North"
  
  if(dir.exists(paste0("outputs/dta.Occurrence/GBIF/",cardinal,'/',key,'/Split'))){
    setwd(paste0("outputs/dta.Occurrence/GBIF/",cardinal,'/',key,'/Split'))
  } else{
    setwd(paste0("outputs/dta.Occurrence/GBIF/",cardinal,'/',key))
  }
  
  lapply(str_subset(list.files(),"csv"),read_file)
}


#####################################
## 1 - assign PA ID to each grid cell ------------------------------------------



## Merged PA more than 5 years old
# merged PAs: ID.new from 1.FilteringWDPA.GridCreation.5yo (we merged PAs less than 1km from each others)
merged.PAs = readRDS("outputs/WDPA.Grid/NewID.MergedPAs.1km.20102023.Over5yo_st_distance.RDS")

Bioregions = read.csv('outputs/Env.Grid/Grid.Bioregions.csv'); names = names(Bioregions)
names(Bioregions)[2:5] = c('Bioregion.ID',"Value",'Bioregion','Bioregions_shortname')

## CREATING SITES (Portected cells + surrounding non protected cells)

## we label the non-protected surrounding grid cell so that they are assigned to a protected area:
## this column will be called NeighbouringPAID.new
## protected cells : ID.new == NeighbouringPAID.new, 
## non protected cells : ID.new = NA and a value for NeighbouringPAID.new

## the list and distance of the 10 000 neighbours within 200 km of each protected area 
dt = readRDS('outputs/WDPA.Grid/distance/NearestNeighbours/NN_PA_POLYGONS_k10e3k_200km_2023.rds') # load in the distance table, 
dt = dt[-which(dt$dist.PA>=100e3),] # keep only the cells within 100km of a PA

summary(dt)

# create a table of grid cells within 100 km of a PA:
dt.neighbours = merge(dt, unique(merged.PAs[,c("WDPAID.new","ID.new")]), by.x = 'ID', by.y = "WDPAID.new") 

# remove protected points (we want to locate non protected grid cells here)
dt.neighbours = dt.neighbours[-which(dt.neighbours$Value %in% merged.PAs$Value),]
# assign to each non protected grid cell the nearest protected area, and store its distance
dt.neighbours = ddply(dt.neighbours, c("Value"), summarise, dist.PA.min = min(dist.PA), ID = ID[which.min(dist.PA)], ID.new = ID.new[which.min(dist.PA)])

dt.neighbours$Value = as.numeric(dt.neighbours$Value)
names(dt.neighbours)[which(names(dt.neighbours) == "ID")] = "NeighbouringPAID" # rename so we don't confuse ID (only protected cells have an ID) and Neighbouring PAID (all cells close enough to a PA have a NeighbouringPAID)
names(dt.neighbours)[which(names(dt.neighbours) == "ID.new")] = "NeighbouringPAID.new"


################################################################################
### now we group together protected and non protected grid cells into sites:
### protected grid cells have both and neighbouringPAID and and an ID, whilst
### non-protected grid cells only have an neighbouring PAID

setdiff(dt.neighbours$NeighbouringPAID, merged.PAs$WDPAID.new) ; setdiff(merged.PAs$WDPAID.new, dt.neighbours$NeighbouringPAID)

## (2) add non protected grid cells to merged.PAs
merged.PAs = merge(merged.PAs[,-which(names(merged.PAs) %in% c('x','y','PageName','Count'))], cells_info, by = 'Value', all = T)
print(paste('Should be TRUE:',anyNA (merged.PAs$ID.new))) 

## (3) add distance info to merged.PAs
setdiff(dt.neighbours$Value,merged.PAs$Value)
merged.PAs = merge(merged.PAs, dt.neighbours, by = "Value", all = T)

# (4) assign a NeighbouringPAID to protected cells (their own ID)
merged.PAs$NeighbouringPAID.new[which(!is.na(merged.PAs$IUCN_CAT))] = merged.PAs$ID.new[which(!is.na(merged.PAs$IUCN_CAT))]
merged.PAs$NeighbouringPAID[which(!is.na(merged.PAs$IUCN_CAT))] = merged.PAs$WDPAID.new[which(!is.na(merged.PAs$IUCN_CAT))]

# summary(factor(merged.PAs[,c('NeighbouringPAID.new')]));summary(factor(merged.PAs[,c('NeighbouringPAID')]))
# we remove grid cells too far away from a PA
merged.PAs = merged.PAs[-which(is.na(merged.PAs$NeighbouringPAID.new)),]

# add bioregions info (grid cell level)
merged.PAs$Bioregion = Bioregions$Bioregion[match(merged.PAs$Value, Bioregions$Value)]
merged.PAs = merged.PAs[-which(is.na(merged.PAs$Bioregion)),]
merged.PAs$ID.new.Bioregion = apply(merged.PAs[,c('NeighbouringPAID.new','Bioregion')],1, FUN = function(x) paste(x, collapse = '_'))
merged.PAs$ID.new.Bioregion = str_replace(merged.PAs$ID.new.Bioregion, '   ','')

saveRDS(merged.PAs, "outputs/WDPA.Grid/GridCells.NeighbouringPAID.new.Over5yo_2023.RDS")
st_write(st_as_sf(merged.PAs, coords = c('x','y'), crs = crs(europeRaster)), "outputs/WDPA.Shp/4.EventList/GridCells.NeighbouringPAID.new.Over5yo.shp")


############################################################################################
## 2 - create list of events and remove those that are too old compared to the creation of the PA -----

merged.PAs = readRDS("outputs/WDPA.Grid/NewID.MergedPAs.1km.20012023.Over5yo_st_distance.rds")

# merged.PAs$STATUS_YR = PA.grid.nodups$STATUS_YR[match(merged.PAs$ID,PA.grid.nodups$ID)]
merged.PAs.youngestDate = ddply(merged.PAs, "ID.new", summarise, youngestDate = max(STATUS_YR))
merged.PAs = merge(merged.PAs, merged.PAs.youngestDate, by.y = "ID.new", by.x = "ID.new")
rm(merged.PAs.youngestDate)


library(sf)
## read in the event lists created above
eventList = data.table::rbindlist(lapply(paste0("outputs/EventList/",
                                                str_subset(list.files("outputs/EventList"),'csv')), read.csv), use.names=TRUE, fill = T)
eventList = eventList[,-c(1,6)]

setDT(eventList)
# Count repetitions for each category
eventList[datasetKey != 'eBird', origin := "GBIF"]; eventList[datasetKey == 'eBird',origin := "eBird"] 
class_counts <- eventList[, .N, by = .(class, origin)]
round(class_counts$N[which(class_counts$origin == 'GBIF')]/sum(class_counts$N[which(class_counts$origin == 'GBIF')]),3)*100
ggplot() + geom_bar(data = class_counts, aes(x = class, y = N, fill = origin), stat = 'identity') 

setDF(eventList)
eventList = subset(eventList, class == 'Aves')

saveRDS(eventList, "outputs/EventList/EventListMerged/EventList.GBIF.eb.effort80.birds.RDS")

## add bioregion
eventList = readRDS("outputs/EventList/EventListMerged/EventListMerged/EventList.GBIF.eb.effort80.birds.RDS")
eventList$Bioregion = merged.PAs[match(eventList$df_grid, merged.PAs$Value),'Bioregion']
eventList$ID.new.Bioregion = merged.PAs[match(eventList$df_grid, merged.PAs$Value),'ID.new.Bioregion']

nrow(eventList)
eventList = eventList[-which(is.na(eventList$Bioregion)),]
nrow(eventList)

# ## keep only birds for now
# eventList = subset(eventList, class == 'Aves'); dim(eventList)

## add protection info, date of creation and distance to PA
eventList$ID.new = merged.PAs$ID.new[match(eventList$df_grid,merged.PAs$Value)] # PA ID for protected cells only
eventList$NeighbouringPAID.new = merged.PAs$NeighbouringPAID.new[match(eventList$df_grid,merged.PAs$Value)] # site ID 
eventList$youngestDate = merged.PAs$youngestDate[match(eventList$df_grid,merged.PAs$Value)] 
eventList$dist.PA.min = merged.PAs$dist.PA.min[match(eventList$df_grid,merged.PAs$Value)]
eventList$dist.PA.min[which(!is.na(eventList$ID.new))] = 0 # protected grid cells have a distance of zero here
## keep only grid cells neighbouring a Protected Area
eventList = subset(eventList, !is.na(NeighbouringPAID.new)); dim(eventList)
summary(eventList$youngestDate) # should be no PAs after 2015

# two different formats for dates depending on whether the data is from ebird or GBIF 
## we extract the year of the event to remove those too old compared to the creation of the PA
eventList$Year = NA
# /!\ takes 20 mins to run
# GBIF
eventList$Year[which(!eventList$datasetKey %in% c('eBird'))] = as.numeric(sapply(eventList$eventDate[which(!eventList$datasetKey %in% c('eBird'))], function(x) unlist(str_split(x,"-"))[1]))
# ebird
eventList$Year[which(eventList$datasetKey == 'eBird')] = as.numeric(sapply(eventList$eventDate[which(eventList$datasetKey == 'eBird')], function(x) unlist(str_split(x," "))[2]))

summary(eventList$Year); print(paste("Should be FALSE: ", anyNA(eventList$Year)))
dim(eventList)

## changed the bioregion term so that it is one per polygon 
saveRDS(eventList, "outputs/EventList/EventListMerged/EventList.GBIF.eb.effort80.birds.Year_2023.RDS")

eventList = readRDS("outputs/EventList/EventListMerged/EventList.GBIF.eb.effort80.birds.Year_2023.RDS")

## remove events that are too long ago
nrow(eventList)
eventList = unique(eventList);nrow(eventList)
eventList = subset(eventList, Year >= youngestDate)
dim(eventList)

############################################################################################
## 3 - number of visits per grid cell (= survey effort) -----

## survey effort per guild

setDT(eventList)
effort = eventList[, .(effort = uniqueN(.SD)), by = df_grid, .SDcols = c("eventDate", "datasetKey")]
summary(effort)

eventList$effort_Aves = effort$effort[match(eventList$df_grid,effort$df_grid)]
summary(eventList$effort_Aves)



####################################################################################################################
## 4 - keep grid cells with larger survey effort than threshold (=50) and at least 15 protected ------
## AND 15 non protected cells and rarefy

eventList = subset(eventList, effort_Aves >= 50); dim(eventList)

# 20/01/2023 over5yo : only protected areas that were created before 2015
# 23/01/2023 after desig : filtered eventList for events only recorded AFTER the designation of the PA and didn't filter ebird by 2km dist
# saveRDS(eventList, "outputs/EventList/EventListMerged/EventList.GBIF.eb.effort80.NeighbouringPAID.new.BIRDS.Over5yo.AfterDesig.NewEffort.RDS")
# saveRDS(eventList, "outputs/EventList/EventListMerged/EventList.GBIF.eb.effort60.NeighbouringPAID.new.BIRDS.Over5yo.AfterDesig.NewEffort.RDS")
saveRDS(eventList, "outputs/EventList/EventListMerged/EventList.GBIF.eb.effort50.NeighbouringPAID.new.BIRDS.Over5yo.AfterDesig.NewEffort.RDS")



eventList = readRDS("outputs/EventList/EventListMerged/EventList.GBIF.eb.effort50.NeighbouringPAID.new.BIRDS.Over5yo.AfterDesig.NewEffort.RDS")

merged.PAs = readRDS("outputs/WDPA.Grid/GridCells.NeighbouringPAID.new.Over5yo.RDS")
merged.PAs$ID.new.Bioregion = apply(merged.PAs[,c('NeighbouringPAID.new','Bioregion')],1, FUN = function(x) paste(x, collapse = '_'))

## number of protected + non protected grid cells per site
# eventList$Bioregion = merged.PAs$Bioregion[match(eventList$df_grid,merged.PAs$Value)]
# eventList_breedingGBIF = subset(eventList, (datasetKey != "eBird"))
# eventList_breedingGBIF = subset(eventList_breedingGBIF, month(as.POSIXct(eventDate)) %in% c(4,5,6,7,8))

eventList = subset(eventList, datasetKey == "eBird")
eventList = rbind(eventList, eventList_breedingGBIF)

effort.min = ddply(eventList, c('NeighbouringPAID.new','Bioregion'), summarise, 
                   nb.inside = length(unique(df_grid[which(!is.na(ID.new))])), 
                   nb.outside = length(unique(df_grid[which(is.na(ID.new))])), 
                   effort.min = min(effort_Aves)) ## to rarefy: minimum survey effort > 80 per site

setDT(eventList)
setDT(effort.min)
eventList = merge(eventList, effort.min, by = c('NeighbouringPAID.new','Bioregion'))
dim(eventList)
eventList.min = subset(eventList, nb.outside >= 15 & nb.inside >= 15); dim(eventList.min)
eventList.min$eventID = seq(1:nrow(eventList.min))
eventList.min$ID.new.Bioregion = apply(eventList.min[,c('NeighbouringPAID.new','Bioregion')],1, 
                                       FUN = function(x) paste(x, collapse = '_'))

effort.min$ID.new.Bioregion = apply(effort.min[,c('NeighbouringPAID.new','Bioregion')],1, 
                                    FUN = function(x) paste(x, collapse = '_'))

length(unique(eventList.min$ID.new.Bioregion)) 
summary(eventList.min$effort.min)

### Rarefying - for each site we rarefy the number of event for each local community to the minimum of that site
### We do that 100 times to assess how the rarefying could affect the results

count = 0
while(count<100){
  count = count + 1
  
  # sample effort.min number of events among the nb events from eventList min for each couple of PAID and Bioregion
  sample = rbindlist(apply(unique(eventList.min[,c('df_grid','effort.min')]), 1, 
                            FUN = function(x){ cat(as.numeric(x['df_grid'])) ; 
                              sub = subset(eventList.min, df_grid == as.numeric(x['df_grid']))
                              sample = sub[sample(1:nrow(sub), x['effort.min']),]
                              return(sample)}))
  print(paste('Should be true : ',nrow(sample) == 
                sum(effort.min$effort.min[which(effort.min$ID.new.Bioregion %in% eventList.min$ID.new.Bioregion)]*
                      ddply(eventList.min, 'ID.new.Bioregion', summarise, nb = length(unique(df_grid)))$nb)))
  
  saveRDS(sample, paste0("outputs/EventList/EventListMerged/Rarefying/EventListRarefied.Over5yo.NewEffort.", count ,'.rds'))

}

summary(eventList.min.sub$nb.inside);summary(eventList.min.sub$nb.outside); summary(eventList.min.sub$effort.min)
length(unique(eventList.min.sub$ID.new.Bioregion)) 

setDT(sample)
sample[, .(effort = uniqueN(.SD)), by = df_grid, .SDcols = c("eventDate", "datasetKey")] # should be equal to effort min



#################################################
## 4 - subset the GBIF and ebird files using eventList.min to get rarefied dataset -----

eventList.min.sub = readRDS('outputs/EventList/EventListMerged/Rarefying/EventListRarefied.Over5yo.NewEffort.1.rds')


# key = all.south.keys[3]
# cardinal = "South"

################################################################################
## TAXONOMIC HOMOGENEISATION -----
## 1. Get the list of all the species in western EU from all the GBIF occurrence files
df = lapply(list.files('outputs/dta.Occurrence/GBIF/Filtered.effort/'), 
            function(x) unique(fread(paste0('outputs/dta.Occurrence/GBIF/Filtered.effort/',x), 
                                     select = c('species','order','class','taxonKey'))))
speciesList = unique(rbindlist(df))
write.csv(speciesList, 'outputs/SpeciesList_eBird_GBIF_FilteredEffort.csv')


speciesList.birds = subset(speciesList, class == "Aves")
m_info.birds = subset(m_info, Class == "Aves")

################################################################################
## 2. Merge the gbif species with the metaweb names

library(rgbif)
## use rgbif function to match as many metaweb species to the GBIF taxonomy
gbif_names = name_backbone_checklist(m_info.birds$Species) 
gbif_names = subset(gbif_names, class == "Aves")
# verbatim_name are the queried names (i.e. master names here)


## 1. match master codes to each gbif species
gbif_names = merge(gbif_names, m_info[,c('Species','Species.ID')], by.x = 'verbatim_name', by.y = 'Species')
names(gbif_names)[which(names(gbif_names)%in%c('Species.ID','verbatim_name'))] = c('masterSpecies','masterCode')

anyNA(gbif_names$masterCode) ## we have all species.
nrow(speciesList.birds) - nrow(gbif_names) # we ignore 247 species who are not in the EU metaweb

write.csv(gbif_names, 'outputs/Species_Taxonomy_Metaweb_birds.csv')
# 509 species

#######################################################################################################

speciesList = read.csv('outputs/Species_Taxonomy_Metaweb_birds.csv')
# speciesList = read.csv('E:/TheseSwansea/dta.Occurrence/GBIF/Europe/MetawebSpecies/GBIFSpecies_BackboneTaxonomy_minfo.withEast.csv')

library(stringr)
  
## keep only 'events' that are within PAs with enough grid cells (15 outside + 15 inside) and sufficient survey effort
sf_df.sub = rbindlist(lapply(paste0("outputs/dta.Occurrence/GBIF/Filtered.effort/",
                                     str_subset(list.files("outputs/dta.Occurrence/GBIF/Filtered.effort"),'csv')), function(x) {
    sf_df = fread(x, colClasses = list(character = 'eventDate'))
    dim.initial = nrow(sf_df)
    
    sf_df = subset(sf_df, class == 'Aves' & 
                     eventDate %in% eventList.min.sub$eventDate & 
                     datasetKey %in% eventList.min.sub$datasetKey & 
                     df_grid %in% eventList.min.sub$df_grid)
    sf_df$masterCode = speciesList[match(sf_df$taxonKey, speciesList$usageKey),'masterCode']
    print(paste('We keep', nrow(sf_df)/dim.initial, '% of the observations'))
    return(sf_df)}), fill=TRUE)


print(anyNA(sf_df.sub$masterCode))
  
saveRDS(sf_df.sub, paste0("outputs/dta.Occurrence/Rarefied/RarefiedOccData.Over5yo_gbif.NewEffort1.rds"))


## do the same with the eb data

files = paste0("data/eBird/Europe/Exported/0.Merged/chlist/eBird.Grid/",
              list.files("data/eBird/Europe/Exported/0.Merged/chlist/eBird.Grid/"))
ind = 0
for (f in files){
  ind = ind + 1
  eb = readRDS(f)
  print(dim(eb))
  
  if(!'date' %in% names(eb)){
    eb$date = paste(eb$day, eb$year)
  }
  
  eb.sub = subset(eb, date %in% eventList.min.sub$eventDate[which(eventList.min.sub$datasetKey == 'eBird')] & 
                    Value %in% eventList.min.sub$df_grid[which(eventList.min.sub$datasetKey == 'eBird')])
  
  print(paste('We keep', nrow(eb.sub)/nrow(eb), '% of the observations for file', f))
  saveRDS(eb.sub, paste0("outputs/dta.Occurrence/Rarefied/eBird.RarefiedOccData.Over5yo.NewEffort.", ind,'.rds'))
  
}


#################################################
## 5 - compute network metrics for each grid cell -----

library(igraph)
library(plyr)
metaweb_europe <- read.graph('data/metaweb-europe.graphml', format = 'graphml')
m <- as_adjacency_matrix(metaweb_europe, attr = 'copresence', sparse=FALSE)
m[is.nan(m)] <- 0
source("data/Galiana2021_network-area-europe-master/whois-function.R")

m_code = m # transform species names into species codes from TETRA EU
rownames(m_code) <- unlist((apply(as.matrix(rownames(m_code)), 1, function(x){ n <- whois(SPPNAME=x)[[1]]; if(length(n) == 0) {return (x)} else{return (n)} }  ) ))
colnames(m_code) <- unlist((apply(as.matrix(colnames(m_code)), 1, function(x){ n <- whois(SPPNAME=x)[[1]]; if(length(n) == 0) {return (x)} else{return (n)} }  ) ))

## source code for localMetrics (FUNCTION THAT COMPUTES NETWORK METRICS)
source('Script/MetricCalculationFunction.R')

eventList.min.sub = readRDS('outputs/EventList/EventListMerged/Rarefying/EventListRarefied.Over5yo.NewEffort.1.rds')


## we now have to (1) get the occurrence data corresponding to each event in eventList.min
## and (2) combine into a dataframe and (3) compute network metrics

eventList.min.sub = subset(eventList.min.sub, class == 'Aves')
species.dt = data.frame()
nb = c()
for(f in str_subset(list.files('outputs/dta.OCcurrence/Rarefied/'), 'rds')){
  df = readRDS(paste0('outputs/dta.OCcurrence/Rarefied/',f))
  names(df)[which(names(df) == 'Value')] = 'df_grid'
  names(df)[which(names(df) == 'date')] = 'eventDate'
  if(!'datasetKey' %in% names(df)){
    df$datasetKey = 'eBird'
  }
  
  setDT(df)
  setDT(eventList.min.sub)
  dim.in = nrow(df)
  df1 = df
  df = merge(df[,c('masterCode','df_grid','eventDate','datasetKey')], 
             eventList.min.sub[,c('df_grid','eventDate','datasetKey','eventID')], 
             by = c('eventDate','datasetKey','df_grid'))
  df = df[-which(duplicated(df)),]
  
  nb = rbind(nb, unique(df[,c('eventDate','datasetKey','df_grid','eventID')]))
  species.dt = rbind(species.dt,ddply(df, 'df_grid', summarise, species = unique(masterCode)))
}

# only birds and no ornitho
paste('Should be TRUE : ', nrow(unique(nb)) == nrow(eventList.min.sub[which(eventList.min.sub$datasetKey != 'ornitho.ch')]))
species.dt$pres = 1

# regroup everything into a list
if(anyDuplicated(species.dt[,c("df_grid","species")])>0){
  species.dt = species.dt[-which(duplicated(species.dt[,c("df_grid","species")])),] # remove duplicates from different files
}

print(summary(ddply(species.dt, 'df_grid', summarise, nb = length(species))))
print(summary(ddply(nb, c('df_grid'), summarise, nb = length(eventID))))
print(summary(ddply(eventList.min.sub, c('df_grid'), summarise, nb = length(eventID))))
summary(ddply(eventList.min.sub, c('Bioregion',"NeighbouringPAID.new"), summarise, 
                    nb.in = length(unique(df_grid[which(!is.na(ID.new))])),
                    nb.out = length(unique(df_grid[which(is.na(ID.new))]))
                    ))
setdiff(eventList.min.sub$df_grid, nb$df_grid)
setdiff(eventList.min.sub$df_grid, species.dt$df_grid)

## pivot into grid
grid = tidyr::pivot_wider(species.dt[,c("df_grid","species","pres")], names_from = species, values_from = pres, values_fill = 0)
paste0('should be TRUE : ', nrow(grid) == length(unique(species.dt$df_grid)))
grid = merge(grid, cells_info[,c('x','y','PageName','Value')], by.x = 'df_grid', by.y = 'Value', all.x = T)
dim(grid)

codes = colnames(grid)[which(!colnames(grid)%in%c('x','y','Value','df_grid','X','PageName','country','Count'))]
hist(rowSums(grid[,codes])) ## number of species per grid cell


# check
print(ggplot() + geom_sf(data = europe) + geom_point(data = grid, aes(x = x, y = y, colour = B207), size = .5) +
        xlim(min(grid$x, na.rm = T), max(grid$x, na.rm = T)) + ylim(min(grid$y, na.rm = T), max(grid$y, na.rm = T))
)


### getting body mass data
temp <- tempfile()
download.file("https://ndownloader.figshare.com/files/5631081", temp)
ET <- read.table(temp, header = TRUE, fill  = TRUE, quote = "\"", stringsAsFactors = FALSE,sep = "\t")
unlink(temp)

codes_df = data.frame(masterCode = unlist(lapply(codes,function(x) str_replace(whois(x), "_", " "))))
codes_df$ET_names = ET$Scientific[match(codes_df$masterCode,ET$Scientific)]

missing = codes_df$masterCode[which(is.na(codes_df$ET_names))]

library(rgbif) ## for getting species names from taxonomy
library(readxl)
TaxMatch <- read_excel("data/TaxonomyHomogeneisationEltonTraits.xlsx")
names = name_backbone_checklist(missing)
names = subset(names, species %in% ET$Scientific)
codes_df$ET_names[match(names$verbatim_name, codes_df$masterCode)] = names$species
codes_df$ET_names[match(TaxMatch$masterSpecies,codes_df$masterCode)] = TaxMatch$EltonTraits

bodySize = merge(ET,codes_df, by.x = "Scientific", by.y = "ET_names")
hist(bodySize$BodyMass.Value)


options(cheddarMaxQueue = 0)

## compute network metrics for each grid cell
res = localMetrics(master_reduced = grid, iterations = nrow(grid), m_code = m_code, bodySize = bodySize)

res_traits = get_species_traits(master_reduced = grid, iterations = nrow(grid), m_code = m_code, bodySize = bodySize)
saveRDS(res_traits, paste0("outputs/1.Metrics/Species_traits_grid_4662cells.rds"))


anyDuplicated(res)
anyDuplicated(res$Value)
setdiff(res$PageName, grid$PageName)
summary(res)

if(length(which(rowSums(grid[,codes]) < 2)) > 0){
  print(paste0('Should be true : ', nrow(res) == nrow(grid[-which(rowSums(grid[,codes]) < 2),])))
}  else{
  print(paste0('Should be true : ', nrow(res) == nrow(grid)))
}

res = merge(res, cells_info, by = 'PageName', all.x = T)
res = merge(res, unique(eventList.min.sub[,c('df_grid','NeighbouringPAID.new','Bioregion','effort.min','nb.inside','nb.outside','dist.PA.min','effort_Aves')]), by.x = 'Value', by.y = 'df_grid', all.x = T)

## save
## re ran with new taxonomy 
saveRDS(res, paste0("outputs/1.Metrics/NetworkMetrics.Rarefied.nb15.effort50.Over5yo.AfterDesig.NewEffort.BodyMass.rds"))
saveRDS(species.dt, paste0("outputs/1.Metrics/Sample/OccurrenceDataSample.Rarefied.nb15.effort50.Over5yo.AfterDesig.NewEffort.rds"))


# #### TEMP FOR MISSING CELLS AFTER CHIANGING BIOREGION
# # SHOULD RERUN ALL
# res = readRDS(paste0("E:/TheseSwansea/dta.Occurrence/GBIF/Europe/MetawebSpecies/Over5yo/Metrics/NetworkMetrics.Rarefied.nb15.effort80.Over5yo.AfterDesig.23012023.rds"))
# res1 = readRDS(paste0("E:/TheseSwansea/dta.Occurrence/GBIF/Europe/MetawebSpecies/Over5yo/Metrics/NetworkMetrics.Rarefied.nb15.effort80.Over5yo.AfterDesig.20022023.missingCells.rds"))
# 
# res = res[-which(res$Value %in% extra.grids),]
# res = rbind(res[names(res1)],res1)
# 
# saveRDS(res, paste0("E:/TheseSwansea/dta.Occurrence/GBIF/Europe/MetawebSpecies/Over5yo/Metrics/NetworkMetrics.Rarefied.nb15.effort80.Over5yo.AfterDesig.20022023.MissingCells.Combined.rds"))


################################################################################
################################################################################
### plotting the networks 
  cell.protected = which(grid$df_grid == '104125') # 68 mediterranean
  cell.nonprotected = which(grid$df_grid == '103752')
  
  pixel.p = grid[cell.protected,] ## local community (sp that occur within this 10x10km grid cell)
  pixel.np = grid[cell.nonprotected,] ## local community (sp that occur within this 10x10km grid cell)
  
  sppSTRING.p <- names(pixel.p)[which(pixel.p==1 & str_detect(names(pixel.p), pattern = "[:digit:]"), arr.ind=FALSE)]
  sppSTRING.np <- names(pixel.np)[which(pixel.np==1 & str_detect(names(pixel.np), pattern = "[:digit:]"), arr.ind=FALSE)]
  local_web.p = m_code[sppSTRING.p,sppSTRING.p]
  local_web.np = m_code[sppSTRING.np,sppSTRING.np]
  
  colnames(local_web.p) = rownames(local_web.p) = unlist(lapply(colnames(local_web.p), whois))
  colnames(local_web.np) = rownames(local_web.np) = unlist(lapply(colnames(local_web.np), whois))

  # creating a cheddar community to look for isolated nodes
  community.p = cheddar::Community(nodes=data.frame(node=colnames(local_web.p)), trophic.links=PredationMatrixToLinks(local_web.p), properties=list(title='Community'));
  community.np = cheddar::Community(nodes=data.frame(node=colnames(local_web.np)), trophic.links=PredationMatrixToLinks(local_web.np), properties=list(title='Community'));
  
  # new local web that only keeps connected nodes
  to_keep.p = setdiff(colnames(local_web.p),IsolatedNodes(community.p)) # species that are not 'isolated' i.e. have at least 1 prey or 1 pred
  to_keep.np = setdiff(colnames(local_web.np),IsolatedNodes(community.np)) # species that are not 'isolated' i.e. have at least 1 prey or 1 pred
  
  local_web.p = local_web.p[to_keep.p,to_keep.p]
  local_web.np = local_web.np[to_keep.np,to_keep.np]
  
  community.p = cheddar::Community(nodes=data.frame(node=colnames(local_web.p)), trophic.links=PredationMatrixToLinks(local_web.p), properties=list(title='Community'));
  community.np = cheddar::Community(nodes=data.frame(node=colnames(local_web.np)), trophic.links=PredationMatrixToLinks(local_web.np), properties=list(title='Community'));
  
  plot(community.p, show.nodes.as='labels', node.labels='node')
  plot(community.np, show.nodes.as='labels', node.labels='node')
  
  par(mfrow = c(1,2))
  plot(community.p, main = 'protected network - 68 Mediterranean bioregion')
  plot(community.np, main = 'non protected network - 68 Mediterranean bioregion')
  