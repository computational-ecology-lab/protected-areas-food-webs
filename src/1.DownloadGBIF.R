###############################################################################
### Code for Protected areas with clear management goals enhance avian food webs
## Downloading GBIF data 
###############################################################################

setwd('E:/TheseSwansea/WDPA')

library(rgbif)
library(rnaturalearth)
library(ggplot2)

europe = ne_countries(continent = "Europe", returnclass = "sf", scale = 'large')


m_info = read.csv("data/Galiana2021_network-area-europe-master/species_codes_and_taxonomy.csv", sep = ";")
gbif_names = name_backbone_checklist(m_info$Species) 
gbif_taxon_keys <- gbif_names$usageKey

## set up credentials https://docs.ropensci.org/rgbif/articles/gbif_credentials.html

## North 
eu_poly1 = 'POLYGON((-11.42578 51.1946,3.00000 51.1946,3.00000 73.13515,-11.42578 73.13515,-11.42578 51.1946))'
eu_poly2 = "POLYGON((15.0000 51.1946,29.03027 51.1946,29.03027 73.13515,15.0000 73.13515,15.0000 51.1946))" ##
eu_poly3.1 = 'POLYGON((3.00000 51.1946,9.0000 51.1946,9.0000 73.13515,3.00000 73.13515,3.00000 51.1946))' 
eu_poly3.2 = 'POLYGON((9.00000 51.1946,15.0000 51.1946,15.0000 73.13515,9.00000 73.13515,9.00000 51.1946))' 
eu_poly4 = 'POLYGON((29.03027 51.1946,41 51.1946,41 73.13515,29.03027 73.13515,29.03027 51.1946))'
## South 
eu_poly5 = 'POLYGON((-11.42578 30.17552,3.00000 30.17552,3.00000 51.1946,-11.42578 51.1946,-11.42578 30.17552))'
eu_poly6 = 'POLYGON((3.00000 30.17552,15.0000 30.17552,15.0000 51.1946,3.00000 51.1946,3.00000 30.17552))'
eu_poly7 = 'POLYGON((15.0000 30.17552,29.03027 30.17552,29.03027 51.1946,15.0000 51.1946,15.0000 30.17552))'
eu_poly8 = 'POLYGON((29.03027 42,41 42,41 51.1946,29.03027 51.1946,29.03027 42))'


usethis::edit_r_environ()
# GBIF_USER="jwaller"
# GBIF_PWD="safe_fake_password_123"
# GBIF_EMAIL="jwaller@gbif.org"


### Download and put into North or South folder
res <- occ_download(type = 'and', 
                      #pred("country", code),
                      pred_within(eu_poly8),
                      pred_in("taxonKey", gbif_taxon_keys), # only sp in the metaweb
                      pred_in("basisOfRecord", c('HUMAN_OBSERVATION')),
                      pred("hasCoordinate", TRUE),
                      pred("hasGeospatialIssue", FALSE),
                      ### removes the issues:
                      # Zero coordinate; Country coordinate mismatch; Coordinate invalid; Coordinate out of range
                      pred_gte("year", 2000), ## greater than 2000
                      #pred_not(pred('collectionCode','EBIRD')), ## dowloaded ebird dta seperately
                      #pred_not(pred('institutionCode','EOD - eBird Observation Dataset')), 
                      #pred_not(pred("establishmentMeans","MANAGED")), ## species in zoo
                      pred_lte("coordinateUncertaintyInMeters",10000),
                      pred_notnull("coordinateUncertaintyInMeters"),
                      format = "SIMPLE_CSV")

## 0222868-210914110416597 ## poly1                  
## 0222874-210914110416597 ## poly5                    
## 0222877-210914110416597 ## poly6
## 0222905-210914110416597 ## poly3.1
## 0222907-210914110416597 ## poly3.2
## 0222933-210914110416597 ## poly7
## 0282245-210914110416597 ## poly8

### checking that the boxes are somewhat equal
NW = data.frame(x=c(-11.42578, 3.00000, 3.00000, -11.42578, -11.42578),
                y=c(51.1946, 51.1946, 73.13515, 73.13515, 51.1946))
NMW = data.frame(x=c(3.00000, 9.0000, 9.0000, 3.00000, 3.00000), 
                 y=c(51.1946, 51.1946, 73.13515, 73.13515, 51.1946))
NME = data.frame(x=c(9.00000, 15.0000, 15.0000, 9.00000, 3.00000), 
                 y=c(51.1946, 51.1946, 73.13515, 73.13515, 51.1946))
NE = data.frame(x=c(15.0000, 29.03027, 29.03027, 15.0000, 15.0000), ##
                y=c(51.1946, 51.1946, 73.13515, 73.13515, 51.1946)) 
NE1 = data.frame(x=c(29.03027, 41, 41, 29.03027, 29.03027), 
                y=c(51.1946, 51.1946, 73.13515, 73.13515, 51.1946)) 

SW = data.frame(x=c(-11.42578, 3.00000, 3.00000, -11.42578, -11.42578),
                y=c(30.17552, 30.17552, 51.1946, 51.1946, 30.17552))
SM = data.frame(x=c(3.00000, 15.0000, 15.0000, 3.00000, 3.00000),
                y=c(30.17552, 30.17552, 51.1946, 51.1946, 30.17552))
SE = data.frame(x=c(15.0000, 29.03027, 29.03027, 15.0000, 15.0000),
                y=c(30.17552, 30.17552, 51.1946, 51.1946, 30.17552))
SE1 = data.frame(x=c(29.03027, 41, 41, 29.03027, 29.03027),
                y=c(30.17552, 30.17552, 51.1946, 51.1946, 30.17552))

ggplot(data = europe) + geom_sf() + geom_point(data = NW, aes(x = x, y = y), size = 5, shape = 3) +
  geom_point(data = NE, aes(x = x, y = y), size = 5, shape = 3) +
  geom_point(data = NMW, aes(x = x, y = y), size = 5, shape = 3) +
  geom_point(data = NME, aes(x = x, y = y), size = 5, shape = 3) +
  geom_point(data = SW, aes(x = x, y = y), size = 5, shape = 3) + 
  geom_point(data = SE, aes(x = x, y = y), size = 5, shape = 3) +
  geom_point(data = SE1, aes(x = x, y = y), size = 5, shape = 3) + 
  geom_point(data = NE1, aes(x = x, y = y), size = 5, shape = 3) +
  xlim(c(-15,45)) + ylim(c(25, 80))    

keys = c("0222825-210914110416597", "0222868-210914110416597", "0222905-210914110416597", "0222907-210914110416597", "0282225-210914110416597",
        "0222874-210914110416597", "0222877-210914110416597", "0222933-210914110416597", "0282245-210914110416597")

# get dataset info
library(rgbif)
library(plyr)
for (key in keys[1:9]){
  
  index = 0 ;     
  dta = data.frame()
  while(nrow(dta) %% 1000 == 0){
    dta = rbind.fill(dta,occ_download_datasets(key, start = 1000*index, limit = 1000)$results)
    index = index + 1
    print(nrow(dta))
    print(index)
  }
  write.csv(dta,paste0("data/outputs/dta.Occurrence/GBIF/datasets.",key,".csv"))
}

library(data.table)
setwd("data/outputs/dta.Occurrence/GBIF/")
Datasets = rbindlist(lapply(list.files(), read.csv), fill = T)

Datasets = ddply(Datasets, "datasetTitle", summarise, numberRecords = sum(numberRecords))


################################################################################
## Subset each polygon and save into 'Split' folder
################################################################################

## columns names to keep:
names = c("gbifID","occurrenceID",'species','class','order',"scientificName",'datasetKey',"countryCode","occurrenceStatus","individualCount",
          "decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters","eventDate","year" ,"taxonKey","taxonRank",
          "speciesKey","basisOfRecord" , "institutionCode" ,"collectionCode","recordNumber","identifiedBy","issue")

key = '0222905-210914110416597';cardinal = "North" ## SPLIT
key = '0222907-210914110416597';cardinal = "North" ## SPLIT
key = '0282225-210914110416597';cardinal = "North" 
key = '0222868-210914110416597';cardinal = "North" ## SPLIT

key = '0222874-210914110416597';cardinal = "South" ## SPLIT
key = '0222877-210914110416597';cardinal = "South" ## SPLIT
key = '0222933-210914110416597';cardinal = "South"
key = '0282245-210914110416597';cardinal = "South"


# splitting the data and saving
setwd(paste0("data/outputs/dta.Occurrence/GBIF/",cardinal,'/',key))

skip = 0
nb = 3e6 # number of rows per individual splitted file
rows = nb
i = 0
names_index = match(names,names(fread(paste0(key,'.csv'), nrows = 1, header = T))) # need to provide row numbers instead of col names to fread

while (rows == nb){
  df = fread(paste0(key,'.csv'), skip = skip, nrows = nb , select = names_index, nThread = 10, header = F)
  names(df) = names
  write.csv(df, paste0('Split/',key,'.',i,'.csv'))
  
  i = i + 1
  skip = i*nb
  rows = nrow(df)
}





                    