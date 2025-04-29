################################################################################
### Code for Protected areas with clear management goals enhance avian food webs
### Format raw ebird data into eventList (one event = one data collection protocol + one date/time)



setwd("F:/TheseSwansea/WDPA/")
version="May-2022"

library(auk) ; 
library(plyr) ; 
library(readr) ; 
library(nlme) ; 
library(data.table) ; 
library(stringr)
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_)) 

##################
## raw ebird files

files.countries<-list.files("data/eBird",recursive = T)
files.countries = str_subset(files.countries, paste0(version,'.txt'))

################################################################################
## Split data for Spain and GB as very large files

files.countries.ES.GB = str_subset(files.countries, 'ES|GB')

## save colnames 
names = names(data.table::fread(files.countries.ES.GB[1], nrows = 1, nThread = 10, header = T))

for (key in files.countries.ES.GB){
  
  skip = 0
  nb = 3e6 # number of rows per individual splitted file
  rows = nb
  i = 0
  
  while (rows == nb){
    
    df = data.table::fread(key, skip = skip, nrows = nb, nThread = 10, header = F)
    names(df) = names
    
    write.table(df, paste0('data/eBird/Europe/Split/',str_split(key, '_')[[1]][2],'.',i,'.csv'))
    
    i = i + 1
    skip = i*nb
    rows = nrow(df)
  }
}

## update files.countries to include the split files for GB and ES instead of the full ones:
files.countries = c(files.countries[which(!files.countries %in% str_subset(files.countries, 'ES|GB'))],
                    paste0('data/eBird/Europe/Split/',list.files("data/eBird/Europe/")))



for (i in c(1:length(files.countries))){
  
  eb<-read_ebd(files.countries[i], rollup=T, unique=T) 
  
  
  ### Subset data
  # Remove useless columns
  eb[,c("global_unique_identifier", "last_edited_date", "subspecies_common_name", "breeding_bird_atlas_category", "effort_area_ha", "county", "subnational2_code", "iba_code", "bcr_code", "usfws_code", "atlas_block", "locality", "locality_id", "locality_type", "first_name", "last_name", "has_media", "reviewed", "x")]<-NULL
  # Remove domestic species (Columba livia only)
  eb<-subset(eb, eb$category != "domestic")
  # Remove checklists that did not report all species
  eb<-subset(eb, eb$all_species_reported==T)
  # Remove disapproved observations
  eb<-subset(eb, eb$approved==T)
  # Select protocols
  eb<-subset(eb, eb$protocol_type %in% c("Traveling","Stationary","Historical"))
  # Remove data with no duration
  eb<-subset(eb, is.na(eb$duration_minutes)==FALSE)
  # Remove historical protocols with no distance reported (ie only stationary are allowed not to have a distance reported)
  eb[which(eb$protocol_type == 'Stationary'),'effort_distance_km'] = 0
  
  
  ### Format date
  eb$year<-as.numeric(format(eb$observation_date, "%Y"))
  eb$month<-as.numeric(format(eb$observation_date, "%m"))
  eb$day<-as.numeric(format(eb$observation_date, "%j"))
  
  # Remove short or long checklists
  eb<-subset(eb, eb$duration_minutes>30 & eb$duration_minutes<600)
  # Remove old data
  eb<-subset(eb, eb$year>=2005)
  # Remove long travelling distances
  eb<-subset(eb, effort_distance_km<=5)
  eb<-subset(eb, !is.na(effort_distance_km))
  
  if( i <= 45 ){
    saveRDS(eb, file= paste0("Exported/", 
                             paste0(str_split(files.countries[i],'_|\\.', simplify = T)[c(2,5)],collapse = '.'), ".rds", sep=""))
  } else {
    saveRDS(eb, file= paste0("Exported/", 
                             paste0(str_split(files.countries[i],'/|\\.', simplify = T)[c(2,3)],collapse = '.'), ".rds", sep=""))
  }
  cat(i)
  
}


saveRDS(eb, file=paste("data/eBird/Europe/Exported/",version,".Export.", substr(list.files()[i],5,6), ".rds", sep=""))



################################################################################
### 2. Get unique species' names for taxonomy matching and combine the files above into three main rds files
################################################################################

setwd("E:/TheseSwansea/dta.Occurrence/e-birds/ebirdRAWdata/Europe/Exported/")

files.countries.merged = paste("data/eBird/Europe/Exported",str_subset(list.files("data/eBird/Europe/Exported"),'.rds'))
seq = c(0,35,58,88)
for (i in 1:(length(seq)-1)){
  
  df<- rbindlist(lapply(files.countries.merged[(seq[i]+1):seq[i+1]], readRDS), fill = T) 
  print(unique(df$country))
  
  saveRDS(df, paste0("data/eBird/Europe/Exported/0.Merged/Merged.data.", version,'.Europe.',i, ".rds"))
  
  species = unique(df[,c('scientific_name','common_name','taxon_concept_id')])
  write.csv(species, paste0('data/eBird/Europe/Exported/Taxonomy/speciesList.', version,'.Europe.',i, '.csv'))
  
}

################################################################################
### TAXONOMIC MATCHING ebird with master data


library(auk)

m_info = read.csv("data/Galiana2021_network-area-europe-master/species_codes_and_taxonomy.csv", sep = ";")


################################################################################
## 1. Merge the ebird species the metaweb sp names

ebirdSpList = read.csv('data/eBird/Europe/Exported/0.Merged/Taxonomy/speciesList.May-2022.csv')

## extract the equivalent names of the master species in the ebird taxonomy (auk package, ebird_species function)
ebMaster_names = data.frame(m_info[which(m_info$Class == 'Aves'),c('Species','Species.ID')],
                            ebSpecies = ebird_species(m_info$Species[which(m_info$Class == 'Aves')]))

## some species here are from other continents / are probably mistakes or exceptions so 
## would not be expected to be included in the metaweb
## check with the birdlife international dist BOTW 2021
# write.csv(missing, 'missingSp.ebirdFeb-2022.csv')
## run on qgis expression selector of BOTW maps: array_contains(aggregate('missingSp.ebirdFeb-2022', 'array_agg', "x"), "sci_name") 
## selects only the dist maps of missing species that occur in europe (just frew a rectangle over europe)

missingINeurope = read.csv('data/eBird/Europe/Exported/0.Merged/Taxonomy/missingSp.ebirdFeb-2022.Europe.csv') ## species missing that actually have an AOO that overlaps with europe
missingINBOTW = read.csv('data/eBird/Europe/Exported/0.Merged/Taxonomy/missingSp.ebirdFeb-2022.InBOTW.csv') ## all species did not have BOTW maps

outsideEU = setdiff(missingINBOTW$sci_name, missingINeurope$sci_name)
intersect(outsideEU, ebMaster_names$Species) ; intersect(outsideEU, ebMaster_names$ebSpecies)
ebirdSpList = ebirdSpList[-which(ebirdSpList$scientific_name %in% outsideEU),]

ebMaster_names[which(ebMaster_names$ebSpecies %in% setdiff(ebMaster_names$ebSpecies,ebirdSpList$scientific_name)),'ebSpecies'] = NA
names(ebMaster_names) = c('masterSpecies','masterCode','ebSpecies')

missing_auk = ebMaster_names[which(is.na(ebMaster_names$ebSpecies)),] # species in the master data for whom we have no eb equivalent
# 94 species unmatched

missing = setdiff(ebirdSpList$scientific_name, ebMaster_names$ebSpecies) # species from the ebird european occurrence dataset for whom we have not master equivalent
# 106 species unmatched

## for the remaining species we need to find whether they have synonym names that 
# would match a species in the metaweb

## create a column to remember which species names we inferred
ebMaster_names$inferred = NA

## Tring different methods :
################################################################################
## library for taxonomic queries/synonyms
library(taxize)
library(data.table)

# https://github.com/hannahwauchope/PAImpact/blob/master/1_CollateCleanRawData.R
#The aim here is to get all possible synonyms for each species
#Get alternative species names using taxize

### 1. TOL
## a. Species from ebird
TreeofLifeCleaning <- function(TreeofLifeID, SpeciesList, SpeciesToMatch){
  #Pull out those with only one match, get all synonyms
  TreeofLifeIDOneMatch <- TreeofLifeID[lapply(TreeofLifeID, nrow)==1]
  TreeofLifeIDOneMatch <- rbindlist(lapply(1:length(TreeofLifeIDOneMatch), function(x){
    TOL <- TreeofLifeIDOneMatch[[x]]
    TOL$CountsSpecies <- names(TreeofLifeIDOneMatch)[x]
    names(TOL) <- c("UniqueName", "MatchedName", "IsApproxMatch", "IsSynonym", "NomenclatureCode","Score", "Flags",'IsSuppressedFromSynth', "OTTID", "Rank",'Source', "CountsSpecies")
    #str_split(names(TreeofLifeID[[1]]),'\\.',simplify = T)[,3]
    return(TOL)
  }), fill = T)
  
  TreeofLifeIDOneMatch <- rbindlist(lapply(unique(TreeofLifeIDOneMatch$CountsSpecies), function(x){
    ok <- subset(TreeofLifeIDOneMatch, CountsSpecies == x)
    Names <- as.data.frame(t(as.data.frame(c(unique(c(ok$CountsSpecies, ok$UniqueName, ok$MatchedName))))))
    row.names(Names) <- NULL
    if(ncol(Names)==1){
      names(Names) <- "Species"
    } else {
      names(Names) <- c("Species", paste0("Synonym_", 1:(ncol(Names)-1)))
    }
    return(Names)
  }), fill=TRUE)
  
  # Pull out those with many matches
  TreeofLifeIDManyMatch <- TreeofLifeID[lapply(TreeofLifeID, nrow)>1]
  
  if(length(TreeofLifeIDManyMatch)>0){
    TreeofLifeIDManyMatch <- rbindlist(lapply(1:length(TreeofLifeIDManyMatch), function(x){
      TOL <- TreeofLifeIDManyMatch[[x]]
      TOL$CountsSpecies <- names(TreeofLifeIDManyMatch)[x]
      names(TOL) <- c("UniqueName", "MatchedName", "IsSynonym", "Score", "NomenclatureCode", "IsApproxMatch", "Flags", "OTTID", "Rank", "CountsSpecies")
      return(TOL)
    }))
    
    TreeofLifeIDSynonyms <- subset(TreeofLifeIDManyMatch, IsApproxMatch=="FALSE")
    TreeofLifeIDSynonyms <- TreeofLifeIDSynonyms[ifelse(TreeofLifeIDSynonyms$UniqueName==TreeofLifeIDSynonyms$CountsSpecies,FALSE,TRUE),c(1,2)]
    names(TreeofLifeIDSynonyms) <- c("Synonym_1", "Species")
    TreeofLifeIDSynonyms <- rbindlist(lapply(unique(TreeofLifeIDSynonyms$Species), function(x){
      if(nrow(subset(TreeofLifeIDSynonyms, Species==x))==1){
        return(subset(TreeofLifeIDSynonyms, Species==x)[,c(2,1)])} else{
          Synonyms <- dcast(subset(TreeofLifeIDSynonyms, Species==x), Species~Synonym_1, value.var="Synonym_1")
          names(Synonyms) <- c("Species", paste0("Synonym_", 1:(ncol(Synonyms)-1)))
          return(Synonyms)
        }
    }), fill=TRUE)
    
    #Leave the approx matches to be manually checked
    
    ### Bring all the synonyms together
    SpeciesSynonyms <- rbind(TreeofLifeIDOneMatch, TreeofLifeIDSynonyms, fill=TRUE)}
  else {SpeciesSynonyms <- TreeofLifeIDOneMatch}
  
  SpeciesNoMatch <- as.data.frame(SpeciesList[!SpeciesList %in% SpeciesSynonyms$Species])
  names(SpeciesNoMatch) <- "Species"
  SpeciesSynonyms <- rbind(SpeciesSynonyms, SpeciesNoMatch, fill=TRUE)
  
  if(nrow(SpeciesSynonyms)!=length(SpeciesList)){stop("Somethings gone wrong with synonym matching!")}
  
  #Find all cases where original name, updated name, TOLID_ID, or a synonym matches a shapefile
  SpeciesSynonyms$FinalName <- apply(SpeciesSynonyms, 1, function (x){
    if(x["Species"] %in% SpeciesToMatch){ # species you want to match to
      x["Species"]
    } else if(x["Synonym_1"] %in% SpeciesToMatch){
      x["Synonym_1"]
    } else if(x["Synonym_2"] %in% SpeciesToMatch){
      x["Synonym_2"]
    } else {"NoMatch"}
  })
  
  SpeciesSynonyms$FinalName1 <- unlist(apply(SpeciesSynonyms, 1, function (x){
    if(x["FinalName"] != "NoMatch"){ 
      if (!all(is.na(x[c("Synonym_1","Synonym_2")]))){
        if(x["Synonym_1"] %in% SpeciesToMatch & !x["Synonym_1"] %in% x["FinalName"]){
          x["Synonym_1"]
        } else if(!is.na(x["Synonym_1"]) & x["Synonym_2"] %in% SpeciesToMatch & !x["Synonym_2"] %in% x["FinalName"]){
          x["Synonym_2"]
        } else {"NoMatch"} 
      } else {"NoMatch"}
    } else {"NoMatch"}}))
  
  return(SpeciesSynonyms)
}

## a. Species from ebird
# using taxize to get synonyms of species names
TreeofLifeID <- lapply(missing, function (x) tryCatch(as.data.frame(get_tolid_(as.character(x))), error=function(e){"NoMatch"}))
names(TreeofLifeID) <- missing

# Get any synonyms, cross check all names with master dta, return a list of either the correct name (with original name + matched name, which will either be the same or a synonym), or just the incorrect name with "nomatch"
SpeciesSynonyms <- TreeofLifeCleaning(TreeofLifeID, missing, m_info$Species)
to_replace = SpeciesSynonyms[which(SpeciesSynonyms$FinalName != 'NoMatch'),]
ebMaster_names[match(to_replace$FinalName,ebMaster_names$masterSpecies, nomatch = 0),c('ebSpecies','inferred')] =
  cbind(to_replace$Species,'ebird_TOL') # inferred 16 ebird names from TOL taxonomy


## b. Species from master dta
# using taxize to get synonyms of species names
TreeofLifeID1 <- lapply(missing_auk$masterSpecies, function (x) tryCatch(as.data.frame(get_tolid_(as.character(x))), error=function(e){"NoMatch"}))
names(TreeofLifeID1) <- missing_auk$masterSpecies

# looking for the master names in the ebird species list
SpeciesSynonyms1 <- TreeofLifeCleaning(TreeofLifeID1, missing_auk$masterSpecies, ebirdSpList$scientific_name)
to_replace1 = SpeciesSynonyms1[which(SpeciesSynonyms1$FinalName != 'NoMatch'),]
ebMaster_names[match(to_replace1$Species,ebMaster_names$masterSpecies, nomatch = 0),c('ebSpecies','inferred')] =
  cbind(to_replace1$FinalName,'ebird_TOL')

missing_auk = ebMaster_names[which(is.na(ebMaster_names$ebSpecies)),]
## 67 species unmatched

missing = setdiff(ebirdSpList$scientific_name, ebMaster_names$ebSpecies)
## still 79 species unmatched

################################################################################
## 2. with synonyms
# ITIS: Integrated Taxonomic Information Service 
# TSN : taxonomix serial number for itis

## a. Species from ebird
syn = synonyms(missing,db = "itis")
syn_df = synonyms_df(syn)
syn_df = subset(syn_df, acc_name%in%m_info$Species)[,c(".id","acc_name")]
ebMaster_names[match(syn_df$acc_name,ebMaster_names$masterSpecies, nomatch = 0),c('ebSpecies','inferred')] =
  cbind(syn_df$.id,'ebird_itis')
## 1 species
missing = missing[-which(missing%in%syn_df$.id)]


## b. Species from master dta
syn1 = synonyms(missing_auk$masterSpecies,db = "itis")
syn_df1 = synonyms_df(syn1)
syn_df1 = subset(syn_df1, syn_name%in%m_info$Species)[,c(".id","syn_name")]
# nothing


################################################################################
### 3. by hand
# looking for potential matches in the master scientific names
# check among the genus and identical species [2nd word]
taxonomyCleaning = function(names, GBIFTaxonomy){
  
  library(stringr)
  
  tax = c('class','order')
  taxM = c('Class','Order')
  
  #setDF(df)
  score = data.frame()
  splitNames = str_split(names, ' ', simplify = T) # split into genus + species
  
  ## first names is the one to be matches and subsequent names are possibilities
  matches = apply(splitNames, 1 , function(x) 
    c(paste0(x, collapse = ' '),str_subset(GBIFTaxonomy,x[1]),str_subset(GBIFTaxonomy,x[2]))) ## those species are not known to the resolution of the species
  
  for(species in 1:length(matches)){
    row = matches[[species]]
    if(length(row)>1){
      cand = unique(row[2:length(row)]) # candidate names
      nb = length(cand)
      
      score = rbind(score, 
                    data.frame(species = rep(row[1], nb), cand = cand, score = rep(0, nb)))
      #lookup = unique(subset(df, species == row[1])[,tax]) 
      # lookup = unique(subset(m_info, Species == row[1])[,taxM]) 
      #   for (level in tax){
      #   for (levelM in taxM){
      #   for (candidate in cand){
      #     if(lookup[levelM] == unique(subset(gbif_names, species == candidate)[level])){
      #         score[score$cand == candidate,'score'] = 
      #           score[score$cand == candidate,'score'] + 1 ## the closer in terms of taxonomy the better
      #     }
      #   }
      # }}
    } else {score = rbind(score, 
                          data.frame(species = row[1], cand = NA, score = NA))}
  }
  
  return(score)
}


## result of the 'by hand' matching, searching Avibase for synonyms one by one
library(readxl)
TaxMatch <- read_excel("data/eBird/Europe/Exported/0.Merged/Taxonomy/TaxonomyHomogeneisation.xlsx")

to_replace = ebMaster_names[match(TaxMatch$master,ebMaster_names$masterSpecies, nomatch = 0),]
to_replace = to_replace[apply(to_replace[,c('ebSpecies','inferred')], 1, function(x) all(is.na(x))),]
ebMaster_names[match(to_replace$masterSpecies,ebMaster_names$masterSpecies, nomatch = 0),c('ebSpecies','inferred')] = cbind(TaxMatch$ebird[match(to_replace$masterSpecies, TaxMatch$master)],'manual')

TaxMatch = TaxMatch[-which(TaxMatch$master%in%to_replace$masterSpecies),]

to_add = ebMaster_names[match(TaxMatch$master,ebMaster_names$masterSpecies, nomatch = 0),]
to_add$inferred = 'manual'
to_add[,'ebSpecies'] = TaxMatch[match(to_add$masterSpecies,TaxMatch$master, nomatch = 0),'ebird']
ebMaster_names = rbind(ebMaster_names, to_add)

missing = setdiff(ebirdSpList$scientific_name, ebMaster_names$ebSpecies)

# for cases like this one, all we need is to exchange species in 'a' w/ sp in 'us'
#colnames(master)[colnames(master)== 'Podarcis taurica'] = 'Podarcis tauricus'

score = taxonomyCleaning(missing, m_info$Species) 
score = taxonomyCleaning(missing, missing_auk$masterSpecies) 

setdiff(ebirdSpList$scientific_name, ebirdSpList1)%in%score$species

## trying some "automatic" matching which obv does not add anything
for (i in unique(score$species[which(!is.na(score$cand))])){
  cand = subset(score, species == i)[,'cand']
  name = cand[sapply(cand, function(x) substr(i,1,nchar(i)-1) == substr(x,1,nchar(x)-2))]
  
  if (length(name)==1){
    print(paste(i, name))
    #colnames(master)[colnames(master)== i] = name
    next
  } 
  cand1 = cand[nchar(cand) == nchar(i)]
  if (length(cand1)>0){
    for (candidate in cand1){
      mat = table(unlist(str_split(candidate,'')),unlist(str_split(i,'')))
      if(sum(mat)-sum(diag(mat)) <= 3){
        print(paste(candidate, i))
        #colnames(master)[colnames(master)== i] = candidate
        next
      }
    }
    
  } 
  # if (i %in% gbif_names$canonicalName){
  #   candidate = gbif_names[gbif_names$canonicalName == i,'species']
  #   print(paste(candidate,i))
  # }
}

################################################################################
# save the result

write.csv(ebMaster_names,'data/eBird/Europe/Exported/0.Merged/Taxonomy/ebirdMay-2022.SpeciesMatched.master.csv')

paste0('There are ',length(unique(ebMaster_names$ebSpecies)), ' unique ebird species for whom we have master codes')
paste0('There are ',length(which(is.na(ebMaster_names$ebSpecies))), ' master species for whom we have no correspondance to ebird species')
# those seem to mostly occur outside of non-russia europe though


################################################################################
# 3.1 Keep only metaweb species and save eb data and create data frame of all the eb checklists 
# 3.2 Extract grid cell ID for each checklist
# 3.3 Merge grid cell ID with the ebird occurrence data using checklist ID and subset breeding season only
################################################################################

# 3.1 Keep only metaweb species and save eb data and create data frame of all the eb checklists 

europe = rnaturalearth::ne_countries(continent = "Europe", returnclass = "sf", scale = 'large')
europe = europe[-which(europe$admin == 'Russia'),]

speciesList = read.csv('data/eBird/Europe/Exported/0.Merged/Taxonomy/ebirdMay-2022.SpeciesMatched.master.csv')
metawebSp = speciesList[!is.na(speciesList$ebSpecies),] # species whom I managed to match
anyNA(metawebSp$masterCode) # should be false

files = str_subset(list.files("data/eBird/Europe/Exported/0.Merged/"),'Merged') ## 3 files


## function that lists all the ebird checklists and their characteristics (used in loop below)
## Need to make sure that all the data for each country is togetehr so that no checklists are split
## START FUNCTION
chlist_creation = function(eb, n = 4){ #splitting the data so that it doesn't reach memory limit
  
  chlist = data.frame()
  list = unique(eb$checklist)
  pos = round(seq(0,length(list),length.out = n))
  
  for (i in 1:(length(pos)-1)){
    
    start = pos[i] + 1
    stop = pos[i + 1]
    sublist = list[start:stop]
    
    chlist<-rbind(chlist,
                  ddply(subset(eb, checklist %in% sublist), ~ checklist, function(x){
                    data.frame(duration=mean(x$duration_minutes), distance=mean(x$effort_distance_km), year=mean(x$year), day=mean(x$day), 
                    lat=mean(x$latitude), lon=mean(x$longitude), protocol=x$protocol_type[1], Country=x$country_code[1],
                    Observers=paste(unique(unlist(strsplit(x$observer_id, ","))), collapse=";"), richness = length(unique(x$masterCode)))}))
    
  }
  return(chlist)
} ## END FUNCTION




for (f in 1:length(files)){
  
  eb = readRDS(paste("data/eBird/Europe/Exported/0.Merged/",files[f]))
  
  intersect(metawebSp$ebSpecies, unique(eb$scientific_name))
  setdiff(metawebSp$ebSpecies, unique(eb$scientific_name)) 
  # there are 4 birds from the metaweb that were removed when filtering above
  
  print(paste('dimension before removing sp not in metaweb', dim(eb)))
  print(dim(eb[which(eb$scientific_name %in% setdiff(unique(eb$scientific_name), metawebSp$ebSpecies))]))
  eb = merge(eb, metawebSp, by.x = 'scientific_name', by.y = 'ebSpecies')
  print(paste('dimension after removing species not in metaweb', dim(eb)))
  
  if(anyDuplicated(eb)>0){
    eb = eb[-which(duplicated(eb)),]
    print(paste('dimension after removing duplicates', dim(eb)))
    }
  
  ################################################################################
  ### create data frame with info on all checklists (inspired from 2.Format.eBird.Cazalis.R)
  
  names(eb)[names(eb) == "checklist_id"] = "checklist"
  chlist = chlist_creation(eb) # function that summarizes the info about each individual checklist event
  colnames(chlist)[1]<-"Liste"
  print(paste('Should be TRUE:', nrow(chlist) == length(unique(eb$checklist))))
  
  ## save ebird occurrence data
  saveRDS(eb, paste0("data/eBird/Europe/Exported/0.Merged/masterSpecies/Merged.data.",version,'.',f,'.MetawebSp.rds'))
  ## save info about checklists
  write.csv(chlist, paste0('data/eBird/Europe/Exported/0.Merged/chlist/checklist.ebird-',version,'.',f,'.Richness.csv'))
}


################################################################################
# 3.2 Extract grid cell ID (10x10km) for each checklist

setwd("data/eBird/Europe/Exported/0.Merged/chlist")

files = list.files()[which(str_detect(list.files(),'checklist'))][1:3]

### combine checklists
chlist = rbindlist(lapply(files, read.csv))
# chlist = readRDS(paste0('checklist.ebird-',version,'.Richness.Dist2km.rds'))

## grid of europe at a 10x10km resolution
europeRaster <- raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")
cells_info <- read.dbf('data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img.vat.dbf')
cells_info = cbind(cells_info, coordinates(europeRaster)[which(!is.na(raster::values(europeRaster))),])

## into a spatial feature
sf_df = st_as_sf(x = chlist, coords = c("lon","lat"), crs = 4326) ## spatial points for the occurrence data GBIF
# europeRaster.4326 = projectRaster(europeRaster, crs = 4326)
sf_df = st_transform(sf_df, crs = crs(europeRaster))

grid = raster::extract(europeRaster, sf_df, df = T) # extract the 10x10 km grid cell in which the checklist was made
sf_df$Value = grid$reference_grid_10km 

paste0('there are ',length(which(is.na(grid$reference_grid_10km))),' checklist that fall outside of the master grid, or ',round(length(which(is.na(grid$reference_grid_10km)))/nrow(chlist),3),'% of the checklists')

# remove the chelists with missing values
sf_df = subset(sf_df, !is.na(Value))

## save the file 
## 23/01/2023 re ran the code without filtering on distance between survey events on the same day to match what is done with the GBIF data (above)
write.csv(st_drop_geometry(sf_df), paste0('data/eBird/Europe/Exported/0.Merged/chlist/checklist.ebird-',version,'.GridValues.23012023.csv'))

### plots
europe = st_transform(europe, crs = crs(europeRaster))

ggplot() + geom_sf(data = sf_df[which(is.na(sf_df$Value)),]) + geom_sf(data = europe) +
  xlim(extent(sf_df)[c(1:2)]) + ylim(extent(sf_df)[c(3:4)])# the NA data seem to fall outside of the continent
ggplot() + geom_sf(data = europe) + geom_sf(data = sf_df[sample(nrow(sf_df), nrow(sf_df)*0.1),]) +
  xlim(extent(sf_df)[c(1:2)]) + ylim(extent(sf_df)[c(3:4)]) 

################################################################################
# 3.3 Merge grid cell ID with the ebird occurrence data using checklist ID and subset breeding season only

## eb = rbindlist(lapply(str_subset(list.files(), 'Merged')[1:3], readRDS), fill = T)
## too big to put in one go
sf_df = read.csv(paste0('data/eBird/Europe/Exported/0.Merged/chlist/checklist.ebird-',version,'.GridValues.23012023.csv'))

files = str_subset(list.files("data/eBird/Europe/Exported/0.Merged/masterSpecies"), 'Merged')

months = c(4,5,6,7,8) ## breeding season
nb = 0
df_unique = data.frame()
season = 'breeding'

for (f in 1:length(files)){
  
  eb = readRDS(files[f])
  print(dim(eb))
  # add grid cell info to ebird data
  eb <- merge(eb, st_drop_geometry(sf_df[,c('Liste','Value')]), by.x = 'checklist', by.y = 'Liste')
  print(dim(eb))
  
  ### create eventList
  eb$date = paste(eb$day, eb$year)
  eventList = unique(eb[,c('date','Value')])
  names(eventList)[which(names(eventList) == 'date')] = "eventDate"
  names(eventList)[which(names(eventList) == 'Value')] = "df_grid"
  eventList$datasetKey = 'eBird';eventList$class = 'Aves'
  
  write.csv(eventList, paste0("outputs/EventList/EventList.ebird",version,".",f,".23012023.csv"))
  ### Save tables 
  saveRDS(eb, paste0('data/eBird/Europe/Exported/0.Merged/chlist/eBird.Grid/eBird.data.',version,'.',f,'.23012023.rds'))
  
  ### create species checklist per grid cell for breeding season
  df1 = plyr::ddply(subset(eb, month %in% months), 'Value', summarise, Species = unique(masterCode))
  df_unique = rbind(df_unique,df1)
  
  ## check
  nb = c(nb, unique(subset(eb, month %in% months)$Value))
  print(paste0('Should be TRUE: ', length(unique(nb))-1 == length(unique(df_unique$Value))))
  
}

df_unique$pres = 1

# list of species
write.csv(df_unique, paste0('data/eBird/Europe/Exported/0.Merged/LongSpeciesList.masterValues/SpeciesList.Values.', season,'.23012023.csv'))
# 458 ebird species
