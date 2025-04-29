### code for 'Mixed effect of protected areas on avian food webs'

# Creating figures and variable selection of environmental drivers of PA effectiveness
## 05/12/2022

## Loading packages ----
library(data.table)
library(stringr)
library(plyr)
library(raster)
library(mgcv)
library(MASS)
library(plyr)
library(corrplot)
library(DHARMa)
library(car)
library(lme4)
library(ggplot2)
library(ggpubr)
library(vegan)
library(sf)

# colour palette
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

setwd("E:/TheseSwansea/WDPA")
source('Script/FunctionsAnalysis.R') # source functions for analysis


## Loading data ----
## 'res' is the dataframe with network metrics calculated at the grid cell level 
res = readRDS("outputs/1.Metrics/NetworkMetrics.Rarefied.nb15.effort50.Over5yo.AfterDesig.NewEffort.rds")
res = merge(res, readRDS("outputs/1.Metrics/NetworkMetrics.Rarefied.nb15.effort50.Over5yo.AfterDesig.NewEffort.BodyMass.rds"))

# communities with no top species or intermediate species
res$mass_top[which(is.na(res$mass_top))] = 0
res$mass_int[which(is.na(res$mass_int))] = 0

#### FW metrics of interest
Indice.interest = c("S","C","L.S", "L", "modularity.2", "mfcl", "omnivory.1", "Gen", 
                    "Vul", "TL.mean","mass_top", "mass_basal","mass_int",# "Fbasal","FInt","Ftop",
                    "S_basal","S_intermediate","S_top")


# library(corrplot)
# library(RColorBrewer)
# M <- cor(res[,Indice.interest], use = 'complete.obs')
# corrplot(M, type="upper", order="hclust", col=brewer.pal(n=8, name="RdYlBu"))
# # corrplot(M[c("Gen","mfcl","S","FInt"),c("Gen","mfcl","S","FInt")], type="upper", order="hclust", col=brewer.pal(n=8, name="RdYlBu"))

## list of rarefied survey events that were kept for analysis - 
## this is each checklist of bird observation:
## (1) that falls within on of our PA networks
## (2) that had sufficient survey effort (>50) 
## (3) that was picked in the rarefaction
eventList.min.sub = readRDS('outputs/EventList/EventListMerged/Rarefying/EventListRarefied.Over5yo.NewEffort.1.rds')
eventList.min.sub = subset(eventList.min.sub, class == 'Aves')
eventList.min.sub = subset(eventList.min.sub, ID.new.Bioregion != "2_Continental Bio-geographical Region")

length(unique(eventList.min.sub$ID.new.Bioregion)) ## should be 45


### Adding environmental variables to food web metrics -----

## raster of extent of TETRA EU database
europeRaster <- raster::raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")
## info about each cell from TETRA EU + combine with raster
cells_info <- foreign::read.dbf('data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img.vat.dbf')
cells_info = cbind(cells_info, coordinates(europeRaster)[which(!is.na(raster::values(europeRaster))),])

# get shp of europe and transform to right crs
europe = rnaturalearth::ne_countries(returnclass = "sf", scale = 'large')
europe = sf::st_transform(europe, crs(europeRaster))
europe = sf::st_crop(europe, extent(europeRaster))

# from 5.Analysis - environmental variables formatted to the grid cell level 
load("outputs/ContextVariables.Over5yo.RData")



## We devise a coarse quantification of the level of protection based on IUCN PA categories
levels = data.frame(IUCN_CAT = levels(factor(PA.grid.nodups$IUCN_CAT)), score = c(10,9,8,7,6,1,1,1,5,4))
PA.grid.nodups = merge(PA.grid.nodups, levels, by = 'IUCN_CAT')
PA.grid.nodups$Directed_at_birds = 
  ifelse(PA.grid.nodups$DESIG_ENG %in% c("Ramsar Site, Wetland of International Importance", "Special Protection Area (Birds Directive)"),
         1, 0)

## CREATE A TABLE WITH ALL THE CARACTERISTICS OF EACH GRID PA
## merged PAs hasd the info about which grid cells fall into which PA network
merged.PAs = readRDS("outputs/GridCells.NeighbouringPAID.new.Over5yo.RDS")
cells.to.seperate = c('555534657.1','555534704.1','555592567.1')
merged.PAs$NeighbouringPAID.new[which(merged.PAs$NeighbouringPAID %in% cells.to.seperate)] = '5.1'

## Need to label non protected grid cells in 'res' 
res$NeighbouringPAID = merged.PAs$NeighbouringPAID[match(res$Value,merged.PAs$Value)]

## ADD environmental variables to food web metrics
# Only grid cells with occurrence data (and network metrics):
res = add_context_var(res = res)

setDF(res)
res = res[,-which(str_detect(names(res), 'NA'))] # remove NA columns (NA landcover types)
res = subset(res, STATUS_YR != 0 | is.na(STATUS_YR)) # remove grid cells with no info on year of implementation ('0') of protection (but keep NAs for non-protected grid cells)
res = res[-which(res$S == 0),] # remove one community which had non interacting birds

# Look at map of remoteness:
ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = Remoteness), size = 0.8) + scale_colour_gradientn(colors = rainbow(5))

# ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = m.gHP), size = 0.8) + scale_colour_gradientn(colors = rainbow(5))
# ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = m.THPI))
# ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = score), size = 0.8) + scale_colour_gradientn(colors = rainbow(5))


## We calculate PA characteristics using all grid cells within 100 km of the PA 
## (not just those with sufficiently high survey effort, this is why we don't use 'res' above)

## add env variables to the grid cells (all grid cells in a site, not just surveyed grid cells with bird occurrence data)
merged.PAs = add_context_var(res = merged.PAs[,c('Value','dist.PA.min','Bioregion','NeighbouringPAID.new','x','y')])

## remove missing values for the env variables
merged.PAs = subset(merged.PAs, !(is.na(LC) | LC == 'lc.NA'))
merged.PAs = subset(merged.PAs, !(is.na(m.elev)|is.na(m.slope)|is.na(m.gHP)|is.na(Remoteness)))
merged.PAs = subset(merged.PAs, STATUS_YR != 0 | is.na(STATUS_YR))
merged.PAs = merged.PAs[,-which(str_detect(names(merged.PAs), 'NA'))]

## keep 100 km closest cells and only bioregions with sufficient data
merged.PAs = subset(merged.PAs, dist.PA.min <= 100e3 & ID.new.Bioregion %in% res$ID.new.Bioregion)
merged.PAs$IncludedInAnalysis = ifelse(merged.PAs$Value %in% res$Value, "Yes", "No")

## 2 continental: this site is in the east and CLC doesn't cover that far east
removed = subset(merged.PAs, ID.new.Bioregion == "2_Continental Bio-geographical Region")
merged.PAs = subset(merged.PAs, ID.new.Bioregion != "2_Continental Bio-geographical Region")


length(unique(merged.PAs$ID.new.Bioregion)) 

## Calculate the site-level context metrics ----

summary(merged.PAs$STATUS_YR[which(merged.PAs$PA == 1)]) ## Should be NO 'NAs'

# anthropised land cover types 
human_hab = c("lc.Continuous urban fabric", "lc.Discontinuous urban fabric", "lc.Industrial or commercial units",
              "lc.Road and rail networks and associated land", "lc.Port areas","lc.Airports","lc.Mineral extraction sites",
              "lc.Dump sites" ,"lc.Construction sites","lc.Green urban areas","lc.Sport and leisure facilities", 
              "lc.Non-irrigated arable land", "lc.Permanently irrigated land", "lc.Rice fields", "lc.Vineyards" , 
              "lc.Fruit trees and berry plantations", "lc.Olive groves", "lc.Pastures","lc.Annual crops associated with permanent crops",
              "lc.Complex cultivation patterns", "lc.Land principally occupied by agriculture, with significant areas of natural vegetation",
              "lc.Agro-forestry areas")
natural_hab = str_subset(names(merged.PAs),'lc\\.')[which(!str_subset(names(merged.PAs),'lc\\.') %in% human_hab)]

urban_hab = c("lc.Continuous urban fabric", "lc.Discontinuous urban fabric", "lc.Industrial or commercial units",
              "lc.Port areas","lc.Airports","lc.Construction sites","lc.Green urban areas", "lc.Sport and leisure facilities")

agricultural_hab = c(
  "lc.Non-irrigated arable land", "lc.Permanently irrigated land", "lc.Rice fields", "lc.Vineyards" , 
  "lc.Fruit trees and berry plantations", "lc.Olive groves", "lc.Pastures","lc.Annual crops associated with permanent crops",
  "lc.Complex cultivation patterns", "lc.Land principally occupied by agriculture, with significant areas of natural vegetation",
  "lc.Agro-forestry areas")



# summarise the proportion of land cover for each category of habitat from lists above
merged.PAs$human_hab = rowSums(merged.PAs[, human_hab])
merged.PAs$urban_hab = rowSums(merged.PAs[, urban_hab])
merged.PAs$agricultural_hab = rowSums(merged.PAs[, agricultural_hab])
merged.PAs$natural_hab = rowSums(merged.PAs[, natural_hab])

####

# agri = tidyr::pivot_longer(data = merged.PAs[,c(str_subset(names(merged.PAs), "lc\\."),"PA")], 
#                            cols = str_subset(names(merged.PAs), "lc\\."))
# agri = ddply(agri, c("name","PA"), summarise, value = mean(value))
# agri = ddply(agri, 'PA', summarise, cum = rev(cumsum(rev(value))), name = name, value = value)

# ggplot(data = agri, aes(y = value, x = PA, fill = name, label = name)) + 
#   # theme(axis.text.y=element_text(colour = colour)) +
#   geom_bar(stat = 'identity') + 
#   geom_text(data = agri, aes(y = (cum + cum - value)/2, x = PA, label = ifelse(value > 0.03,paste(name, round(value,2)), "")), size = 5) +
#   xlab("") + ylab("Proportion of land cover") 
# 
# ggsave("Plots/PropLC.png", width = 20, height = 10)


agri = tidyr::pivot_longer(merged.PAs[,c(agricultural_hab,"PA")], cols = all_of(agricultural_hab))
agri = ddply(agri, c("name","PA"), summarise, value = mean(value))
agri = ddply(agri, 'PA', summarise, cum = rev(cumsum(rev(value))), name = name, value = value)

# ggplot(data = agri, aes(y = value, x = PA, fill = name, label = name)) + 
#   # theme(axis.text.y=element_text(colour = colour)) +
#   geom_bar(stat = 'identity') + 
#   geom_text(data = agri, aes(y = (cum + cum - value)/2, x = PA, label = ifelse(value > 0.01,paste(name, round(value,2)), "")), size = 5) +
#   xlab("") + ylab("Proportion of land cover") 
# 
# ggsave("Plots/PropLC.png", width = 20, height = 10)


setDF(merged.PAs)
# calculate land cover diversity
merged.PAs$LC.natural.diversity = vegan::diversity(merged.PAs[,natural_hab], MARGIN = 1)
merged.PAs$LC.diversity = vegan::diversity(merged.PAs[,str_subset(names(merged.PAs),'lc\\.')], MARGIN = 1)

names(merged.PAs) = str_replace_all(names(merged.PAs), ' ', '')

drivers = ddply(merged.PAs, "ID.new.Bioregion", summarise, NeighbouringPAID.new = unique(NeighbouringPAID.new), 
                 nb.inside.all = length(which(PA == 1)), nb.outside.all = length(which(PA == 0)), nb.all = unique(nb.inside.all + nb.outside.all), 
                 prop.area.protected = nb.inside.all/nb.all,
                 
                 
                 x.mean = mean(x), y.mean = mean(y), nb.PA = length(unique(ID[which(!is.na(ID))])),
                 
                 gHP.m = mean(m.gHP, na.rm = T), 
                 gHP.m.outside = mean(m.gHP[which(PA == 0)], na.rm = T),
                 gHP.m.inside = mean(m.gHP[which(PA == 1)], na.rm = T),
                 gHP.diff = gHP.m.inside - gHP.m.outside,
                 
                 gHP.sd = sd(m.gHP, na.rm = T),
                 gHP.sd.outside = sd(m.gHP[which(PA == 0)], na.rm = T),
                 gHP.sd.inside = sd(m.gHP[which(PA == 1)], na.rm = T),
                 
                 HD.m = mean(m.HD, na.rm = T), 
                 HD.m.outside = mean(m.HD[which(PA == 0)], na.rm = T),
                 HD.m.inside = mean(m.HD[which(PA == 1)], na.rm = T),
                 HD.diff = HD.m.inside - HD.m.outside,
                 
                 remoteness.mean = mean(Remoteness), remoteness.sd = sd(Remoteness), 
                 remoteness.m.inside = mean(Remoteness[which(PA == 1)]), remoteness.m.outside = mean(Remoteness[which(PA == 0)]),
                 remoteness.diff = remoteness.m.inside - remoteness.m.outside,
                 
                 # elevation
                 elev.mean = mean(m.elev), slope.mean = mean(m.slope,na.rm = T), 
                 elev.mean.inside = mean(m.elev[which(PA == 1)]), 
                 elev.mean.outside = mean(m.elev[which(PA == 0)]), 
                 elev.diff = elev.mean.inside - elev.mean.outside,
                 

                 # Slope
                 slope.mean.inside = mean(m.slope[which(PA == 1)]), 
                 slope.mean.outside = mean(m.slope[which(PA == 0)]), 
                 slope.diff = slope.mean.inside - slope.mean.outside,
                 
                 # Human hab 
                 area.human.hab = mean(human_hab), sd.area.human.hab = sd(human_hab),
                 area.human.hab.inside = mean(human_hab[which(PA == 1)]), 
                 area.human.hab.outside = mean(human_hab[which(PA == 0)]),
                 area.human.hab.diff = area.human.hab.inside - area.human.hab.outside,
                 
                 # urban hab 
                 area.urban.hab = mean(urban_hab), 
                 area.urban.hab.inside = mean(urban_hab[which(PA == 1)]), 
                 area.urban.hab.outside = mean(urban_hab[which(PA == 0)]),
                 area.urban.hab.diff = area.urban.hab.inside - area.urban.hab.outside,
                 
                 # natural hab 
                 area.natural.hab = mean(natural_hab), 
                 area.natural.hab.inside = mean(natural_hab[which(PA == 1)]), 
                 area.natural.hab.outside = mean(natural_hab[which(PA == 0)]),
                 area.natural.hab.diff = area.natural.hab.inside - area.natural.hab.outside,
                 
                 # agricultural hab 
                 area.agricultural.hab = mean(agricultural_hab), 
                 area.agricultural.hab.inside = mean(agricultural_hab[which(PA == 1)]), 
                 area.agricultural.hab.outside = mean(agricultural_hab[which(PA == 0)]),
                 area.agricultural.hab.diff = area.agricultural.hab.inside - area.agricultural.hab.outside,
                 
                 ## Forest
                 area.forest = mean(lc1.Forestandseminaturalareas),
                 area.forest.outside = mean(lc1.Forestandseminaturalareas[which(PA == 0)]),
                 area.forest.inside = mean(lc1.Forestandseminaturalareas[which(PA == 1)]),
                 area.forest.diff = area.forest.inside - area.forest.outside,
                 
                 # Waterbodies
                 area.Waterbodies.outside = mean(lc1.Waterbodies[which(PA == 0)]),
                 area.Waterbodies.inside = mean(lc1.Waterbodies[which(PA == 1)]),
                 area.Waterbodies.diff = area.Waterbodies.inside - area.Waterbodies.outside,
                 
                 # Wetlands
                 area.Wetlands.outside = mean(lc1.Wetlands[which(PA == 0)]),
                 area.Wetlands.inside = mean(lc1.Wetlands[which(PA == 1)]),
                 area.Wetlands.diff = area.Wetlands.inside - area.Wetlands.outside)

# PA characteristics
PA.chara = unique(merged.PAs[,c("ID.new.Bioregion","ID","area","score","Directed_at_birds","STATUS_YR")])
PA.chara = PA.chara[-which(is.na(PA.chara$ID)),]
PA.chara = ddply(PA.chara, "ID.new.Bioregion", summarise, 
                 area.sum = sum(area, na.rm = T), ## some mega PAs are only constituted of 1 PA so no sd 
                 area.mean = mean(area, na.rm = T),
                 year.mean = mean(STATUS_YR, na.rm = T), # year.sd = sd(STATUS_YR, na.rm = T),
                 score.m = mean(score, na.rm = T), #score.sd = sd(score, na.rm = T), 
                 prop.forBirds = length(unique(ID[which(Directed_at_birds == 1)]))/length(unique(ID[which(!is.na(ID))])))

drivers = merge(drivers, PA.chara, by = "ID.new.Bioregion")

## calculating land cover diversity inside and outside PAs for natural 
## land cover
dt_div = setDT(merged.PAs)[, lapply(.SD, sum), by = c('ID.new.Bioregion','PA'), .SDcols = str_subset(names(merged.PAs), 'lc')]
setDF(dt_div)
dt_div$diversity.natural[which(dt_div$PA == 1)] = vegan::diversity(dt_div[which(dt_div$PA == 1),str_replace_all(natural_hab,' ', '')], MARGIN = 1)
dt_div$diversity.natural[which(dt_div$PA == 0)] = vegan::diversity(dt_div[which(dt_div$PA == 0),str_replace_all(natural_hab,' ', '')], MARGIN = 1)
dt_div = tidyr::pivot_wider(dt_div[,c('ID.new.Bioregion','PA','diversity.natural')], names_from = PA, values_from = diversity.natural)
names(dt_div) = c('ID.new.Bioregion','diversity.natural.LC.outside','diversity.natural.LC.inside')
drivers = merge(dt_div, drivers, by = 'ID.new.Bioregion', all.y = T)

## calculating land cover diversity for the whole site
dt_div1 = setDF(setDT(merged.PAs)[, lapply(.SD, sum), by = ID.new.Bioregion, .SDcols = str_subset(names(merged.PAs), 'lc')])
dt_div1$diversity.LC = vegan::diversity(dt_div1[,str_subset(names(dt_div1),'lc\\.')], MARGIN = 1)
drivers = merge(dt_div1[,c('ID.new.Bioregion','diversity.LC')], drivers, by = 'ID.new.Bioregion', all.y = T)

print(paste("Should be False :", anyNA(drivers)))


drivers.res = ddply(res, "ID.new.Bioregion", summarise, 
                 nb.inside = length(which(PA == 1)), nb.outside = length(which(PA == 0)), nb = unique(nb.inside + nb.outside))

drivers = merge(drivers, drivers.res, by = 'ID.new.Bioregion')
drivers = subset(drivers, nb.inside >= 15 & nb.outside >= 15)
## 5 boreal : very large PA with very few protected grid cells
## 188 : small PA (protecting a lake) with few protected points and many non-protected points

res = subset(res, ID.new.Bioregion %in% drivers$ID.new.Bioregion)
summary(res)
length(unique(res$ID.new.Bioregion))
 
drivers.list = c('year.mean','prop.forBirds','log.area.sum','frag',
                 'log.remoteness.m','area.urban.hab.outside','area.urban.hab.diff',
                 'area.agricultural.hab.outside','area.agricultural.hab.diff',
                 'log.HD.inside','log.HD.outside','HD.diff',
                 # 'gHP.diff', 'gHP.m.inside','gHP.m.outside',
                 'log.area.forest.inside','area.forest.diff','log.area.forest.outside',
                 'diversity.LC','log.slope','elev.mean','score.m',"log.area.mean","log.nb")

ggplot(data = res, aes(y = S/effort_Aves, x = paste(Bioregion,PA), fill = factor(PA))) + 
  geom_bar(stat = 'summary', fun = 'mean') + 
  scale_x_discrete(labels = c("Alpine Bio-geographical Region", '',"Boreal Bio-geographical Region",'',
                              "Atlantic Bio-geographical Region", '', "Continental Bio-geographical Region", '',
                              "Mediterranean Bio-geographical Region",'')) 



## Importing outputs from the GAMs (6.Analysis.WDPA) - difference in food web metrics ----

## from Script/proceedings/6.....MixedEffects
contrast_res = read.csv('outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/randomEffect-FullDataset-RandomSlope.csv')
names(contrast_res)[names(contrast_res) == 'group'] = 'ID.new.Bioregion'
contrast_res = subset(contrast_res, effect == 'PA') ## keep only effects related to protection effects (ignoring all the control variables)
contrast_res = subset(contrast_res, Model == 'Species richness' | (Model =='No species richness' & Indice == 'S'))

fixed_effect = read.csv('outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/fixedEffect-FullDataset-RandomSlope.csv')
fixed_effect = subset(fixed_effect, Model == 'Species richness' | (Model =='No species richness' & Indice == 'S'))
fixed_effect = subset(fixed_effect, term == 'PA')

names(fixed_effect)[names(fixed_effect) == 'value'] = 'value_fixed_effect'
names(fixed_effect)[names(fixed_effect) == "p_value"] = 'p_value_fixed_effect'
names(fixed_effect)[names(fixed_effect) == "upper_97.5"] = "upper_97.5_fixed_effect"
names(fixed_effect)[names(fixed_effect) == "lower_2.5"] = "lower_2.5_fixed_effect"

contrast_res = merge(fixed_effect[,c('value_fixed_effect','p_value_fixed_effect',
                                     'Indice','Model',"upper_97.5_fixed_effect","lower_2.5_fixed_effect")], 
                     contrast_res, by = c('Indice','Model'))

## random slope 
## se and CI are extracted from the covariance matrix
contrast_res$absolute_value = contrast_res$value_fixed_effect + contrast_res$value

contrast_res$rd_slope_significance = ifelse((contrast_res$lower_2.5 + contrast_res$value_fixed_effect) *
                                              (contrast_res$upper_97.5 + contrast_res$value_fixed_effect) > 0, '***', '')

## check that we still have all 45 sites and 21 food web metrics
length(unique(contrast_res$ID.new.Bioregion))
ddply(contrast_res, "ID.new.Bioregion", summarise, nb = length(Indice))[c("ID.new.Bioregion","nb")]

## keep only 16 FW metrics of interest
contrast_res = subset(contrast_res, Indice %in% Indice.interest)
## and 45 PA networks
contrast_res = subset(contrast_res, ID.new.Bioregion %in% drivers$ID.new.Bioregion)



## Create clean names for all variables ----

library(tidyr)
setDF(contrast_res)

## 'translate' variable names for paper :
clean_Indices = list("Gen"="Mean generality", "FInt"="Fraction of intermediate sp", "S" = 'Species richness',"mfcl" = "Mean food chain length",
                     "L.S" = "Linkage density", "TL.mean" = 'Mean trophic position', "Fbasal" = "Fraction of basal sp", "SDVul" = "Standard deviation of Vulnerability", 
                     "L" = "Number of links", "omnivory.1" = 'Omnivory', "Vul" = "Mean vulnerability", "C" = "Connectance", "modularity.2" = "Modularity",
                     "SDGen" = "Standard deviation of generality","Ftop" = "Fraction of top sp","S_top" = "Number of top species",
                     "S_basal" = "Number of basal species",'S_intermediate' = "Number of intermediate species",
                     "mass_basal" = "Mean mass of basal species", "mass_top" = "Mean mass of top species",
                     "mass_int" = "Mean mass of intermediate species")

clean_Drivers = list("(Intercept" = "Intercept", "log.N" = "Number of communities (logged)","gHP.diff" = "Difference in human pressure (inside - outside))",
                     "diversity.LC" = "Diversity of land cover types","score.m" = "Protection level", "frag" = "Protected area continuity",
                     "log.area.sum" = "Total protected area (logged)", "frag.LC.inside" = "Land cover continuity", "gHP.m.outside" = "Mean human pressure outside",
                     "log.area.forest.inside" = "Proportion of forested area inside (logged)", "elev.mean" = "Mean site elevation", "gHP.m.inside" = "Mean human pressure inside",
                     "log.area.forest.outside" = "Proportion of forested area outside (logged)", "year.mean" = "Mean age of the protected areas", 
                     "area.agricultural.hab.diff" = "Difference in proportion of agricultural area", "log.remoteness.m" = "Mean remoteness of the site (logged)",
                     "area.agricultural.hab.outside" = "Proportion of agricultural area outside", "log.slope" = "Mean slope of the site (logged)",
                     "area.forest.diff" = "Difference in the proportion of forested area", 
                     "area.urban.hab.outside" = "Proportion of urban habitat outside",
                     "prop.forBirds" = "Proportion of protected areas managed for birds", 
                     "area.urban.hab.diff" = "Difference in proportion of urban habitat","log.HD.inside" = "Human density inside (logged)",
                     "log.HD.outside" = "Human density outside (logged)", "HD.diff" = "Difference in human density",
                     "prop.area.protected" = "Proportion of protected area per site")

summary(contrast_res)

contrast_res$clean_Indice = unlist(clean_Indices[contrast_res$Indice])


# Plotting random effects #####################################################

## point plots with confidence interval
## Fixed effect on FW metrics
ggplot(data = unique(contrast_res[!str_detect(contrast_res$Indice, 'mass'),
                                  c('value_fixed_effect','Indice',
                                    "upper_97.5_fixed_effect",
                                    "lower_2.5_fixed_effect")]), 
                     aes(y = value_fixed_effect, x = Indice)) + 
  geom_point(stat = 'identity') +
  geom_errorbar(aes(ymin = lower_2.5_fixed_effect, ymax = upper_97.5_fixed_effect)) +
  geom_hline(yintercept = 0, col = 'red') + 
  coord_flip() + 
  ggtitle('Fixed effect')
  
# body mass
ggplot(data = unique(contrast_res[str_detect(contrast_res$Indice, 'mass'),
                                  c('value_fixed_effect','Indice',
                                    "upper_97.5_fixed_effect",
                                    "lower_2.5_fixed_effect")]), 
       aes(y = value_fixed_effect, x = Indice)) + geom_point() +
  geom_errorbar(aes(ymin = lower_2.5_fixed_effect, ymax = upper_97.5_fixed_effect)) +
  geom_hline(yintercept = 0, col = 'red') + 
  coord_flip()


### plots of random effects with shaded fixed effect in the middle
for(indice in unique(contrast_res$Indice)){
  
  sub = subset(contrast_res, Indice == indice)
  
  ggplot(data = sub, aes(x = absolute_value, y = reorder(ID.new.Bioregion, absolute_value))) +
    ## create box with CI around fixed effect
    geom_rect(aes(
      xmin = unique(lower_2.5_fixed_effect),
      xmax = unique(upper_97.5_fixed_effect),
      ymin = -Inf,
      ymax = Inf
    ), alpha = 0.1) +
    ## add line with fixed effect
    geom_vline(aes(xintercept = unique(value_fixed_effect)),
               color = "pink",
               size = 1) +
    ## random effects
    geom_errorbar( ## errorbar around absolute value of random effect
      aes(
        xmin = lower_2.5 + value_fixed_effect,
        xmax = upper_97.5 + value_fixed_effect
      ),
      width = 0,
      size = 1
    ) +
    geom_point(
      size = 3,
      shape = 21,
      color = "black",
      fill = "white"
    ) +
    geom_vline(
      xintercept = 0,
      color = "red",
      size = 1,
      linetype = "dashed"
    ) +
    
    labs(x = "Intercept",
         y = "Subject ID", 
               title = paste("Random Slopes (PA effect) - ", indice))
  
  ggsave(paste0("Plots/",indice,"_ranef-RandomSlope-NoRandomIntercept.png"), width = 10, height = 8, dpi = 1000)

}


# Figure 3 -----
### Plot of random effect as distributions ----

contrast_res$fixed_effect_significance = ifelse(contrast_res$p_value_fixed_effect < 0.001, 'p-value < 0.001', 
                                                   ifelse(contrast_res$p_value_fixed_effect < 0.01, 'p-value < 0.01',
                                                          ifelse(contrast_res$p_value_fixed_effect < 0.05, 'p-value < 0.05',
                                                                 ifelse(contrast_res$p_value_fixed_effect < 0.1, 'p-value < 0.1',''))))


g1 = ggplot(data = droplevels.data.frame(contrast_res[!contrast_res$Indice %in% 
                                                           c('mass_int','mass_basal','mass_top','omnivory.1','S_basal'),]), 
            aes(x = absolute_value, fill = rd_slope_significance)) +
  geom_vline(aes(xintercept = 0),
             color = "black", linewidth = 1, linetype = 3) + ## dotted line for vertical line in zero
  geom_histogram(stat = 'bin') +
  geom_rect(aes(
    xmin = lower_2.5_fixed_effect,
    xmax = upper_97.5_fixed_effect,
    ymin = -Inf,
    ymax = Inf
  ), fill = "grey", alpha = 0.01, show.legend = FALSE) +
  geom_vline(aes(xintercept = value_fixed_effect),
             color = "grey", linewidth = 1.5) +
  facet_wrap(~clean_Indice, scales = 'free', ncol = 3) + 
  theme_classic() + 
  xlab('') + 
  scale_fill_manual(values= c("#DDCC77", "#88CCEE")) +
  labs(fill = 'Significance') +
  # Add text only to panels where fixed_effect_significance is significant (you might need to adjust based on how significance is coded)
  geom_text(aes(x = Inf, y = Inf, label = fixed_effect_significance), 
            hjust = 1.1, vjust = 1.1, 
            color = "black", size = 3, fontface = "bold")


g2 = ggplot(data = droplevels.data.frame(contrast_res_og[contrast_res_og$Indice %in% c('mass_int','mass_basal','mass_top'),]), 
            aes(x = absolute_value, fill = rd_slope_significance)) +
  # geom_density(stat = 'bin') +
  geom_vline(aes(xintercept = 0),
             color = "black", linewidth = 1, linetype = 3) + ## dotted line for vertical line in zero
  geom_histogram(stat = 'bin') +
  geom_vline(aes(xintercept = value_fixed_effect),
             color = "grey",, linewidth = 1.5) + ## full line for fixed effect slope (main effect across all PAs)
  geom_rect(aes(
    xmin = lower_2.5_fixed_effect,
    xmax = upper_97.5_fixed_effect,
    ymin = -Inf,
    ymax = Inf
  ), fill = 'grey', alpha = 0.01) +
  facet_wrap(~clean_Indice, scales = 'free', ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values= c("#DDCC77", "#88CCEE")) +
  xlab('Random slope - Site-level effects of protection') + 
  labs(fill = 'Significance difference from zero') +
  # Add text only to panels where fixed_effect_significance is significant (you might need to adjust based on how significance is coded)
  geom_text(aes(x = Inf, y = Inf, label = fixed_effect_significance), 
            hjust = 1.1, vjust = 1.1, 
            color = "black", size = 3, fontface = "bold")

g = ggpubr::ggarrange(g1,g2, ncol = 1, heights = c(0.75,0.25), common.legend = TRUE)

ggsave(g, filename = 'Plots/Figures-mixedModels/RandomSlopesDistributions-MainAnalysis.png',
       dpi = 1000, width = 8, height = 7)


## Table with fixed effect across all sites ----

to_save = unique(contrast_res[,c('clean_Indice','value_fixed_effect',
                                 'p_value_fixed_effect',
                                 'upper_97.5_fixed_effect',
                                 'lower_2.5_fixed_effect')])
to_save$significance = ifelse(to_save$p_value_fixed_effect < 0.001, '***', 
                              ifelse(to_save$p_value_fixed_effect < 0.01, '**',
                              ifelse(to_save$p_value_fixed_effect < 0.05, '*',
                              ifelse(to_save$p_value_fixed_effect < 0.1, '.', ''))))

names(to_save) = c('Food web metrics', 
                   'Fixed protection effect',
                   'P-value', 
                   'Upper 97.5 confidence interval',
                   'Lower 2.5 confidence interval',
                   'Significance code')

write.csv(to_save, 'Plots/Figures-mixedModels/FixedEffects-Protection.csv')





## Compute fragmentation for each PA ----

# We take the PA raster and for each PA on by one we use the 
# 'landscapemetrics' package to compute fragmentation measure 

europeRaster <- raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")

merged.PAs$ID.site = as.factor(merged.PAs$ID.new.Bioregion)
merged.PAs = droplevels(merged.PAs)
levels(merged.PAs$ID.site) = paste('site',1:46)
merged.PAs$ID.new.Bioregion.num = as.numeric(merged.PAs$ID.site)

res(europeRaster)
ncol(europeRaster) *  nrow(europeRaster)



### RUN IF NEED TO RECREATE THE PA RASTER, IF NOT: SKIP TO NEXT ----
# raster resolution

## create raster layers for PA and ID new Bioregion

# europeRaster <- raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")
# 
# europeRaster$PA = merged.PAs$PA[match(values(europeRaster),merged.PAs$Value)]
# europeRaster$ID.new.Bioregion.num = factor(merged.PAs$ID.new.Bioregion.num[match(values(europeRaster[[1]]),merged.PAs$Value)])
# europeRaster$LC = factor(merged.PAs$LC[match(values(europeRaster[[1]]),merged.PAs$Value)])
# 
# plot(europeRaster)
# 
# writeRaster(europeRaster, "outputs/WDPA.Grid/rasterWDPA.Over5yo.45PAs.NewPAID.LC.mergedPAs.tif", overwrite = T)

### Calculate PA fragmentation measure per site ----


library(landscapemetrics)
library(terra)
## created just above
PA.raster = rast("outputs/WDPA.Grid/rasterWDPA.Over5yo.45PAs.NewPAID.LC.mergedPAs.tif")
names(PA.raster) = c('Value','PA','ID.new.Bioregion.num','LC')
PA.raster$ID.new.Bioregion.num[which(is.na(terra::values(PA.raster$ID.new.Bioregion.num)))] = -10

# https://r-spatialecology.github.io/landscapemetrics/
landscapemetrics::check_landscape(subset(PA.raster, 2))
landscapemetrics::list_lsm(level = 'class', simplify = T)
landscapemetrics::list_lsm(level = 'landscape', simplify = T)

# landscapemetrics::list_lsm(level = 'class', type = 'cohesion')

## %: Approaches COHESION = 0 if patches of class i become more isolated. Increases if patches of class i become more aggregated.
plot(PA.raster[[2]])
lsm_c_cohesion(PA.raster[[2]], directions = 4)
## Equals 0 for maximally disaggregated and 100 for maximally aggregated classes.
## AI is an 'Aggregation metric'. It equals the number of like adjacencies divided 
## by the theoretical maximum possible number of like adjacencies for that class. 
## The metric is based on the adjacency matrix and the the single-count method.
lsm_c_ai(PA.raster[[2]])$value[1]


drivers$frag = NA
for(i in unique(drivers$ID.new.Bioregion)){
  
  ## extract the numerical ID for i
  ID.num = unique(merged.PAs$ID.new.Bioregion.num[which(merged.PAs$ID.new.Bioregion == i)])
  ## subset the raster (this only crops it to the extent of the PA + surroundings)
  PA.raster.sub = PA.raster[which(terra::values(PA.raster$ID.new.Bioregion.num) == ID.num), drop = F]
  ## subset to only grid cells within the PA by setting all others to NA
  terra::values(PA.raster.sub$PA)[which(terra::values(PA.raster.sub$ID.new.Bioregion.num) != ID.num)] = rep(NA, length(which(terra::values(PA.raster.sub$ID.new.Bioregion.num) != ID.num)))
  terra::values(PA.raster.sub$ID.new.Bioregion.num)[which(terra::values(PA.raster.sub$ID.new.Bioregion.num) != ID.num)] = rep(NA, length(which(terra::values(PA.raster.sub$ID.new.Bioregion.num) != ID.num)))
  
  ## extract fragmentation measure for protected grid cells
  drivers$frag[which(drivers$ID.new.Bioregion == i)] = lsm_c_ai(PA.raster.sub$PA)$value[2]
  
  plot(PA.raster.sub['PA'], main = paste(i, round(drivers$frag[which(drivers$ID.new.Bioregion == i)])))
}



# drivers$cat.sampled = ifelse(drivers$nb.inside > drivers$nb.outside, 'Inside','Outside')
# summary(factor(drivers$cat.sampled))

contrast_res = merge(contrast_res, drivers, by = 'ID.new.Bioregion', all.x = T)
# contrast_res = merge(contrast_res, sign, by = c('Indice','ID.new.Bioregion'))
contrast_res$ID.new.Bioregion = str_remove(contrast_res$ID.new.Bioregion, 'Bio-geographical ')
contrast_res$Bioregion = apply(contrast_res, 1, function(x) str_split(x['ID.new.Bioregion'], '_')[[1]][2])
# contrast_res$nb.Bioregion = paste('n =', contrast_res$nb.Bioregion)
# names(contrast_res) = str_replace_all(names(contrast_res), ' ', '')

# contrast_res$nb


## Computing the proportion of significant contrast ----


contrast_res$sign = ifelse(contrast_res$lower_2.5*contrast_res$upper_97.5 > 0 & contrast_res$value > 0, "Positive", 
                   ifelse(contrast_res$lower_2.5*contrast_res$upper_97.5 > 0 & contrast_res$value < 0, 'Negative',
                          'Overlapping zero'))



## Correlation of mfcl and nb of top species
cor(res$mfcl, res$S_top)




## Figure 1 ----


### Figure 1.A ----
# ## plotting the protected area map
bioregions <- sf::st_read("data/Galiana2021_network-area-europe-master/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp")
europeRaster <- raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")
bioregions = sf::st_transform(bioregions, crs = europeRaster@crs)

g_map = ggplot() + geom_sf(data = europe, fill = "white") +
  geom_sf(data = bioregions, fill = '#e6e6e6', alpha = 0.5) +
  geom_point(data = merged.PAs, aes(x = x, y = y, colour = factor(PA)), alpha = 0.30, size = 0.01, shape = 15) +
  geom_point(data = res, aes(x = x, y = y, colour = factor(PA)), size = 0.01, shape = 15) +
  scale_colour_manual(values = c('#FAA43A','forestgreen'), paste0("Protection","\n","status")) +
  xlim(2500000.1, 6200000.1) + ylim(1401766.1, 5231766.1) +
  guides(color = 'none', alpha = 'none')  +
  # guides(colour = guide_legend(override.aes = list(size = 2)), alpha = guide_legend(override.aes = list(size = 2))) +
  # labs(alpha = paste0("Included in","\n", "contrast","\n","analysis")) +
  ylab("") + xlab("") + 
  theme_bw(base_size = 20) +
  theme(panel.background = element_rect(fill = "azure1")) + 
  ggspatial::annotation_scale() + ggspatial::annotation_north_arrow(location = 'br')


### Side boxplots for Figure 1 -----
## proportion of agricultural area, forested areas,
## human density and protected areas directed at birds in protected vs non 
## protected communities across bioregions

merged.PAs$PA = as.factor(merged.PAs$PA)


mod.HD = lm(data = merged.PAs, log(m.HD) ~ Bioregion * PA)
x = visreg::visreg(mod.HD,'PA' ,by ='Bioregion', plot = F)

g1 = ggplot() + 
  geom_jitter(data = merged.PAs, aes(y = log(m.HD), x = paste(Bioregion, PA), colour = PA), alpha = 0.5) +
  # geom_line(data = x$fit, aes(x = paste(Bioregion,PA), y = visregFit)) +
  geom_crossbar(data = x$fit, aes(x = paste(Bioregion,PA), y = visregFit, ymin = visregLwr, ymax = visregUpr), fill = 'lightgrey', colour = 'darkgrey', alpha = 0.3) +
  geom_vline(xintercept = seq(2.5,8.5,2), colour = 'grey', linetype = 'longdash') +
  scale_colour_manual(values = c('#FAA43A','forestgreen'), paste0("Protection","\n","status")) + 
  scale_x_discrete(labels = rep('',10)) +
  ggtitle('Mean human density (logged)') + xlab('') + ylab('') + theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(size=10), legend.text = element_text(size = 10), plot.margin = unit(c(0,0,0,0),'lines'))

merged.PAs$lc1.Forestandseminaturalareas_plot = round(merged.PAs$lc1.Forestandseminaturalareas,4)
mod.forest = glm(data = merged.PAs, lc1.Forestandseminaturalareas_plot ~ Bioregion * PA, family = binomial)
x = visreg::visreg(mod.forest,'PA', by = 'Bioregion', scale = 'response')

g2 = ggplot() + 
  geom_jitter(data = merged.PAs, aes(y = lc1.Forestandseminaturalareas, x = paste(Bioregion, PA), colour = PA), alpha = 0.5) +
  # geom_line(data = x$fit, aes(x = paste(Bioregion,PA), y = visregFit)) +
  geom_crossbar(data = x$fit, aes(x = paste(Bioregion,PA), y = visregFit, ymin = visregLwr, ymax = visregUpr), fill = 'lightgrey', colour = 'darkgrey', alpha = 0.3) +
  geom_vline(xintercept = seq(2.5,8.5,2), colour = 'grey', linetype = 'longdash') +
  scale_colour_manual(values = c('#FAA43A','forestgreen'), paste0("Protection","\n","status")) + 
  scale_x_discrete(labels = rep('',10)) +
  ggtitle('Proportion of forested area') + xlab('') + ylab('') + theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(size=10), legend.text = element_text(size = 10), plot.margin = unit(c(0,0,0,0),'lines'))

merged.PAs$lc1.Agriculturalareas_plot = round(merged.PAs$lc1.Agriculturalareas,4)
mod.agri = glm(data = merged.PAs, lc1.Agriculturalareas_plot ~ Bioregion * PA, family = binomial)
x = visreg::visreg(mod.agri,'PA', by = 'Bioregion', scale = 'response')

g3 = ggplot() + 
  geom_jitter(data = merged.PAs, aes(y = lc1.Agriculturalareas, x = paste(Bioregion, PA), colour = PA), alpha = 0.5) +
  # geom_line(data = x$fit, aes(x = paste(Bioregion,PA), y = visregFit)) +
  geom_crossbar(data = x$fit, aes(x = paste(Bioregion,PA), y = visregFit, ymin = visregLwr, ymax = visregUpr), fill = 'lightgrey', colour = 'darkgrey', alpha = 0.3) +
  geom_vline(xintercept = seq(2.5,8.5,2), colour = 'grey', linetype = 'longdash') +
  scale_colour_manual(values = c('#FAA43A','forestgreen'), paste0("Protection","\n","status")) + 
  scale_x_discrete(labels = rep('',10)) +
  ggtitle('Proportion of agricultural area') + xlab('') + ylab('') + theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(size=10), legend.text = element_text(size = 10), plot.margin = unit(c(0,0,0,0),'lines'))

SE = broom::tidy(mod.agri)$std.error

X <- model.matrix(mod.agri)
p <- fitted(mod.agri) 
W <- diag(p*(1-p)) # matrix with diagonal == estimates
# this is the covariance matrix (inverse of Fisher information)
V <- solve(t(X)%*%W%*%X) # cov(x,y) = E[(x - E(x))*(y - E(y))]
all.equal(vcov(mod.agri), V)
#> [1] "Mean relative difference: 1.066523e-05"
# close enough

# these are the standard errors: take square root of diagonal 
all.equal(SE, sqrt(diag(V)))


mod.forBirds = glm(data = merged.PAs, Directed_at_birds ~ Bioregion, family = 'binomial')
x = visreg::visreg(mod.forBirds, 'Bioregion', scale = 'response', partial = F)
x$fit$PA = 1

g4 = ggplot() + 
  geom_jitter(data = merged.PAs, aes(y = Directed_at_birds, x = paste(Bioregion, PA), colour = PA), alpha = 0.5, height = 0.05) +
  # geom_line(data = x$fit, aes(x = paste(Bioregion,PA), y = visregFit)) +
  geom_crossbar(data = x$fit, aes(x = paste(Bioregion,PA), y = visregFit, ymin = visregLwr, ymax = visregUpr), width = 0.7, fill = 'lightgrey', colour = 'darkgrey', alpha = 0.3) +
  geom_vline(xintercept = seq(2.5,8.5,2), colour = 'grey', linetype = 'longdash') +
  scale_colour_manual(values = c('#FAA43A','forestgreen'), paste0("Protection","\n","status")) + 
  scale_x_discrete(labels = rep('',10)) +
  ggtitle('Protected area directed at birds') + xlab('') + ylab('') + theme_minimal(base_size = 10) + 
  scale_x_discrete(labels = c('', "Alpine Bio-Region", '',"Atlantic Bio-Region",'',
                                                             "Boreal Bio-Region", '', "Continental Bio-Region", '',
                                                             "Mediterranean Bio-Region")) +
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1), legend.text = element_text(size = 10), plot.margin = unit(c(0,0,0,0),'lines'), plot.title = element_text(size=10)) # change orientation of Bioregion labels and adjust plot margin so that there is not space between plots
  


### MAP for Figure 1 - Snowdonia site ----

library(sf)

PA.grid.nodups.shp = sf::st_read("outputs/WDPA.Shp/NewID.MergedPAs.1km.20012023.Over5yo.shp")
PA.grid.nodups.shp$ID.new.Bioregion = merged.PAs$ID.new.Bioregion[match(PA.grid.nodups.shp$Value,merged.PAs$Value)]
PA.grid.nodups.shp.eu = st_transform(PA.grid.nodups.shp, crs = crs(europeRaster))

PAID = "133_Atlantic Bio-geographical Region"

## cells for environmental chara
context.sf = st_as_sf(st_drop_geometry(merged.PAs[which(merged.PAs$ID.new.Bioregion == PAID),]), coords = c("x","y"), crs = crs(europeRaster))
context.sf = st_transform(context.sf, crs = crs(PA.grid.nodups.shp))

## cells for analysis
res.sf = st_as_sf(st_drop_geometry(res[which(res$ID.new.Bioregion == PAID),]), coords = c("x","y"), crs = crs(europeRaster))
res.sf = st_transform(res.sf, crs = crs(PA.grid.nodups.shp))

dt.sf = rbind(context.sf[,'Value'], res.sf[,'Value'])

# extract coordinates
res.sf$lat = st_coordinates(res.sf)[,'Y']
res.sf$lon = st_coordinates(res.sf)[,'X']

context.sf$lat = st_coordinates(context.sf)[,'Y']
context.sf$lon = st_coordinates(context.sf)[,'X']

# sub_WDPA = st_crop(PA.grid.nodups.shp,context.sf)

raster_sub = raster::crop(europeRaster, st_transform(dt.sf, crs(europeRaster)))
raster_to_pol_sub = st_as_sf(raster::rasterToPolygons(raster_sub))
raster_to_pol_sub_Context = merge(raster_to_pol_sub, st_drop_geometry(context.sf), by.x = 'PageName', by.y = 'Value', all.x = T)

raster_to_pol_sub_Analysis = merge(raster_to_pol_sub, st_drop_geometry(res.sf), by.x = 'PageName', by.y = 'Value')

sub_WDPA = PA.grid.nodups.shp.eu[which(PA.grid.nodups.shp.eu$ID.new.Bioregion == PAID),]
europe_sub = st_crop(europe, st_transform(raster_to_pol_sub, crs(europe)))

map = ggplot() + 
  geom_sf(data = europe_sub, fill = '#e6e6e6') + 
  geom_sf(data = sub_WDPA, fill = 'darkolivegreen3') +
  geom_sf(data = raster_to_pol_sub_Context, fill = NA) + 
  geom_sf(data = raster_to_pol_sub_Analysis, aes(fill = factor(PA)), shape = 15, size = 5, alpha = 0.6) +
  scale_fill_discrete(type = list(c('#FAA43A','forestgreen'))) +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = "azure1")) + xlab('') + ylab('') + 
  ggspatial::annotation_scale() 

# ggsave(plot = map, paste0('Plots/Figures/Figure1_map.png'), dpi = 200,scale = 1.5)


Fig1 = ggarrange(ggarrange(g_map, map, ncol = 1, heights = c(1, 0.8), widths = c(1,0.8), labels = c('A', 'B')) ,
                 ggarrange(g1,g2,g3,g4, heights = c(1,1,1,1.2), align = 'v', nrow = 4, common.legend = T, labels = c('C', 'D', 'E', 'F')), ncol = 2)

ggsave(plot = Fig1, "Plots/Figures-mixedModels/Figure1.png", width = 10, height = 8, dpi = 1000)
ggsave(plot = Fig1, "Plots/Figures-mixedModels/Figure1.pdf", width = 10, height = 8)


# Supp - All maps of PA networks ----
# for supplementary

for (PAID in unique(res$ID.new.Bioregion)){
  
  ## cells for environmental chara
  context.sf = st_as_sf(st_drop_geometry(merged.PAs[which(merged.PAs$ID.new.Bioregion == PAID),]), coords = c("x","y"), crs = crs(europeRaster))
  context.sf = st_transform(context.sf, crs = crs(PA.grid.nodups.shp))
  
  ## cells for analysis
  res.sf = st_as_sf(st_drop_geometry(res[which(res$ID.new.Bioregion == PAID),]), coords = c("x","y"), crs = crs(europeRaster))
  res.sf = st_transform(res.sf, crs = crs(PA.grid.nodups.shp))
  
  dt.sf = rbind(context.sf[,'Value'], res.sf[,'Value'])
  
  # sub_WDPA = st_crop(PA.grid.nodups.shp,context.sf)
  
  raster_sub = raster::crop(europeRaster, st_transform(dt.sf, crs(europeRaster)))
  raster_to_pol_sub = st_as_sf(raster::rasterToPolygons(raster_sub))
  raster_to_pol_sub_Context = merge(raster_to_pol_sub, st_drop_geometry(context.sf), by.x = 'PageName', by.y = 'Value', all.x = T)
  
  raster_to_pol_sub_Analysis = merge(raster_to_pol_sub, st_drop_geometry(res.sf), by.x = 'PageName', by.y = 'Value')
  
  sub_WDPA = PA.grid.nodups.shp.eu[which(PA.grid.nodups.shp.eu$ID.new.Bioregion == PAID),]
  europe_sub = st_crop(europe, st_transform(raster_to_pol_sub, crs(europe)))
  
  map = ggplot() + geom_sf(data = europe_sub, fill = 'white') + geom_sf(data = sub_WDPA, fill = 'darkolivegreen3') +
    geom_sf(data = raster_to_pol_sub_Context, fill = NA) + 
    geom_sf(data = raster_to_pol_sub_Analysis, aes(fill = factor(PA)), shape = 15, size = 5, alpha = 0.6) +
    scale_fill_discrete(type = list(c('#FAA43A','forestgreen'))) +
    ggtitle(PAID) + labs(colour = 'Protection status') + 
    theme(panel.background = element_rect(fill = "azure1")) + xlab('') + ylab('') + 
    ggspatial::annotation_scale()
  
  ggsave(plot = map, paste0('Plots/Figures/',PAID,'_map.png'), dpi = 200,scale = 1.5)
  
  ## version with google maps
  # p <- ggmap_bbox(p)
  # 
  # g1 = p +
  #   geom_point(data = context.sf, aes(x = lon, y = lat, colour = factor(PA), alpha = 0.3), shape = 15) +
  #   geom_point(data = res.sf, aes(x = lon, y = lat, colour = factor(PA)), shape = 15) +
  #   ggtitle(PAID) +
  #   geom_sf(data = PA.grid.nodups[which(PA.grid.nodups$ID.new.Bioregion == PAID),])
  # 
  # map = ggmap::get_map(location = unname(st_bbox(context.sf)), crop = F, source = 'google')
  # plot(st_transform(sub_WDPA,3857)[1], bgMap = map)
  # 
  # 
  # p <- ggmap::ggmap(ggmap::get_googlemap(center = c(lon = mean(extent(context.sf)[1:2]), 
  #                                                   lat = mean(extent(context.sf)[3:4])),
  #                                        zoom = 6, scale = 2,
  #                                        maptype ='satellite',
  #                                        color = 'color'))
  
}







# IDENTIFY DRIVERS  ----

## Final model for drivers of contrast in food web metrics ----


to_save = contrast_res[,c("ID.new.Bioregion","Indice",'Coef','Coef.S',"x.mean" ,"y.mean",'prop.area.protected',
                          intersect(str_replace(drivers.list, 'log.',''), names(contrast_res)),'nb.PA',
                          'remoteness.mean', 'HD.m.inside', 'HD.m.outside', 'slope.mean')]
to_save = merge(to_save, sign[,c('ID.new.Bioregion','Indice','sign')], by = c('ID.new.Bioregion','Indice'))

write.csv(to_save, 'outputs/2.Drivers/MainContrastTable_45PAs.csv')

### Version with scaled variables - BIC as selection criteria ----

bestmodels.BIC_scaled = list()
BIC.step_scaled = data.frame()
coef.combined_scaled = data.frame()

for(i in Indice.interest){
  cat(paste(c('-----------------------', i,'-----------------------'), collapse = "\n"))  
  
  data = subset(contrast_res, Indice == i)  
  
  # data$year.sd.std = data$year.sd/data$year.mean
  data$log.area.forest.inside = log(data$area.forest.inside)
  data$log.area.forest.outside = log(data$area.forest.outside)
  data$log.area.sum = log(data$area.sum)
  data$log.nb = log(data$nb.all)
  data$log.remoteness.m = log(data$remoteness.mean)
  data$log.slope = log(data$slope.mean)
  # data$log.N = log(data$N)
  data$log.HD.inside = log(data$HD.m.inside)
  data$log.HD.outside = log(data$HD.m.outside)
  data$log.HD.inside = log(data$HD.m.inside)
  data$log.area.mean = log(data$area.mean)
  
  data[,c(drivers.list,"prop.area.protected")] = apply(data[,c(drivers.list,"prop.area.protected")], 2, scale)
  
  # data$Coef.LM = data$mean.coef.PA
  data$Coef.LM = data$absolute_value
  data = data[,c('Coef.LM',drivers.list,"prop.area.protected")] 
  
  #### BIC STEP SELECTION
  
  print('BIC step selection --------------------------------------------------')
  
  null = lm(Coef.LM ~ 1 + prop.area.protected, data = data)
  ## doesn't like the formula function for some reason
  ## data was scaled above
  full = lm(Coef.LM ~ prop.area.protected + 
              
              # PA characteristics
              year.mean + prop.forBirds +
              frag + 
              score.m + log.area.mean +
              
              # context variables
              log.remoteness.m +
              area.urban.hab.outside + area.urban.hab.diff +
              area.agricultural.hab.outside + area.agricultural.hab.diff +
              log.HD.outside + HD.diff +
              # gHP.diff + gHP.m.inside + gHP.m.outside + 
              
              # natural landcover types
              log.area.forest.inside + area.forest.diff + log.area.forest.outside +
              diversity.LC +
              log.slope + elev.mean, data = data)
  
  selectionBIC = MASS::stepAIC(null, scope = list(lower = null, upper = full), data = data, direction = 'both',k=log(nrow(data)), trace = F)
  selectionAIC = MASS::stepAIC(null, scope = list(lower = null, upper = full), data = data, direction = 'both', trace = F)
  
  # print(selection$anova)
  # print(summary(selection))
  # simulationOutput <- simulateResiduals(fittedModel = selection, plot = T)
  
  ## get AIC
  nullAIC = lm(Coef.LM ~ 1,  data = data)
  
  ## get AIC
  tLL <- deviance(nullAIC) - 2*logLik(selectionBIC) 
  k <- selectionBIC$df
  
  BIC.step_scaled = rbind(BIC.step_scaled, data.frame(Indice = i, AIC = tLL+2*k, AIC.f = AIC(selectionBIC), mod = 'BIC step selection'))
  
  coefs = data.frame(c = selectionBIC$coefficients, se = summary(selectionBIC)$coefficients[,"Std. Error"], 
                     pvalue = summary(selectionBIC)$coefficients[,"Pr(>|t|)"], 
                     name = names(selectionBIC$coefficients), adjr2 = summary(selectionBIC)$adj.r.squared,
                     r2 = summary(selectionBIC)$r.squared, df = summary(selectionBIC)$df[2], Indice = i,
                     selection = "BIC")  
  
  coefs$name = str_remove_all(names(selectionBIC$coefficients), c('scale\\('))
  coefs$name = str_remove_all(coefs$name, c('\\)'))
  coefs$name = str_remove_all(coefs$name, c('log\\('))
  # coefs = coefs[-which(coefs$name %in% c('y.mean','x.mean')),]
  coefs = coefs[str_which(coefs$name, '\\:', negate = T),]
  
  coef.combined_scaled = rbind(coef.combined_scaled, coefs)
  
  bestmodels.BIC_scaled = append(bestmodels.BIC_scaled, list(selectionBIC))
  
  
  
  # if(mod[model] == "Passeriformes"){
  #   par(mfrow = c(3,3))
  # } else if (mod[model] == "Northern"){
  #   par(mfrow = c(3,2))
  # } else {par(mfrow = c(2,2))}
  
  par(mfrow = c(2,2))
  
  names.all = names(selectionBIC$coefficients)[which(!names(selectionBIC$coefficients) %in% c("(Intercept)","log.N"))]
  if(length(names.all)>0){
    # pdf(file = paste0('Plots/Figures/visreg/',i,'.pdf'), width = 10, height = 8)
    for(co in 1:length(names.all)){
      
      nn = names.all[co]
      print(nn)
      
      visreg::visreg(selectionBIC, nn, xlab = nn, ylab = unlist(i), cex.lab = 1.5)
      
    }
    # dev.off()
  }
  
} ### END FOR LOOP

coef.combined_scaled$Indice = unlist(clean_Indices[coef.combined_scaled$Indice])
coef.combined_scaled$name = unlist(clean_Drivers[coef.combined_scaled$name])
coef.combined_scaled$Indice = factor(coef.combined_scaled$Indice, levels = unique(coef.combined_scaled$Indice[order(coef.combined_scaled$adjr2)]))
# coef.combined_scaled$name = factor(coef.combined_scaled$name, levels = c("Proportion of protected areas managed for birds",
#                                                                          "Difference in proportion of agricultural area",
#                                                                          "Difference in human density", 
#                                                                          "Human density outside (logged)",
#                                                                          "Human density inside (logged)",                                                          "Proportion of urban habitat outside",
#                                                                          "Protected area continuity",
#                                                                          "Mean remoteness of the site (logged)" ,
#                                                                          "Diversity of land cover types",
#                                                                          "Proportion of forested area inside (logged)",
#                                                                          "Proportion of forested area outside (logged)",
#                                                                          "Proportion of protected area per site",             
#                                                                          "Intercept"))

write.csv(coef.combined_scaled,"outputs/MixedModels/2.Drivers/outputs-lm-scaled-MixedModels.csv")

coef.combined_scaled$sign = ifelse(coef.combined_scaled$pvalue<=0.001,"***", 
                                   ifelse(coef.combined_scaled$pvalue<=0.01,"**", 
                                          ifelse(coef.combined_scaled$pvalue<=0.05,"*", 
                                                 ifelse(coef.combined_scaled$pvalue<=0.1,".",""))))
coef.combined_scaled$sign = paste0(round(coef.combined_scaled$c,3),"\n",coef.combined_scaled$sign)


ggplot(coef.combined_scaled[-which(coef.combined_scaled$name %in% c("Intercept","Number of communities (logged)",
                                   "(Intercept","log.N","prop.area.protected","Proportion of protected area per site")),], 
       aes(x = name, y = Indice, fill = cut(c, c(-Inf,0,Inf)), label = sign)) + 
  geom_tile(colour = "white") + stat_bin2d(bins = 20) + geom_text(size = 8, vjust = 0.6) +
  scale_fill_manual(values = c("#df9faa","#cbe9f8"), "Coefficient") + 
  theme_bw(base_size = 28) + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  geom_label(aes(y = Indice, x = 0, label = round(adjr2, 2)), size = 8, fill = "white") + expand_limits(x= c(-0.5, 4)) +
  xlab("Drivers") + ylab("Mean contrast in food web metrics") +
  ggtitle("adjusted R²")

ggsave("Plots/Figures-MixedModels/BICstepSelection_ControllingForPropProtected_withmass_mixedModels.png", width = 20, height = 13)


ggplot(coef.combined_scaled, aes(x = name, y = Indice, fill = cut(c, c(-Inf,0,Inf)), label = sign)) + 
  geom_tile(colour = "white") + stat_bin2d(bins = 20) + geom_text(size = 8, vjust = 0.6) +
  scale_fill_manual(values = c("#df9faa","#cbe9f8"), "Coefficient") + 
  theme_bw(base_size = 28) + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  geom_label(aes(y = Indice, x = 0, label = round(adjr2, 2)), size = 8, fill = "white") + expand_limits(x= c(-0.5, 4)) +
  xlab("Drivers") + ylab("Mean contrast in food web metrics") +
  ggtitle("adjusted R²")

ggsave("Plots/Figures-MixedModels/BICstepSelection_ControllingForPropProtected_withmass_WITHCONTROLS.png", width = 22, height = 20)


ggplot(coef.combined_scaled, 
       aes(x = name, y = Indice, fill = cut(c, c(-Inf,0,Inf)), label = sign)) + 
  geom_tile(colour = "white") + stat_bin2d(bins = 20) + geom_text(size = 8, vjust = 0.6) +
  scale_fill_manual(values = c("#df9faa","#cbe9f8"), "Coefficient") + 
  theme_bw(base_size = 28) + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  geom_label(aes(y = Indice, x = 0, label = round(adjr2, 2)), size = 8, fill = "white") + expand_limits(x= c(-0.5, 4)) +
  xlab("Controls") + ylab("Mean contrast in food web metrics") +
  ggtitle("adjusted R?")

ggsave("Plots/Figures-MixedModels/BICstepSelection_CONTROLS_ControllingForPropProtected.png", width = 20, height = 20)


### Non-scaled version ----

## save plots
bestmodels.BIC_non_scaled = list()
BIC.step_non_scaled = data.frame()
coef.combined_non_scaled = data.frame()
plots_save = list()

for(i in Indice.interest){
  cat(paste(c('-----------------------', i,'-----------------------'), collapse = "\n"))  
  
  data = subset(contrast_res, Indice == i)  
  
  # data$year.sd.std = data$year.sd/data$year.mean
  data$log.area.forest.inside = log(data$area.forest.inside)
  data$log.area.forest.outside = log(data$area.forest.outside)
  data$log.area.sum = log(data$area.sum)
  data$log.nb = log(data$nb.all)
  data$log.remoteness.m = log(data$remoteness.mean)
  data$log.slope = log(data$slope.mean)
  #data$log.N = log(data$N)
  data$log.HD.inside = log(data$HD.m.inside)
  data$log.HD.outside = log(data$HD.m.outside)
  data$log.HD.inside = log(data$HD.m.inside)
  data$log.area.mean = log(data$area.mean)
  
  
  # data$Coef.LM = data$mean.coef.PA
  data$Coef.LM = data$absolute_value
  data = data[,c('Coef.LM',drivers.list,"prop.area.protected",'ID.new.Bioregion')] 
  
  
  #### BIC STEP SELECTION -----
  
  print('BIC step selection --------------------------------------------------')
  
  null = lm(Coef.LM ~ 1 + prop.area.protected, data = data)
  ## doesn't like the formula function for some reason
  ## data was scaled above
  full = lm(Coef.LM ~ prop.area.protected + 
              
              # PA characteristics
              year.mean + prop.forBirds +
              frag + 
              score.m + log.area.mean +
              
              # context variables
              log.remoteness.m +
              area.urban.hab.outside + area.urban.hab.diff +
              area.agricultural.hab.outside + area.agricultural.hab.diff +
              log.HD.outside + HD.diff +
              # gHP.diff + gHP.m.inside + gHP.m.outside + 
              
              # natural landcover types
              log.area.forest.inside + area.forest.diff + log.area.forest.outside +
              diversity.LC +
              log.slope + elev.mean, data = data)
  
  selectionBIC = MASS::stepAIC(null, scope = list(lower = null, upper = full), data = data, direction = 'both',k=log(nrow(data)), trace = F)
  
  # print(selection$anova)
  # print(summary(selection))
  # simulationOutput <- simulateResiduals(fittedModel = selection, plot = T)
  
  ## get AIC
  nullAIC = lm(Coef.LM ~ 1,  data = data)

  ## get AIC
  tLL <- deviance(nullAIC) - 2*logLik(selectionBIC) 
  k <- selectionBIC$df
  
  BIC.step_non_scaled = rbind(BIC.step_non_scaled, data.frame(Indice = i, AIC = tLL+2*k, AIC.f = AIC(selectionBIC), mod = 'BIC step selection'))
  
  coefs = data.frame(c = selectionBIC$coefficients, se = summary(selectionBIC)$coefficients[,"Std. Error"], 
                     pvalue = summary(selectionBIC)$coefficients[,"Pr(>|t|)"], 
                     name = names(selectionBIC$coefficients), adjr2 = summary(selectionBIC)$adj.r.squared,
                     r2 = summary(selectionBIC)$r.squared, df = summary(selectionBIC)$df[2], Indice = i,
                     selection = "BIC")  
  
  coefs$name = str_remove_all(names(selectionBIC$coefficients), c('scale\\('))
  coefs$name = str_remove_all(coefs$name, c('\\)'))
  coefs$name = str_remove_all(coefs$name, c('log\\('))
  # coefs = coefs[-which(coefs$name %in% c('y.mean','x.mean')),]
  coefs = coefs[str_which(coefs$name, '\\:', negate = T),]
  
  coef.combined_non_scaled = rbind(coef.combined_non_scaled, coefs)
  
  bestmodels.BIC_non_scaled = append(bestmodels.BIC_non_scaled, list(selectionBIC))
  
  par(mfrow = c(2,2))
  
  names.all = names(selectionBIC$coefficients)[which(!names(selectionBIC$coefficients) %in% c("(Intercept)", "prop.area.protected"))]
  if(length(names.all)>0){
    # pdf(file = paste0('Plots/Figures/visreg_NonScaled/',i,'.pdf'), width = 10, height = 8)
    # png(file = paste0('Plots/Figures/visreg_NonScaled/',i,'.png'), width = 800, height = 700, res = 100)
    # par(mar = c(4,4,4,4), mfrow = c(2,2))
    for(co in 1:length(names.all)){
      
      # if (length(names.all)>1){
      #   resid = sum(coefs[names.all,]$c * apply(data[,names.all],2,median))
      # } else {resid = coefs[names.all,]$c * median(data[,names.all])}
      
      var_interest = names.all[co] ## variable to plot
      print(var_interest)
      
      ## all variables but the variable to plot
      reduced_names = c(names.all[-co], "prop.area.protected")
      data_plot = data
      
      ## colour of the sign of the contrast
      colour = list('Negative' = "#CC6677", 'Positive' = "#88CCEE", 'Overlapping zero' = "#DDCC77")
      data_plot$col = unlist(colour[data_plot$sign])
      
      # newdata = data
      # if(length(reduced_names) > 1) {
      #   newdata[,reduced_names] = unlist(apply(data[,reduced_names], 2, function(x) rep(median(x), nrow(data))))
      # } else {
      #   newdata[,reduced_names] = rep(median(data[,reduced_names]), nrow(data))
      # }
      
      
      # data$resid = residuals(selectionBIC) + stats::predict(selectionBIC, newdata = newdata)
      
      newdata = data.frame(init = rep(1,100))
      newdata[,var_interest] = seq(min(data_plot[,var_interest]), max(data_plot[,var_interest]), length.out = 100)
      if(length(reduced_names) > 1) {
        newdata[,reduced_names] = unlist(apply(data_plot[,reduced_names], 2, function(x) rep(median(x), 100)))
      } else {
        newdata[,reduced_names] = rep(median(data_plot[,reduced_names]), 100)
      }
      
      # predict_CI1 = coefs['(Intercept)','c'] + as.matrix(newdata[,c(reduced_names, var_interest)]) %*% coefs[c(reduced_names, var_interest),'c']
      predict_CI = predict(selectionBIC, newdata, interval = "confidence") # get confidence interval
      intercept = coefs['(Intercept)',]$c + if(length(reduced_names) > 1) sum(apply(data_plot[,reduced_names],2,median)*coefs[reduced_names,]$c) else (median(data_plot[,reduced_names])*coefs[reduced_names,]$c) 
      
      # pred = coefs['(Intercept)',]$c + if(length(reduced_names) > 1) rowSums(data[,reduced_names]*coefs[reduced_names,]$c) else (data[,reduced_names]*coefs[reduced_names,]$c) ## one Marco thinks is best
      # pred = predict(lm(data = data, as.formula(paste0('Coef.LM ~', paste(reduced_names, collapse = '+')))))
      # data$resid = data$Coef.LM - pred
      
      g = ggplot() +
        geom_point(data = data_plot, aes_string(y = 'Coef.LM', x = var_interest), size = 2.5) +
        geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +
        # geom_smooth(data = data_plot, aes(y = resid, x = data_plot[,var_interest]), inherit.aes = F, method = 'lm', colour = 'grey') +
        xlab(clean_Drivers[var_interest][[1]]) + ylab(paste('Contrast in ','\n', clean_Indices[unlist(i)][[1]])) + 
        theme_bw(base_size = 18) + theme(plot.margin = unit(c(2, 3, 2, 0),'lines')) +
        # geom_abline(intercept = intercept, slope = coefs[var_interest,]$c) +
        geom_line(data = data.frame(predict_CI), aes_string(y = 'fit', x = newdata[,var_interest])) +
        geom_ribbon(data = data.frame(predict_CI), aes_string(ymin = 'lwr', ymax = 'upr', x = newdata[,var_interest]), alpha = 0.1)
      
      # print(g)
      
      plots_save[[paste(c(i, var_interest), collapse = '.')]] = g

    }
  }
  
} ### END FOR LOOP

names(bestmodels.BIC_non_scaled) = paste(Indice.interest, ', BIC non scaled')

names(bestmodels.BIC_non_scaled)
vars.order = c("log.N", drivers.list)
stargazer(bestmodels.BIC_non_scaled, type = 'html', column.labels = names(bestmodels.BIC_non_scaled), align = T, model.numbers = FALSE, intercept.bottom = FALSE, 
          omit.stat = c("rsq", "f"), 
          star.char = c("+", "*", "**", "***"), # change significance levels
          star.cutoffs = c(.1, .05, .01, .001), 
          order = paste0("^", vars.order , "$"))

coef.combined_non_scaled$Indice = unlist(clean_Indices[coef.combined_non_scaled$Indice])
coef.combined_non_scaled$name = unlist(clean_Drivers[coef.combined_non_scaled$name])


coef.combined_non_scaled$Indice = factor(coef.combined_non_scaled$Indice, levels = unique(coef.combined_non_scaled$Indice[order(coef.combined_non_scaled$adjr2)]))
coef.combined_non_scaled$name = factor(coef.combined_non_scaled$name, levels = c("Proportion of protected areas managed for birds",
                                                                         "Difference in proportion of agricultural area",
                                                                         "Difference in human density", 
                                                                         "Human density outside (logged)",
                                                                         "Human density inside (logged)",                                                          "Proportion of urban habitat outside",
                                                                         "Protected area continuity",
                                                                         "Mean remoteness of the site (logged)" ,
                                                                         "Diversity of land cover types",
                                                                         "Proportion of forested area inside (logged)",
                                                                         "Proportion of forested area outside (logged)",
                                                                         "Proportion of protected area per site",             
                                                                         "Intercept"))

write.csv(coef.combined_non_scaled,"outputs/MixedModels/2.Drivers/outputs-lm-non-scaled-MixedModels.csv")

coef.combined_non_scaled$sign = ifelse(coef.combined_non_scaled$pvalue<=0.001,"***", 
                                       ifelse(coef.combined_non_scaled$pvalue<=0.01,"**", 
                                              ifelse(coef.combined_non_scaled$pvalue<=0.05,"*", 
                                                     ifelse(coef.combined_non_scaled$pvalue<=0.1,".",""))))
coef.combined_non_scaled$sign = paste0(round(coef.combined_non_scaled$c,3),"\n",coef.combined_non_scaled$sign)


ggplot(coef.combined_non_scaled[-which(coef.combined_non_scaled$name %in% c("Intercept","Number of communities (logged)","Proportion of protected area per site")),], 
       aes(x = name, y = Indice, fill = cut(c, c(-Inf,0,Inf)), label = sign)) + 
  geom_tile(colour = "white") + stat_bin2d(bins = 20) + geom_text(size = 8, vjust = 0.6) +
  scale_fill_manual(values = c("#df9faa","#cbe9f8"), "Coefficient") + 
  theme_bw(base_size = 28) + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  geom_label(aes(y = Indice, x = 0, label = round(adjr2, 2)), size = 8, fill = "white") + expand_limits(x= c(-0.5, 4)) +
  xlab("Drivers") + ylab("Mean contrast in food web metrics") +
  ggtitle("adjusted R²")

# ggsave("Plots/Figures/BICstepSelection_NewEffort.png", width = 20, height = 10)
ggsave("Plots/Figures-MixedModels/BICstepSelection_NewEffort_NOTSCALED_final.png", width = 20, height = 12)


### FIGURE 4 - regression lines for a selection of environmental drivers ----

Fig4 = ggarrange(plots_save$mfcl.prop.forBirds, plots_save$mass_top.diversity.LC,
          plots_save$mfcl.elev.mean, plots_save$S_intermediate.elev.mean,
          plots_save$S_intermediate.frag, plots_save$C.area.urban.hab.diff,
          ncol = 3, nrow = 2, common.legend = T)
ggsave(plot = Fig4, filename = 'Plots/Figures-MixedModels/Figure4.png', width = 18, height = 10)
ggsave(plot = Fig4, filename = 'Plots/Figures-MixedModels/Figure4.pdf', width = 18, height = 10)


### Supp figs 1-3 - rest of environmental drivers for supplementary ----
SuppFig1 = ggarrange(plots_save$S.area.urban.hab.diff, plots_save$C.area.urban.hab.diff, 
                     plots_save$C.log.slope, plots_save$L.S.area.urban.hab.diff, 
                     plots_save$L.S.log.slope, plots_save$L.frag, 
                     ncol = 3, nrow = 2, common.legend = T)
ggsave(plot = SuppFig1, filename = 'Plots/Figures-MixedModels/SuppFigure1.png', width = 18, height = 10)
ggsave(plot = SuppFig1, filename = 'Plots/Figures-MixedModels/SuppFigure1.pdf', width = 18, height = 10)

SuppFig2 = ggarrange(plots_save$modularity.2.log.slope,
                     plots_save$modularity.2.HD.diff, plots_save$mfcl.elev.mean,
                     plots_save$mfcl.prop.forBirds, plots_save$mfcl.area.urban.hab.diff,
                     plots_save$Vul.prop.forBirds, 
                     ncol = 3, nrow = 2, common.legend = T)
ggsave(plot = SuppFig2, filename = 'Plots/Figures-MixedModels/SuppFigure2.png', width = 18, height = 10)
ggsave(plot = SuppFig2, filename = 'Plots/Figures-MixedModels/SuppFigure2.pdf', width = 18, height = 10)

SuppFig3 = ggarrange(plots_save$Vul.frag,
                     plots_save$mass_top.diversity.LC, plots_save$mass_int.area.urban.hab.diff,
                     plots_save$S_intermediate.frag, plots_save$S_intermediate.elev.mean,
                     plots_save$S_intermediate.area.agricultural.hab.outside, 
                     ncol = 3, nrow = 2, common.legend = T)
ggsave(plot = SuppFig3, filename = 'Plots/Figures-MixedModels/SuppFigure3.png', width = 18, height = 10)
ggsave(plot = SuppFig3, filename = 'Plots/Figures-MixedModels/SuppFigure3.pdf', width = 18, height = 10)

SuppFig4 = ggarrange(plots_save$S_top.elev.mean,
                     plots_save$S_top.HD.diff, plots_save$S_top.prop.forBirds,
                     ncol = 3, nrow = 2, common.legend = T)
ggsave(plot = SuppFig4, filename = 'Plots/Figures-MixedModels/SuppFigure4.png', width = 18, height = 10)
ggsave(plot = SuppFig4, filename = 'Plots/Figures-MixedModels/SuppFigure4.pdf', width = 18, height = 10)



# Plotting the food webs - code for 'grid' dataset from 2.gridCreationGBIF.EventList.Over5yo.R -----

setwd('F:/TheseSwansea/dta.Occurrence/GBIF/Europe/MetawebSpecies/Over5yo/Rarefying/EventLists')
eventList.min.sub = readRDS('EventListRarefied.Over5yo.1.rds')

library(igraph)
library(plyr)
metaweb_europe <- read.graph('F:/TheseSwansea/TraitStudy/code_Miguel/metaweb-europe.graphml', format = 'graphml')
m <- as_adjacency_matrix(metaweb_europe, attr = 'copresence', sparse=FALSE)
m[is.nan(m)] <- 0
source("F:/TheseSwansea/Galiana2021_network-area-europe-master/whois-function.R")

m_code = m # transform species names into species codes from TETRA EU
rownames(m_code) <- unlist((apply(as.matrix(rownames(m_code)), 1, function(x){ n <- whois(SPPNAME=x)[[1]]; if(length(n) == 0) {return (x)} else{return (n)} }  ) ))
colnames(m_code) <- unlist((apply(as.matrix(colnames(m_code)), 1, function(x){ n <- whois(SPPNAME=x)[[1]]; if(length(n) == 0) {return (x)} else{return (n)} }  ) ))

## we now have to (1) get the occurrence data corresponding to each event in eventList.min
## and (2) combine into a dataframe and (3) compute network metrics
setwd('F:/TheseSwansea/dta.Occurrence/GBIF/Europe/MetawebSpecies/Over5yo/Rarefying/OccData')

eventList.min.sub = subset(eventList.min.sub, class == 'Aves')
species.dt = data.frame()
nb = c()
for(f in str_subset(list.files(), 'rds')){
  df = readRDS(f)
  names(df)[which(names(df) == 'Value')] = 'df_grid'
  names(df)[which(names(df) == 'date')] = 'eventDate'
  if(!'datasetKey' %in% names(df)){
    df$datasetKey = 'eBird'
  }
  
  setDT(df)
  setDT(eventList.min.sub)
  df = merge(df[,c('masterCode','df_grid','eventDate','datasetKey')], 
             eventList.min.sub[,c('df_grid','eventDate','datasetKey','eventID')], 
             by = c('eventDate','datasetKey','df_grid'))
  
  
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

## pivot into grid
grid = tidyr::pivot_wider(species.dt[,c("df_grid","species","pres")], names_from = species, values_from = pres, values_fill = 0)
paste0('should be TRUE : ', nrow(grid) == length(unique(species.dt$df_grid)))
grid = merge(grid, cells_info[,c('x','y','PageName','Value')], by.x = 'df_grid', by.y = 'Value', all.x = T)
dim(grid)

### Food webs plotted as graphs ----------

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
community.p = cheddar::Community(nodes=data.frame(node=colnames(local_web.p)), trophic.links=cheddar::PredationMatrixToLinks(local_web.p), properties=list(title='Community'));
community.np = cheddar::Community(nodes=data.frame(node=colnames(local_web.np)), trophic.links=cheddar::PredationMatrixToLinks(local_web.np), properties=list(title='Community'));

# new local web that only keeps connected nodes
to_keep.p = setdiff(colnames(local_web.p),cheddar::IsolatedNodes(community.p)) # species that are not 'isolated' i.e. have at least 1 prey or 1 pred
to_keep.np = setdiff(colnames(local_web.np),cheddar::IsolatedNodes(community.np)) # species that are not 'isolated' i.e. have at least 1 prey or 1 pred

local_web.p = local_web.p[to_keep.p,to_keep.p]
local_web.np = local_web.np[to_keep.np,to_keep.np]

community.p = cheddar::Community(nodes=data.frame(node=colnames(local_web.p)), trophic.links=cheddar::PredationMatrixToLinks(local_web.p), properties=list(title='Community'));
community.np = cheddar::Community(nodes=data.frame(node=colnames(local_web.np)), trophic.links=cheddar::PredationMatrixToLinks(local_web.np), properties=list(title='Community'));

par(mfrow = c(1,2))
plot(community.p, show.nodes.as='labels', node.labels='node')
plot(community.np, show.nodes.as='labels', node.labels='node')

par(mfrow = c(1,2))
plot(community.p, main = 'protected network - 68 Mediterranean bioregion')
plot(community.np, main = 'non protected network - 68 Mediterranean bioregion')


library(visNetwork)
nodes = community.p$nodes; nodes$id = as.numeric(factor(nodes$node)); names(nodes) = c('label','id')
edges = community.p$trophic.links; names(edges) = c('from','to')
edges$from = nodes$id[match(edges$from,nodes$label)]
edges$to = nodes$id[match(edges$to,nodes$label,)]
edges = data.frame(edges, arrows = c("to"))

sub = sample(1:nrow(nodes), size = 30)
nodes = nodes[sub,]
edges = edges[which(edges$from %in% nodes$id & edges$to %in% nodes$id),]
visNetwork(nodes = nodes, edges = edges)

