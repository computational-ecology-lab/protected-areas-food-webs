### Code for "Mixed effect of protected areas on avian food webs"

## PLOTS WITH SUBSET CONTRAST RESULTS - changing overlap rules ----
## MIXED MODELS version
## 23/01/2025


# We changed overlap rules for selecting grid cells for contrast analysis
# 90% overlap = at least 90% of grid cells must overlap with a PA for them to be considered protected
# 50% overlap = at least 50% of grid cells must overlap with a PA for them to be considered protected
# in both cases non-protected grid cells must have 0% overlap

# in the og dataset, the criteria used was that the center of the grid cell
# has to overlap with the PA for them to be protected, otherwise they are considered unprotected



library(data.table)
library(stringr)
library(dplyr)
library(raster)
library(mgcv)
library(MASS)
# library(plyr)
library(corrplot)
library(DHARMa)
library(car)
library(lme4)
library(ggplot2)
library(ggpubr)
library(vegan)


setwd("E:/TheseSwansea/WDPA")
source('Script/FunctionsAnalysis_MixedModels_RandomSlopeControl.R') # source functions for analysis


## loading network metrics calculated at the grid cell level 
res = readRDS("outputs/1.Metrics/NetworkMetrics.Rarefied.nb15.effort50.Over5yo.AfterDesig.NewEffort.rds")
res = merge(res, readRDS("outputs/1.Metrics/NetworkMetrics.Rarefied.nb15.effort50.Over5yo.AfterDesig.NewEffort.BodyMass.rds"))

# communities with no top species or intermediate species
res$mass_top[which(is.na(res$mass_top))] = 0
res$mass_int[which(is.na(res$mass_int))] = 0

# metrics of interest for this analysis
Indice.interest = c("S","C","L.S", "L", "modularity.2", "mfcl", "omnivory.1", "Gen", 
                    "Vul", "TL.mean","mass_top", "mass_basal","mass_int",# "Fbasal","FInt","Ftop",
                    "S_basal","S_intermediate","S_top")


library(corrplot)
library(RColorBrewer)
M <- cor(res[,Indice.interest], use = 'complete.obs')
corrplot(M, type="upper", order="hclust", col=brewer.pal(n=8, name="RdYlBu"))
# corrplot(M[c("Gen","mfcl","S","FInt"),c("Gen","mfcl","S","FInt")], type="upper", order="hclust", col=brewer.pal(n=8, name="RdYlBu"))

### plotting the proportion of ebird vs gbif data
eventList.min.sub = readRDS('outputs/EventList/EventListMerged/Rarefying/EventListRarefied.Over5yo.NewEffort.1.rds')
eventList.min.sub = subset(eventList.min.sub, class == 'Aves')
eventList.min.sub = subset(eventList.min.sub, ID.new.Bioregion != "2_Continental Bio-geographical Region")

length(unique(eventList.min.sub$ID.new.Bioregion))




### CREATE A TABLE WITH ALL DRIVERS OF DIFFERENCES BETWEEN
### PROTECTED AND NON PROTECTED GRID CELLS

europeRaster <- raster::raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")
cells_info <- foreign::read.dbf('data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img.vat.dbf')
cells_info = cbind(cells_info, coordinates(europeRaster)[which(!is.na(raster::values(europeRaster))),])

europe = rnaturalearth::ne_countries(returnclass = "sf", scale = 'large')
europe = sf::st_transform(europe, crs(europeRaster))
europe = sf::st_crop(europe, extent(europeRaster))

# from 5.Analysis
load("outputs/ContextVariables.Over5yo.RData") # load all environmental variables

levels = data.frame(IUCN_CAT = levels(factor(PA.grid.nodups$IUCN_CAT)), score = c(10,9,8,7,6,1,1,1,5,4))
PA.grid.nodups = merge(PA.grid.nodups, levels, by = 'IUCN_CAT')
PA.grid.nodups$Directed_at_birds = 
  ifelse(PA.grid.nodups$DESIG_ENG %in% c("Ramsar Site, Wetland of International Importance", "Special Protection Area (Birds Directive)"),
         1, 0)

## CREATE A TABLE WITH ALL THE CARACTERISTICS OF EACH GRID PA
merged.PAs = readRDS("outputs/GridCells.NeighbouringPAID.new.Over5yo.RDS")
cells.to.seperate = c('555534657.1','555534704.1','555592567.1')
merged.PAs$NeighbouringPAID.new[which(merged.PAs$NeighbouringPAID %in% cells.to.seperate)] = '5.1'

## network metrics
res$NeighbouringPAID = merged.PAs$NeighbouringPAID[match(res$Value,merged.PAs$Value)]

## adding environmental variables to food web metrics ----
# Only grid cells with occurrence data (and network metrics):
res = add_context_var(res = res)

setDF(res)
res = res[,-which(str_detect(names(res), 'NA'))]
res = subset(res, STATUS_YR != 0 | is.na(STATUS_YR))
res = res[-which(res$S == 0),] # remove one community which had non interacting birds

ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = Remoteness), size = 0.8) + scale_colour_gradientn(colors = rainbow(5))

# ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = m.gHP), size = 0.8) + scale_colour_gradientn(colors = rainbow(5))
# ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = score), size = 0.8) + scale_colour_gradientn(colors = rainbow(5))

# ggplot() + geom_sf(data = europe) + geom_point(data = res, aes(x = x, y = y, colour = m.THPI))

## We calculate PA characteristics using all grid cells within 100 km of the PA 
## (not just those with sufficiently high survey effort)

## add env variables to the grid cells (all grid cells in a site, not just surveyed grid cells with bird occurrence data)
merged.PAs = add_context_var(res = merged.PAs[,c('Value','dist.PA.min','Bioregion','NeighbouringPAID.new','x','y')])

## remove missing values for the env variables
merged.PAs = subset(merged.PAs, !(is.na(LC) | LC == 'lc.NA'))
merged.PAs = subset(merged.PAs, !(is.na(m.elev)|is.na(m.slope)|is.na(m.gHP)|is.na(m.gHP)|is.na(Remoteness)|is.na(temp)))
merged.PAs = subset(merged.PAs, STATUS_YR != 0 | is.na(STATUS_YR))
merged.PAs = merged.PAs[,-which(str_detect(names(merged.PAs), 'NA'))]

## keep 100 km closest cells and only bioregions with sufficient data
merged.PAs = subset(merged.PAs, dist.PA.min <= 100e3 & ID.new.Bioregion %in% res$ID.new.Bioregion)
merged.PAs$IncludedInAnalysis = ifelse(merged.PAs$Value %in% res$Value, "Yes", "No")

removed = subset(merged.PAs, ID.new.Bioregion == "2_Continental Bio-geographical Region")
merged.PAs = subset(merged.PAs, ID.new.Bioregion != "2_Continental Bio-geographical Region")


length(unique(merged.PAs$ID.new.Bioregion))

## calculate the site - level metrics for environmental context analysis:

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
# agri = plyr::ddply(agri, c("name","PA"), summarise, value = mean(value))
# agri = plyr::ddply(agri, 'PA', summarise, cum = rev(cumsum(rev(value))), name = name, value = value)

# ggplot(data = agri, aes(y = value, x = PA, fill = name, label = name)) + 
#   # theme(axis.text.y=element_text(colour = colour)) +
#   geom_bar(stat = 'identity') + 
#   geom_text(data = agri, aes(y = (cum + cum - value)/2, x = PA, label = ifelse(value > 0.03,paste(name, round(value,2)), "")), size = 5) +
#   xlab("") + ylab("Proportion of land cover") 
# 
# ggsave("Plots/PropLC.png", width = 20, height = 10)


agri = tidyr::pivot_longer(merged.PAs[,c(agricultural_hab,"PA")], cols = all_of(agricultural_hab))
agri = plyr::ddply(agri, c("name","PA"), summarise, value = mean(value))
agri = plyr::ddply(agri, 'PA', summarise, cum = rev(cumsum(rev(value))), name = name, value = value)

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
merged.PAs$LC.diversity = diversity(merged.PAs[,str_subset(names(merged.PAs),'lc\\.')], MARGIN = 1)

names(merged.PAs) = str_replace_all(names(merged.PAs), ' ', '')

drivers = plyr::ddply(merged.PAs, "ID.new.Bioregion", summarise, NeighbouringPAID.new = unique(NeighbouringPAID.new), 
                 nb.inside.all = length(which(PA == 1)), nb.outside.all = length(which(PA == 0)), nb.all = unique(nb.inside.all + nb.outside.all), 
                 prop.area.protected = nb.inside.all/nb.all,
                 
                 
                 x.mean = mean(x), y.mean = mean(y), nb.PA = length(unique(ID[which(!is.na(ID))])),
                 ## context variables
                 
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
PA.chara = plyr::ddply(PA.chara, "ID.new.Bioregion", summarise, 
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

### Load overlap file (calculated in 6.AnalysisWDPA.BIRDS.FULLdatraset.Overlap) -----
overlap = read.csv("outputs/WDPA.Grid/CoverageFraction-PA-EuropeGrid.csv")
res = merge(res, overlap[,c('value','coverage_fraction')], by.x = 'Value', by.y = 'value', all.x = T)


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








## importing outputs from the GAMs (6.Analysis.WDPA) - difference in food web metrics ----

res_dir = "outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/90coverage/"
# contrast_res = readRDS(paste0(res_dir,"temps.08-08-2023_04h47.GAMoutputs.rds"))
# if everything ran in one go then the latest file is the only one needed 


contrast_res_90 = read.csv("outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/90coverage/GAMM-randomEffect-FullDataset-RandomSlope-FixedSmooths.csv")
contrast_res_90 = subset(contrast_res_90, effect == 'PA')
contrast_res_90 = subset(contrast_res_90, Model == 'Species richness' | (Model =='No species richness' & Indice == 'S'))
names(contrast_res_90)[names(contrast_res_90) == 'group'] = 'ID.new.Bioregion'

length(unique(contrast_res_90$ID.new.Bioregion)) # should be 16
plyr::ddply(contrast_res_90, "ID.new.Bioregion", summarise, length(Indice)) # should be 21 in each row


fixed_effect_90 = read.csv("outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/90coverage/GAMM-fixedEffect-FullDataset-RandomSlope-FixedSmooths.csv")
fixed_effect_90 = subset(fixed_effect_90, Model == 'Species richness' | (Model =='No species richness' & Indice == 'S'))
fixed_effect_90 = subset(fixed_effect_90, term == 'PA')

names(fixed_effect_90)[names(fixed_effect_90) == 'value'] = 'value_fixed_effect'
names(fixed_effect_90)[names(fixed_effect_90) == "p_value"] = 'p_value_fixed_effect'
names(fixed_effect_90)[names(fixed_effect_90) == "upper_97.5"] = "upper_97.5_fixed_effect"
names(fixed_effect_90)[names(fixed_effect_90) == "lower_2.5"] = "lower_2.5_fixed_effect"

contrast_res_90 = merge(fixed_effect_90[,c('value_fixed_effect','p_value_fixed_effect',
                                           'Indice','Model',"upper_97.5_fixed_effect","lower_2.5_fixed_effect")], 
                        contrast_res_90, by = c('Indice','Model'))

## random slope 
## se and CI are extracted from the covariance matrix
contrast_res_90$absolute_value = contrast_res_90$value_fixed_effect + contrast_res_90$value

contrast_res_90$rd_slope_significance = ifelse((contrast_res_90$lower_2.5) *
                                                 (contrast_res_90$upper_97.5) > 0, '***', '')



res_dir = "outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/50coverage/"
# contrast_res = readRDS(paste0(res_dir,"temps.08-08-2023_04h47.GAMoutputs.rds"))
# if everything ran in one go then the latest file is the only one needed 
contrast_res_50 = read.csv("outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/50coverage/GAMM-randomEffect-FullDataset-RandomSlope-FixedSmooths.csv")
contrast_res_50 = subset(contrast_res_50, effect == 'PA')
contrast_res_50 = subset(contrast_res_50, Model == 'Species richness' | (Model =='No species richness' & Indice == 'S'))
names(contrast_res_50)[names(contrast_res_50) == 'group'] = 'ID.new.Bioregion'

length(unique(contrast_res_50$ID.new.Bioregion))
plyr::ddply(contrast_res_50, "ID.new.Bioregion", summarise, length(Indice)) # should be 21 in each row

fixed_effect_50 = read.csv("outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/50coverage/GAMM-fixedEffect-FullDataset-RandomSlope-FixedSmooths.csv")
fixed_effect_50 = subset(fixed_effect_50, Model == 'Species richness' | (Model =='No species richness' & Indice == 'S'))
fixed_effect_50 = subset(fixed_effect_50, term == 'PA')

names(fixed_effect_50)[names(fixed_effect_50) == 'value'] = 'value_fixed_effect'
names(fixed_effect_50)[names(fixed_effect_50) == "p_value"] = 'p_value_fixed_effect'
names(fixed_effect_50)[names(fixed_effect_50) == "upper_97.5"] = "upper_97.5_fixed_effect"
names(fixed_effect_50)[names(fixed_effect_50) == "lower_2.5"] = "lower_2.5_fixed_effect"

contrast_res_50 = merge(fixed_effect_50[,c('value_fixed_effect','p_value_fixed_effect',
                                     'Indice','Model',"upper_97.5_fixed_effect","lower_2.5_fixed_effect")], 
                     contrast_res_50, by = c('Indice','Model'))

## random slope 
## se and CI are extracted from the covariance matrix
contrast_res_50$absolute_value = contrast_res_50$value_fixed_effect + contrast_res_50$value

contrast_res_50$rd_slope_significance = ifelse((contrast_res_50$lower_2.5) *
                                              (contrast_res_50$upper_97.5) > 0, '***', '')


res_dir = "outputs/1.Metrics/Results_differences_FW_Metrics/"
# contrast_res = readRDS(paste0(res_dir,"temps.08-08-2023_04h47.GAMoutputs.rds"))
contrast_res_og = read.csv('outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/randomEffect-FullDataset-RandomSlope.csv')
contrast_res_og = subset(contrast_res_og, effect == 'PA')
contrast_res_og = subset(contrast_res_og, Model == 'Species richness' | (Model =='No species richness' & Indice == 'S'))
names(contrast_res_og)[names(contrast_res_og) == 'group'] = 'ID.new.Bioregion'

length(unique(contrast_res_og$ID.new.Bioregion))
fixed_effect_res_og = read.csv('outputs/1.Metrics/Results_differences_FW_Metrics/NewEffort/fixedEffect-FullDataset-RandomSlope.csv')
fixed_effect_res_og = subset(fixed_effect_res_og, Model == 'Species richness' | (Model =='No species richness' & Indice == 'S'))
fixed_effect_res_og = subset(fixed_effect_res_og, term == 'PA')

names(fixed_effect_res_og)[names(fixed_effect_res_og) == 'value'] = 'value_fixed_effect'
names(fixed_effect_res_og)[names(fixed_effect_res_og) == "p_value"] = 'p_value_fixed_effect'
names(fixed_effect_res_og)[names(fixed_effect_res_og) == "upper_97.5"] = "upper_97.5_fixed_effect"
names(fixed_effect_res_og)[names(fixed_effect_res_og) == "lower_2.5"] = "lower_2.5_fixed_effect"

contrast_res_og = merge(fixed_effect_res_og[,c('value_fixed_effect','p_value_fixed_effect',
                                           'Indice','Model',"upper_97.5_fixed_effect","lower_2.5_fixed_effect")], 
                        contrast_res_og, by = c('Indice','Model'))

## random slope 
## se and CI are extracted from the covariance matrix
contrast_res_og$absolute_value = contrast_res_og$value_fixed_effect + contrast_res_og$value

contrast_res_og$rd_slope_significance = ifelse((contrast_res_og$lower_2.5) *
                                                 (contrast_res_og$upper_97.5) > 0, '***', '')

contrast_res_og$fixed_effect_significance = ifelse(contrast_res_og$p_value_fixed_effect < 0.001, 'p-value < 0.001', 
                                                   ifelse(contrast_res_og$p_value_fixed_effect < 0.01, 'p-value < 0.01',
                                                          ifelse(contrast_res_og$p_value_fixed_effect < 0.05, 'p-value < 0.05',
                                                                 ifelse(contrast_res_og$p_value_fixed_effect < 0.1, 'p-value < 0.1',''))))


# merge together the results with the overlap rule and the original results
contrast_res_90_comparison = merge(contrast_res_og, contrast_res_90, suffixes = c('_og','_new'), 
                        by = c('ID.new.Bioregion','Indice'), all.y = T)

contrast_res_50_comparison = merge(contrast_res_og, contrast_res_50, suffixes = c('_og','_new'), 
                                   by = c('ID.new.Bioregion','Indice'), all.y = T)


# look at how different they look:
plot(contrast_res_90_comparison$value_og, contrast_res_90_comparison$value_new)
plot(contrast_res_50_comparison$value_og, contrast_res_50_comparison$value_new)


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
length(unique(contrast_res$ID.new.Bioregion))

contrast_res_og$clean_Indice = unlist(clean_Indices[contrast_res_og$Indice])
contrast_res_50$clean_Indice = unlist(clean_Indices[contrast_res_50$Indice])
contrast_res_90$clean_Indice = unlist(clean_Indices[contrast_res_90$Indice])

contrast_res_og = subset(contrast_res_og, Indice %in% Indice.interest)
contrast_res_50 = subset(contrast_res_50, Indice %in% Indice.interest)
contrast_res_90 = subset(contrast_res_90, Indice %in% Indice.interest)



# Plotting random effects #####################################################


### plot of random effect as distributions with tick for each site

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


### Original dataset ----

           
write.csv(unique(contrast_res_og[,c("value_fixed_effect","p_value_fixed_effect","clean_Indice")]),
          "outputs/1.Metrics/FixedEffects-mainAnalysis.csv")



g1 = ggplot(data = droplevels.data.frame(contrast_res_og[!contrast_res_og$Indice %in% 
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

### 50 % overlap ----

unique(contrast_res_50[contrast_res_50$p_value_fixed_effect < 0.1,c("p_value_fixed_effect","Indice")])

write.csv(unique(contrast_res_50[,c("value_fixed_effect","p_value_fixed_effect","clean_Indice")]),
          "outputs/1.Metrics/FixedEffects-50percentOverlap.csv")

g1 = ggplot(data = droplevels.data.frame(contrast_res_50[!contrast_res_50$Indice %in% 
                                                           c('mass_int','mass_basal','mass_top'),]), 
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
  labs(fill = 'Significance')

g2 = ggplot(data = droplevels.data.frame(contrast_res_50[contrast_res_50$Indice %in% c('mass_int','mass_basal','mass_top'),]), 
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
  labs(fill = 'Significance difference from zero')

g = ggpubr::ggarrange(g1,g2, ncol = 1, heights = c(0.8,0.2), common.legend = TRUE)


ggsave(g, filename = 'Plots/Figures-mixedModels/Supp-RandomSlopesDistributions-50percent.png',
       dpi = 1000, width = 8, height = 7)


### 90 % overlap ----

write.csv(unique(contrast_res_90[,c("value_fixed_effect","p_value_fixed_effect","clean_Indice")]),
          "outputs/1.Metrics/FixedEffects-90percentOverlap.csv")

g1 = ggplot(data = droplevels.data.frame(contrast_res_90[!contrast_res_90$Indice %in% 
                                                           c('mass_int','mass_basal','mass_top'),]), 
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
  labs(fill = 'Significance')

g2 = ggplot(data = droplevels.data.frame(contrast_res_90[contrast_res_90$Indice %in% c('mass_int','mass_basal','mass_top'),]), 
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
  labs(fill = 'Significance difference from zero')

g = ggpubr::ggarrange(g1,g2, ncol = 1, heights = c(0.8,0.2), common.legend = TRUE)


ggsave(g, filename = 'Plots/Figures-mixedModels/Supp-RandomSlopesDistributions-90percent.png',
       dpi = 1000, width = 8, height = 7)



## Selection of dataset to plot (90% or 50%? New or original?) -----

#### new datasets ----

### UNCOMMENT THE DATASET YOU WANT TO RUN:


## 90% overlap
contrast_res = contrast_res_90
plot_names = '90percent_overlap'
plot_title = '90% overlap (N = 16 sites)'

## 50% overlap with 90% sites 

contrast_res = contrast_res_50[which(contrast_res_50$ID.new.Bioregion %in% unique(contrast_res_90$ID.new.Bioregion))]
plot_names = '50percent_overlap_subsetted_with_90percent_sites'
plot_title = '50% overlap for sites from 90% overlap dataset (N = 16 sites)'

## 50% overlap 

contrast_res = contrast_res_50
plot_names = '50percent_overlap'
plot_title = '50% overlap (N = 34 sites)'


#### old datasets ----

## sites of 90% overlap with original cells (middle crossing only - no overlap criteria)

contrast_res = subset(contrast_res_og, ID.new.Bioregion %in% unique(contrast_res_90$ID.new.Bioregion))
plot_names = '90percent_overlap_original'
plot_title = 'Original dataset with sites from 90% overlap (N = 16 sites)'

## sites of 50% overlap with original cells (middle crossing only - no overlap criteria)

contrast_res = subset(contrast_res_og, ID.new.Bioregion %in% unique(contrast_res_50$ID.new.Bioregion))
plot_names = '50percent_overlap_original'
plot_title = 'Original dataset with sites from 50% overlap (N = 34 sites)'







## ALL MAPS OF PROTECTED AREAS ----
# for supplementary

library(sf)
europeRaster <- raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")

PA.grid.nodups.shp = sf::st_read("outputs/WDPA.Shp/NewID.MergedPAs.1km.20012023.Over5yo.shp")
PA.grid.nodups.shp$ID.new.Bioregion = merged.PAs$ID.new.Bioregion[match(PA.grid.nodups.shp$Value,merged.PAs$Value)]
PA.grid.nodups.shp.eu = st_transform(PA.grid.nodups.shp, crs = crs(europeRaster))


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
  
  ggsave(plot = map, paste0('Plots/Figures/',plot_names,'/',PAID,'_map.png'), dpi = 200,scale = 1.5)
  
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







## Compute fragmentation for each PA: take the PA raster and subset each PA ----
## use the landscapemetrics package to compute fragmentation measure 

europeRaster <- raster(x="data/Galiana2021_network-area-europe-master/mask10k/reference_grid_10km.img")

merged.PAs$ID.site = as.factor(merged.PAs$ID.new.Bioregion)
merged.PAs = droplevels(merged.PAs)
levels(merged.PAs$ID.site) = paste('site',1:46)
merged.PAs$ID.new.Bioregion.num = as.numeric(merged.PAs$ID.site)

res(europeRaster)
ncol(europeRaster) *  nrow(europeRaster)


### RUN IF NEED TO RECREATE THE PA RASTER, IF NOT: SKIP TO NEXT -----
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



library(landscapemetrics)
library(terra)
## created just above
PA.raster = rast("outputs/WDPA.Grid/rasterWDPA.Over5yo.45PAs.NewPAID.LC.mergedPAs.tif")
names(PA.raster) = c('Value','PA','ID.new.Bioregion.num','LC')
PA.raster$ID.new.Bioregion.num[which(is.na(values(PA.raster$ID.new.Bioregion.num)))] = -10

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
## The metric is based on he adjacency matrix and the the single-count method.
lsm_c_ai(PA.raster[[2]])$value[1]


drivers$frag = NA
for(i in unique(drivers$ID.new.Bioregion)){
  
  ## extract the numerical ID for i
  ID.num = unique(merged.PAs$ID.new.Bioregion.num[which(merged.PAs$ID.new.Bioregion == i)])
  ## subset the raster (this only crops it to the extent of the PA + surroundings)
  PA.raster.sub = PA.raster[which(values(PA.raster$ID.new.Bioregion.num) == ID.num), drop = F]
  ## subset to only grid cells within the PA by setting all others to NA
  values(PA.raster.sub$PA)[which(values(PA.raster.sub$ID.new.Bioregion.num) != ID.num)] = rep(NA, length(which(values(PA.raster.sub$ID.new.Bioregion.num) != ID.num)))
  values(PA.raster.sub$ID.new.Bioregion.num)[which(values(PA.raster.sub$ID.new.Bioregion.num) != ID.num)] = rep(NA, length(which(values(PA.raster.sub$ID.new.Bioregion.num) != ID.num)))
  
  ## extract fragmentation measure for protected grid cells
  drivers$frag[which(drivers$ID.new.Bioregion == i)] = lsm_c_ai(PA.raster.sub$PA)$value[2]
  
  plot(PA.raster.sub['PA'], main = paste(i, round(drivers$frag[which(drivers$ID.new.Bioregion == i)])))
}







## Formating contrast_res for plots ----

drivers$cat.sampled = ifelse(drivers$nb.inside > drivers$nb.outside, 'Inside','Outside')
summary(factor(drivers$cat.sampled))

contrast_res = merge(contrast_res, drivers, by = 'ID.new.Bioregion', all.x = T)
# contrast_res = merge(contrast_res, sign, by = c('Indice','ID.new.Bioregion'))
contrast_res$ID.new.Bioregion = str_remove(contrast_res$ID.new.Bioregion, 'Bio-geographical ')
contrast_res$Bioregion = apply(contrast_res, 1, function(x) str_split(x['ID.new.Bioregion'], '_')[[1]][2])
contrast_res$nb.Bioregion = paste('n =', contrast_res$nb.Bioregion)
names(contrast_res) = str_replace_all(names(contrast_res), ' ', '')

# contrast_res$nb

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
length(unique(contrast_res$ID.new.Bioregion))




### Proportion of significant results ----

## running t test on bootstrapped coefficient for difference in metrics
sign = plyr::ddply(contrast_res[which(!contrast_res$Indice %in% c('S')),], c('Indice','ID.new.Bioregion'), summarize, 
                   mean.coef.PA = round(mean(unlist(Coef.distribution.S)),5), sd.coef.PA = round(sd(unlist(Coef.distribution.S)),3),
                   q075 = quantile(unlist(Coef.distribution.S), 0.075), q975 = quantile(unlist(Coef.distribution.S), 0.975),
                   dev.explained.S = mean(DevExplained.S),dev.explained = mean(DevExplained))

sign.S = plyr::ddply(contrast_res[which(contrast_res$Indice %in% c('S')),], c('Indice','ID.new.Bioregion'), summarise, 
                     mean.coef.PA = round(mean(unlist(Coef.distribution)),5), sd.coef.PA = round(sd(unlist(Coef.distribution)),3),
                     q075 = quantile(unlist(Coef.distribution), 0.075), q975 = quantile(unlist(Coef.distribution), 0.975),
                     dev.explained.S = mean(DevExplained),dev.explained = mean(DevExplained))

sign = rbind(sign, sign.S)

sign$sign = ifelse(sign$q075*sign$q975>0 & sign$mean.coef.PA > 0, "Positive", 
                   ifelse(sign$q075*sign$q975>0 & sign$mean.coef.PA < 0, 'Negative',
                          'Overlapping zero'))


## prop of significant positive and negative effects across sites x network metrics (see also figure S1)
# plyr::ddply(sign, c('Indice'), summarise, prop.sign = length(which(t.test<0.05))/length(t.test), 
#       pos = length(which(t.test<0.05 & mean.coef.PA > 0)), prop.pos = pos/length(t.test), 
#       neg = length(which(t.test<0.05 & mean.coef.PA < 0)), prop.neg = neg/length(t.test))

p = plyr::ddply(sign, c('Indice'), summarise, prop.sign = length(which(q075*q975>0))/length(t.test), 
          pos = length(which(q075*q975>0 & mean.coef.PA > 0)), prop.pos = pos/length(q075), mean.coef.pos = mean(mean.coef.PA[which(q075*q975>0 & mean.coef.PA > 0)]),
          neg = length(which(q075*q975>0 & mean.coef.PA < 0)), prop.neg = neg/length(q075), mean.coef.neg = mean(mean.coef.PA[which(q075*q975>0 & mean.coef.PA < 0)]),
          neutral = length(which(q075*q975<0)), prop.neutral = neutral/length(q075), mean.coef.neutral = mean(mean.coef.PA[which(q075*q975<0)]))

p[order(p$pos, decreasing = T),]
p[order(p$neg, decreasing = T),]

mean(p$prop.neutral); mean(p$neutral)
mean(p$prop.pos); mean(p$pos)
mean(p$prop.neg); mean(p$neg)

## plotting the direction + significance of the differences in network metrics across sites
# with proportions
l = tidyr::pivot_longer(p, cols = c('prop.neg','prop.pos','prop.neutral'))

l$Indice = factor(l$Indice, levels = unique(l$Indice)[order(l$value[which(l$name == "prop.pos")])]) # reorder Indices for plotting
l$name = factor(l$name, levels = c("prop.neutral",'prop.pos','prop.neg')) # reorder levels for plotting
levels(l$name) = unlist(list("prop.neutral" = "Overlapping zero",
                      'prop.neg' = "Negative", 
                      'prop.pos' = "Positive")[levels(l$name)]) # rename levels 

# same with mean coefficient (effect size)
l_mean = tidyr::pivot_longer(p, cols = c('mean.coef.pos','mean.coef.neg','mean.coef.neutral'))
l_mean$name = factor(l_mean$name, levels = c("mean.coef.neutral",'mean.coef.pos','mean.coef.neg')) # reorder levels for plotting
levels(l_mean$name) = unlist(list("mean.coef.neutral" = "Overlapping zero",
                                  'mean.coef.neg' = "Negative", 
                                  'mean.coef.pos' = "Positive")[levels(l_mean$name)]) # rename levels 

l = merge(l, l_mean, by = c('Indice','name'), suffixes = c("_prop","_mean"))

# reorder so that it is in the order - Positive - Overlapping zero - Negative
l = plyr::ddply(l, 'Indice', summarise, name_reordered = c("Positive","Overlapping zero","Negative"), 
           value_reordered = value_prop[match(c("Positive","Overlapping zero","Negative"), name)],
           mean_reordered = value_mean[match(c("Positive","Overlapping zero","Negative"), name)],
           cum = cumsum(value_reordered))

#l = plyr::ddply(l, 'Indice', summarise, cum = rev(cumsum(rev(value))), name = name, value = value)
l$clean_Indice = l$Indice
levels(l$clean_Indice) = unlist(clean_Indices[levels(l$Indice)])



## Figures 2 and 3 -- signs of the differences -----


safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")



### Figure 2 ----
img.over = png::readPNG("Plots/Figures/OverlapDifference.png")
img.over = grid::rasterGrob(img.over, interpolate=TRUE)

img.pos = png::readPNG("Plots/Figures/PositiveDifference.png")
img.pos = grid::rasterGrob(img.pos, interpolate=TRUE)

img.neg = png::readPNG("Plots/Figures/NegativeDifference.png")
img.neg = grid::rasterGrob(img.neg, interpolate=TRUE)

l_FW = l[which(l$Indice %in% Indice.interest[-c(11:13)]),]

g1 = ggplot(data = l_FW, aes(x = value_reordered, y = clean_Indice, fill = name_reordered)) + 
  # theme(axis.text.y=element_text(colour = colour)) +
  geom_bar(stat = 'identity', alpha = 0.7) + 
  geom_text(data = l_FW, aes(x = (cum + cum - value_reordered)/2, y = clean_Indice, label = round(mean_reordered,2)), size = 8) +
  ylab("") + xlab("Proportion of sites") + labs(fill = paste0("Sign of the contrast in food web metric")) +
  scale_fill_manual(values= c("#CC6677", "#DDCC77", "#88CCEE")) + theme_bw(base_size = 25) +
  # theme(plot.margin = unit(c(0,20,0,0), "cm"), legend.position = "top") +
  annotation_custom(img.over, xmin=0.20, xmax=0.35, ymin=-0.5, ymax=3) +
  annotation_custom(img.pos, xmin=0, xmax=0.15, ymin=11.1, ymax=13.4) +
  annotation_custom(img.neg, xmin=0.90, xmax=1.05, ymin=-0.5, ymax=3) +
#annotation_custom(ggplotGrob(g2), xmin=1, xmax=2.15, ymin=0, ymax=7)
  ggtitle(plot_title)


ggsave(plot = g1, paste0('Plots/Figures/test/Figure2-reordered',plot_names,'.pdf'), width = 20, height = 10)
ggsave(plot = g1, paste0('Plots/Figures/test/Figure2-reordered',plot_names,'.png'), width = 20, heigh = 10)

## with body mass 
l_mass = l[which(l$Indice %in% Indice.interest[c(11:13)]),]

g2 = ggplot(data = l_mass, aes(x = value_reordered, y = clean_Indice, fill = name_reordered)) + 
  # theme(axis.text.y=element_text(colour = colour)) +
  geom_bar(stat = 'identity', alpha = 0.7) + theme_bw(base_size = 25) + #theme(legend.position="none") +
  geom_text(data = l_mass, aes(x = (cum + cum - value_reordered)/2, y = clean_Indice, label = round(mean_reordered,2)), size = 8) +
  ylab("") + xlab("Proportion of sites") + labs(fill = paste0("Sign of the contrast in body mass")) + # ggtitle("Contrasts in body mass") +
  scale_fill_manual(values= c("#CC6677", "#DDCC77","#88CCEE")) +
  ggtitle(plot_title)

ggsave(plot = g2, paste0('Plots/Figures/test/Figure3-reordered',plot_names,'.png'), width = 20, height = 4)
ggsave(plot = g2, paste0('Plots/Figures/test/Figure3-reordered',plot_names,'.pdf'), width = 20, height = 4) # smaller height to that it matches Figure 2


