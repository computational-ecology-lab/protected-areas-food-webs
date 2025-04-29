### Code for Thompson et al 2025
### "Mixed effect of protected areas on avian food webs"

## Functions to run the GAMMs models to extract PA effect on FW metrics 
## and to merge environmental variables with FW metrics at the grid cell level




## adding context variables to FW metric dataframe (called "res" throughout) ----

## this funtion uses the pre-extracted values for environmental variables 
## e.g. landcover, human density, elevation, protection status etc 
## and merges those info to the food web metrics dataframe (called 'res' throughout)

## it is used in all analysis scripts from 6. onwards.

## Environmental variables must be loaded into the environment 
## (via the .RData from script 3.ExtractEnvironmentalVariables.R)
add_context_var = function(res){
  
  res = merge(res, PA.grid.nodups[,c("WDPAID",'ID',"IUCN_CAT","STATUS_YR","area","Value",'dist','ISO3','score','Directed_at_birds')], by = 'Value', all.x = T)
  sum(res[which(!is.na(res$dist)),"dist.PA.min"]) == 0
  res[which(!is.na(res$dist)),"dist.PA.min"] = res[which(!is.na(res$dist)),"dist"]
  
  # add column of protection status
  res[is.na(res$IUCN_CAT),'IUCN_CAT'] = 0
  res$PA = ifelse(res$IUCN_CAT != 0, 1, 0)
  summary(factor(res$IUCN_CAT))
  
  
  for (x in list(Remoteness,CLC1,CLC2,CLC3,Elevation,Country,THPI,Slope,gHP,HD,temp,precipitation)){
    print(paste('Should be 0: ',anyDuplicated(x$Value)))
    res = merge(x, res, by = 'Value', all.y = T)
  }
  
  res$LC = apply(res[,str_subset(names(res),'lc\\.')], 1, function(x) ifelse(length(which.max(x))>0, str_subset(names(res),'lc\\.')[which.max(x)], NA))
  res$LC1 = apply(res[,str_subset(names(res),'lc1.')], 1, function(x) ifelse(length(which.max(x))>0, str_subset(names(res),'lc\\.')[which.max(x)], NA))
  res$LC2 = apply(res[,str_subset(names(res),'lc2.')], 1, function(x) ifelse(length(which.max(x))>0, str_subset(names(res),'lc\\.')[which.max(x)], NA))
  
  res$ID.new.Bioregion = paste0(res$NeighbouringPAID.new, '_', res$Bioregion)
  
  return(res)
}


## FUNCTIONS to run GAMMs to extract the difference in FW metrics between protected 
## and non protected grid cells

# this function returns three objects, 
# M.INDICE1 = GAMM without species richness
# M.INDICE1.1 = full GAMM
# M.INDICE1.2 = GAMM without fixed PA effect

function_GAMM <- function(data, Indice, fam = Family){
  
  # preparing the control variables (survey effort needs be logged) and all covariates are scaled
  data$log.effort_Aves = log(data$effort_Aves)
  data[,c('Remoteness.s','m.elev.s','log.effort_Aves.s')] = apply(data[,c('Remoteness','m.elev','log.effort_Aves')], 2, scale)
  
  data$ID.new.Bioregion = factor(data$ID.new.Bioregion)
  data$LC = factor(data$LC)
  
  if (Indice == "S") {
    if(length(levels(droplevels(factor(data$LC)))) > 1){
      M.INDICE1    <-gam(INDICE ~ 
                           ## fixed effects
                           PA + 
                           s(Remoteness.s, k = 4) + 
                           s(m.elev.s, k = 4) +
                           s(log.effort_Aves.s, k = 4) + 
                           
                           # random effects
                           s(Remoteness.s, ID.new.Bioregion, bs = 're') + # random slope
                           s(m.elev.s, ID.new.Bioregion, bs = 're') + # random slope
                           s(LC, bs = 're') + # random intercept
                           s(ID.new.Bioregion, bs = 're') + # random intercept
                           s(PA, ID.new.Bioregion, bs = 're'), # random slope - Main PA term that gives the effect of protection on INDICE
                         data = data, family=fam, method = 'REML')
    } else { 
      M.INDICE1    <-gam(INDICE ~ 
                           PA + 
                           s(Remoteness.s, k = 4) + 
                           s(m.elev.s, k = 4) + 
                           s(log.effort_Aves.s, k = 4) + 

                           s(Remoteness.s, ID.new.Bioregion, bs = 're') + 
                           s(m.elev.s, ID.new.Bioregion, bs = 're') + 
                           s(ID.new.Bioregion, bs = 're') + 
                           s(PA, ID.new.Bioregion, bs = 're'), 
                         data = data, family=fam, method = 'REML')
    }
    # gam.check(M.INDICE1)
    M.INDICE1.1 = M.INDICE1.2 = NULL
  } else {
    if(length(levels(droplevels(factor(data$LC)))) > 1){
      
      # model with no species richness
      M.INDICE1    <-gam(INDICE ~ PA + 
                           #s(S, k = 4) + 
                           s(Remoteness.s, k = 4) + 
                           s(m.elev.s, k = 4) + 
                           s(Remoteness.s, ID.new.Bioregion, bs = 're') + 
                           s(m.elev.s, ID.new.Bioregion, bs = 're') + 
                           s(log.effort_Aves.s, k = 4) + s(LC, bs = 're') + 
                           s(ID.new.Bioregion, bs = 're') + 
                           s(PA, ID.new.Bioregion, bs = 're'), 
                         data = data, family=fam, method = 'REML')
      
      # full model with species richness
      M.INDICE1.1  <-gam(INDICE ~ 
                           PA  +  
                           s(S, k = 4) + 
                           s(Remoteness.s, k = 4) + 
                           s(m.elev.s, k = 4) + 
                           s(Remoteness.s, ID.new.Bioregion, bs = 're') + 
                           s(m.elev.s, ID.new.Bioregion, bs = 're') + 
                           s(log.effort_Aves.s, k = 4) + s(LC, bs = 're') + 
                           s(ID.new.Bioregion, bs = 're') + 
                           s(PA, ID.new.Bioregion, bs = 're'), 
                         data = data, family=fam, method = 'REML')
      
      # model witouth PA fixed effect
      M.INDICE1.2  <-gam(INDICE ~   
                           # PA +
                           s(S, k = 4) + 
                           s(Remoteness.s, k = 4) + 
                           s(m.elev.s, k = 4) + 
                           s(Remoteness.s, ID.new.Bioregion, bs = 're') + 
                           s(m.elev.s, ID.new.Bioregion, bs = 're') + 
                           s(log.effort_Aves.s, k = 4) + 
                           s(LC, bs = 're') + 
                           s(ID.new.Bioregion, bs = 're') + 
                           s(PA, ID.new.Bioregion, bs = 're'), 
                         data = data, family=fam, method = 'REML')
    } else {
      M.INDICE1    <-gam(INDICE ~ PA +
                           # s(S, k = 4) + 
                           s(Remoteness.s, k = 4) + 
                           s(m.elev.s, k = 4) + 
                           s(Remoteness.s, ID.new.Bioregion, bs = 're') + 
                           s(m.elev.s, ID.new.Bioregion, bs = 're') + 
                           s(log.effort_Aves.s, k = 4) + 
                           # s(LC, bs = 're') + 
                           s(ID.new.Bioregion, bs = 're') + 
                           s(PA, ID.new.Bioregion, bs = 're'), 
                         data = data, family=fam, method = 'REML')
      
      M.INDICE1.1  <-gam(INDICE ~ PA +  
                           s(S, k = 4) + 
                           s(Remoteness.s, k = 4) + 
                           s(m.elev.s, k = 4) + 
                           s(Remoteness.s, ID.new.Bioregion, bs = 're') + 
                           s(m.elev.s, ID.new.Bioregion, bs = 're') + 
                           s(log.effort_Aves.s, k = 4) +
                           # s(LC, bs = 're') + 
                           s(ID.new.Bioregion, bs = 're') + 
                           s(PA, ID.new.Bioregion, bs = 're'), 
                         data = data, family=fam, method = 'REML')
      
      M.INDICE1.2  <-gam(INDICE ~ 
                           # PA + 
                           s(S, k = 4) + 
                           s(Remoteness.s, k = 4) + 
                           s(m.elev.s, k = 4) + 
                           s(Remoteness.s, ID.new.Bioregion, bs = 're') + 
                           s(m.elev.s, ID.new.Bioregion, bs = 're') + 
                           s(log.effort_Aves.s, k = 4) +
                           # s(LC, bs = 're') + 
                           s(ID.new.Bioregion, bs = 're') + 
                           s(PA, ID.new.Bioregion, bs = 're'), 
                         data = data, family=fam, method = 'REML')
    }
    # print(gam.check(M.INDICE1.1))
  }
  
  return(list(M.INDICE1, M.INDICE1.1, M.INDICE1.2))
  
}### end function


## start function 
function_GLMM <- function(data, Indice, fam = Family, PA = T){
  
  
  ## Run generalised linear mixed models to compare performance with GAMMs
  
  data$log.effort_Aves = log(data$effort_Aves)
  data[,c('Remoteness.s','m.elev.s','log.effort_Aves.s')] = apply(data[,c('Remoteness','m.elev','log.effort_Aves')], 2, scale)
  
  data$ID.new.Bioregion = factor(data$ID.new.Bioregion)
  data$LC = factor(data$LC)
  
  if (Indice == "S") {
    if(length(levels(droplevels(factor(data$LC)))) > 1){
      M.INDICE1  <- gam(INDICE ~ 
                          PA + 
                          s(Remoteness.s, k = 4) + 
                          s(m.elev.s, k = 4) +     
                          s(Remoteness.s, ID.new.Bioregion, bs = 're') + 
                          s(m.elev.s, ID.new.Bioregion, bs = 're') + 
                          log.effort_Aves.s + 
                          s(LC, bs = 're') + 
                          s(ID.new.Bioregion, bs = 're') + 
                          s(PA, ID.new.Bioregion, bs = 're'), 
                        data = data, family=fam, method = 'REML')
    } else {
      M.INDICE1  <- gam(INDICE ~ 
                          PA + 
                          s(Remoteness.s, k = 4) + 
                          s(m.elev.s, k = 4) +     
                          s(Remoteness.s, ID.new.Bioregion, bs = 're') + 
                          s(m.elev.s, ID.new.Bioregion, bs = 're') + 
                          log.effort_Aves.s + 
                          
                          s(ID.new.Bioregion, bs = 're') + 
                          s(PA, ID.new.Bioregion, bs = 're'), 
                        data = data, family=fam, method = 'REML')
    }
  }
  else {
    if(length(levels(droplevels(factor(data$LC)))) > 1){
      M.INDICE1  <- gam(INDICE ~ PA + s(Remoteness.s, k = 4) + s(m.elev.s, k = 4) + S + s(Remoteness.s, ID.new.Bioregion, bs = 're') + s(m.elev.s, ID.new.Bioregion, bs = 're') + log.effort_Aves.s + s(LC, bs = 're') + s(ID.new.Bioregion, bs = 're') + s(PA, ID.new.Bioregion, bs = 're'), data = data, family=fam, method = 'REML')
    } else {
      M.INDICE1  <- gam(INDICE ~ PA + s(Remoteness.s, k = 4) + s(m.elev.s, k = 4) + S + s(Remoteness.s, ID.new.Bioregion, bs = 're') + s(m.elev.s, ID.new.Bioregion, bs = 're') + log.effort_Aves.s                    + s(ID.new.Bioregion, bs = 're') + s(PA, ID.new.Bioregion, bs = 're'), data = data, family=fam, method = 'REML')
    }
  }
  return(M.INDICE1)
}### end function



## FUNCTIONS TO RUN GAMMS ----


# res resampled can be either a bootstrap of res (if we want equal
# survey design with same nb of protected and non protected grid cells - see function below)
# or the raw res dataframe

run_GAMMs = function(res_resampled, Indices.list){
  
  ranef = fixef = data.frame()
  
  for(IND in 1:length(Indices.list)){
    
    ranef_temp = fixef_temp = data.frame()
    
    Nom_indice<-Indices.list[IND]
    print(Nom_indice)
    res_resampled$INDICE = res_resampled[,Nom_indice]
    hist(res_resampled$INDICE, main = Nom_indice)
    
    if(Nom_indice %in% c("mfcl", "SDGen", "Gen", "SDVul", 
                         "Vul", "TL.mean","L.S",
                         "mass_top","mass_int","mass_basal"))                 { Family ="gaussian" 
    } else if (Nom_indice %in% c("S","S_top","S_basal","S_intermediate"))     { Family = "nb"
    } else if (Nom_indice %in% c("L"))                                        { Family = quasipoisson(link = "log") 
    } else                                                                    { Family = quasibinomial(link = "log")}
    
    print(Family)
    
    ## main GAM
    GAM = function_GAMM(data = res_resampled, fam = Family, Indice = Nom_indice)
    
    M.INDICE2 = GAM[[1]] ## model with PA term
    ranef2 = mixedup::extract_random_effects(M.INDICE2)
    ranef2$Model = 'No species richness'
    fixef2 = mixedup::extract_fixed_effects(M.INDICE2)
    fixef2$Model = 'No species richness'
    fixef2$devExplained = 1 - M.INDICE2$deviance/M.INDICE2$null.deviance
    fixef2$AIC = AIC(M.INDICE2)
    
    ranef_temp= rbind(ranef_temp, ranef2)
    fixef_temp = rbind(fixef_temp, fixef2)
    
    ## main GLM
    M.INDICE1 = function_GLMM(data = res_resampled, fam = Family, Indice = Nom_indice)
    ranefGLMM = mixedup::extract_random_effects(M.INDICE1)
    ranefGLMM$Model = 'GLMM'
    fixefGLMM = mixedup::extract_fixed_effects(M.INDICE1)
    fixefGLMM$Model = 'GLMM'
    fixefGLMM$devExplained = 1 - M.INDICE2$deviance/M.INDICE2$null.deviance
    fixefGLMM$AIC = AIC(M.INDICE2)
    
    ranef_temp= rbind(ranef_temp, ranefGLMM)
    fixef_temp = rbind(fixef_temp, fixefGLMM)

    if (Nom_indice != "S"){
      
      M.INDICE2.1 = GAM[[2]] ## PA + Species richness
      
      ranef2.1 = mixedup::extract_random_effects(M.INDICE2.1)
      ranef2.1$Model = 'Species richness'
      fixef2.1 = mixedup::extract_fixed_effects(M.INDICE2.1)
      fixef2.1$Model = 'Species richness'
      fixef2.1$devExplained = 1 - M.INDICE2.1$deviance/M.INDICE2.1$null.deviance
      fixef2.1$AIC = AIC(M.INDICE2.1)
      
      M.INDICE2.2 = GAM[[3]] ## NO PA, just species richness  
      
      ranef2.2 = mixedup::extract_random_effects(M.INDICE2.2)
      ranef2.2$Model = 'Species richness + no PA'
      fixef2.2 = mixedup::extract_fixed_effects(M.INDICE2.2)
      fixef2.2$Model = 'Species richness + no PA'
      fixef2.2$devExplained = 1 - M.INDICE2.2$deviance/M.INDICE2.2$null.deviance
      fixef2.2$AIC = AIC(M.INDICE2.2)
      
      ranef_temp = rbind(ranef_temp, ranef2.1)
      ranef_temp = rbind(ranef_temp, ranef2.2)
      
      fixef_temp = rbind(fixef_temp, fixef2.1)
      fixef_temp = rbind(fixef_temp, fixef2.2)
    }
    ranef_temp$Indice = Nom_indice
    fixef_temp$Indice = Nom_indice
    
    ranef = rbind(ranef, ranef_temp)
    fixef = rbind(fixef, fixef_temp)
 
  }
      
  return(list(ranef, fixef))
    
}### END OF FUNCTION
  

