## metrics calculation


localMetrics = function(master_reduced, iterations, m_code, bodySize){
  
  ## PARALLEL LOOP
  library(foreach)
  library(doParallel)
  library(doSNOW) # for progression bar
  
  
  cores <- detectCores() - 2
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # this uses the nb of cores set in registerdopareallel
  cur.dat = foreach(cell = 1:iterations, .combine=rbind, .options.snow = opts) %dopar% { 
    # for(cell in 1:iterations){
    ## INITIALISATION
    
    L <- S <- C <- L.S <- Gen <- SDGen <- MeanGen <- Vul <- SDVul <- MeanVul <- Fbasal <- Ftop <- FInt <- Fother <- S_basal <- S_top <- S_intermediate <- overlap <- omnivory <- mfcl <- 
      modularity <- TL.mean <- isolated <- mass_int <- mass_basal <- mass_top <- NULL
    
    modularity.1 = modularity.2 = modules.1 = modules.2 = omnivory.1 = NULL
    
    library(stringr)
    library(igraph)
    library(cheddar)
    source("F:/TheseSwansea/Galiana2021_network-area-europe-master/utils.R")
    source("data/Galiana2021_network-area-europe-master/whois-function.R")
    
    
    pixel = master_reduced[cell,] ## local community (sp that occur within this 10x10km grid cell)
    sppSTRING <- names(pixel)[which(pixel==1 & str_detect(names(pixel), pattern = "[:digit:]"), arr.ind=FALSE)]
    
    
    ## CREATION OF LOCAL WEB
    if (length(sppSTRING)>1){
      
      local_web = m_code[sppSTRING,sppSTRING]
      sum(diag(local_web)) == 0
      
      # creating a cheddar community to look for isolated nodes
      community = cheddar::Community(nodes=data.frame(node=colnames(local_web)), trophic.links=PredationMatrixToLinks(local_web), properties=list(title='Community'));
      
      # new local web that only keeps connected nodes
      to_keep = setdiff(colnames(local_web),IsolatedNodes(community)) # species that are not 'isolated' i.e. have at least 1 prey or 1 pred
      
      #### networks with no links: 
      if (length(to_keep)==0){ 
        
        S <- 0 # nb of sp
        S.total <- length(sppSTRING) # nb of species not necessarily involved in the bird food web
        
        L <- C <- L.S <- Gen <- SDGen <- MeanGen <- Vul <- SDVul <- MeanVul <- Fbasal <- Ftop <- FInt <- Fother <- S_basal <- S_top <- S_intermediate <- overlap <- omnivory <- mfcl <- 
          modularity <- TL.mean <- mass_int <- mass_basal <- mass_top <-  NA
        
        modularity.1 = modularity.2 = modules.2 = omnivory.1 = NA
        
        isolated = T
        
        test = cbind(pixel[c("PageName")],cbind(L,S,S.total,C,L.S,Gen,SDGen,Vul,SDVul,Fbasal,Ftop,FInt,Fother,S_basal,S_top,S_intermediate,overlap,omnivory,mfcl,
                                                modularity,TL.mean, isolated, modularity.2, modules.1, modules.2, omnivory.1, mass_int, mass_top, mass_basal))
        
        print(test)
        
      } else {
        
        isolated = F
        
        local_web = local_web[to_keep,to_keep]
        community = cheddar::Community(nodes=data.frame(node=colnames(local_web)), trophic.links=PredationMatrixToLinks(local_web), properties=list(title='Community'))
        network <- graph.adjacency(local_web) ## igraph object     
        
        L <- sum(local_web) # nb of links
        S <- dim(local_web)[1] # nb of sp
        S.total <- length(sppSTRING) # nb of species not necessarily involved in the bird food web
        C <- L/(S**2) # connectance
        L.S <- L/S # nb of links per species

        # # degree distribution
        # DD <- rowSums(local_web) + colSums(local_web)
        # h = hist(DD, breaks = seq(1, S+1, 0.5))
        # nonZeroFreqs = which(h$counts > 0)
        # x = (h$breaks[nonZeroFreqs] + h$breaks[nonZeroFreqs + 1])/2
        # y = h$counts[nonZeroFreqs]
        # gauss = function(x,a,b,m,s){a*exp(-((x - m)^2/(2*s^2))) + b}
        # power = function(x,c,n){c*x^n}
        # 
        # GaussFit = nls(y~gauss(x,a,b,m,s), data = data.frame(x,y), start = list(a = 0.05, b = 0.01, m = mean(h$counts), 
        #                                                                     s = sqrt(var(DD))))
        # 
        # powerFit = nls(log(y)~power(x,c,n), data = data.frame(x,y), start = list(c = 1,n = 1))
        # powerFit.SSreg = sum((predict(powerFit) - mean(y))^2)
        # powerFit.SStot = sum((y - mean(y))^2)
        # powerFit.R = powerFit.SSreg/powerFit.SStot
        
        basal_sp = names(which(degree(network,mode='in') == 0 & degree(network,mode='out') > 0))
        top_sp = names(which(degree(network,mode='out') == 0 & degree(network,mode='in') > 0))
        int_sp = names(which(degree(network,mode='out') > 0 & degree(network,mode='in') > 0))
        
        mass_basal = mean(bodySize$BodyMass.Value[which(bodySize$Scientific %in% unlist(lapply(basal_sp, function(x) str_replace(whois(x), "_", " "))))], na.rm = T)
        mass_top = mean(bodySize$BodyMass.Value[which(bodySize$Scientific %in% unlist(lapply(top_sp, function(x) str_replace(whois(x), "_", " "))))], na.rm = T)
        mass_int = mean(bodySize$BodyMass.Value[which(bodySize$Scientific %in% unlist(lapply(int_sp, function(x) str_replace(whois(x), "_", " "))))], na.rm = T)
        
        ##### fraction of trophic levels
        Fbasal = length(which(degree(network,mode='in') == 0 & degree(network,mode='out') > 0))/S
        Ftop = length(which(degree(network,mode='out') == 0 & degree(network,mode='in') > 0))/S
        FInt = length(which(degree(network,mode='out') > 0 & degree(network,mode='in') > 0))/S
        Fother = 1 - (Fbasal + Ftop + FInt) # species that have no prey and no predators : should mow always be 0

        ##### number sp per trophic levels
        S_top <- NumberOfTop(local_web)
        S_intermediate <- NumberOfIntermediate(local_web)
        S_basal <- NumberOfBasal(local_web)

        Gen = Generality(local_web) # average nb of prey per sp: total nb of prey links / nb of species who have prey
        SDGen = SDGenerality(local_web)

        Vul = Vulnerability(local_web) # total nb of predator links / nb of species who are predators
        SDVul = SDVulnerability(local_web)

        out <- MeanFoodChainLength(network) # calculates both TL and mfcl (utils.r)
        mfcl = out$mfcl
        TL = data.frame(TL = out$ChainAveragedTrophicLevel)
        omnivory.1 = out$omnivory

        TL.mean = mean(TL$TL)

        cur_net <- graph_from_adjacency_matrix(local_web)
        modularity <- tryCatch({
          max(walktrap.community(cur_net)$modularity)
        }, warning = function(w) {NA}, error = function(e) {NA}, finally = {})

        # require(rnetcarto)
        # if(length(to_keep) <= 2){
        #   modularity.1 = 0
        #   modules.1 = 1
        # } else {
        #   mod.1 = netcarto(as.matrix(as_adjacency_matrix(cur_net)))
        #   modules.1 = length(unique(mod.1[[1]]$module))
        #   modularity.1 = mod.1[[2]]
        # }

        mod.2 = cluster_louvain(as.undirected(cur_net)) # igraph
        modularity.2 = modularity(mod.2)
        modules.2 = length(mod.2)

        # fraction of predatory links shared by predators
        overlap <- CalculatePredatorOverlap(local_web)
        # Omnivores are nodes that consume two or more species and have a non-integer trophic level (Polis 1991).
        # feed on more than one trophic level
        omnivory <- Omnivory(local_web)
        
        ### RETURN ALL METRICS
        test = cbind(pixel[c("PageName")],cbind(L,S,S.total,C,L.S,Gen,SDGen,Vul,SDVul,Fbasal,Ftop,FInt,Fother,S_basal,S_top,S_intermediate,overlap,omnivory,mfcl,
                                                modularity,TL.mean, isolated, modularity.2, modules.1, modules.2, omnivory.1, mass_int, mass_top, mass_basal))
        
        print(test)
      }
      
    }
  }
  close(pb)
  stopCluster(cl) 
  
  return(cur.dat)}








## Function that returns list of species and their trophic level inside the community


get_species_traits = function(master_reduced, iterations, m_code, bodySize){
  
  ## PARALLEL LOOP
  library(foreach)
  library(doParallel)
  library(doSNOW) # for progression bar
  
  
  cores <- detectCores() - 2
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # this uses the nb of cores set in registerdopareallel
  cur.dat = foreach(cell = 1:iterations, .combine=rbind, .options.snow = opts) %dopar% { 
    # for(cell in 1:iterations){
    ## INITIALISATION
    
    species_code <- TL_ordered <- trophic_class <- body_mass <- NULL
    

    library(stringr)
    library(igraph)
    library(cheddar)
    source("E:/TheseSwansea/Galiana2021_network-area-europe-master/utils.R")
    source("data/Galiana2021_network-area-europe-master/whois-function.R")
    
    
    
    pixel = master_reduced[cell,] ## local community (sp that occur within this 10x10km grid cell)
    sppSTRING <- names(pixel)[which(pixel==1 & str_detect(names(pixel), pattern = "[:digit:]"), arr.ind=FALSE)]
    
    
    
    ## CREATION OF LOCAL WEB
    if (length(sppSTRING)>1){
      
      local_web = m_code[sppSTRING,sppSTRING]
      sum(diag(local_web)) == 0
      
      # creating a cheddar community to look for isolated nodes
      community = cheddar::Community(nodes=data.frame(node=colnames(local_web)), trophic.links=PredationMatrixToLinks(local_web), properties=list(title='Community'));
      
      # new local web that only keeps connected nodes
      to_keep = setdiff(colnames(local_web),IsolatedNodes(community)) # species that are not 'isolated' i.e. have at least 1 prey or 1 pred
      
      #### networks with no links: 
      if (length(to_keep)==0){ 
        
        S <- 0 # nb of sp
        S.total <- length(sppSTRING) # nb of species not necessarily involved in the bird food web
        
        species_code <- TL_ordered <- trophic_class <- body_mass <- NA
        
        isolated = T
        
        test = cbind(pixel[c("PageName")],cbind(species_code, TL_ordered, trophic_class, body_mass))
        
        print(test)
        
      } else {
        
        isolated = F
        
        local_web = local_web[to_keep,to_keep]
        community = cheddar::Community(nodes=data.frame(node=colnames(local_web)), trophic.links=PredationMatrixToLinks(local_web), properties=list(title='Community'))
        network <- graph.adjacency(local_web) ## igraph object     
        
        species_code = colnames(local_web)
        
        basal_sp = names(which(degree(network,mode='in') == 0 & degree(network,mode='out') > 0))
        top_sp = names(which(degree(network,mode='out') == 0 & degree(network,mode='in') > 0))
        int_sp = names(which(degree(network,mode='out') > 0 & degree(network,mode='in') > 0))

        trophic_class = ifelse(species_code %in% basal_sp, 'basal_sp', ifelse(species_code %in% top_sp, 'top_sp','int_sp'))
        
        # make sure to re order with match()        
        body_mass = bodySize$BodyMass.Value[match(unlist(lapply(species_code, function(x) str_replace(whois(x), "_", " "))), bodySize$Scientific)]
        
        out <- MeanFoodChainLength(network) # calculates both TL and mfcl (utils.r)
        TL = data.frame(TL = out$ChainAveragedTrophicLevel)
        TL_ordered = TL[match(species_code, rownames(TL)),] # re order so that it matches order of species_code

        test = cbind(pixel[c("PageName")],cbind(species_code, TL_ordered, trophic_class, body_mass))
        print(test)
      }
      
    }
  }
  close(pb)
  stopCluster(cl) 
  
  return(cur.dat)}
