# function that perform a test for brownian phylogenetic signal for traits 
# (binary and/or continuous) to perform the better correct the CWM analyses 
# employing a GLMM model.

# comm = community data
# traits = matrix with the functional traits
# tree = the phylogenetic tree
# envir = environmental data
# binary = logical argument, describe if there are binary traits in the matrix,
  # if binary = TRUE each set of trait is treated distinctly; if binary = FALSE
  # all traits is consider continuous.
# AsFactors = which environmental traits are factors
# runs = number of permutations for null models (CWM and phylogentic signal)
# predictor = a string with the name of predictor variable
# random.term = a string with the name of the random term

physig.cwm.function <- function(comm, traits, tree, envir,
                               binary, runs, AsFactors, 
                               predictor, random.term){
  
  x1 <- which(colnames(envir) == predictor)
  r1 <- which(colnames(envir) == random.term)
  
  if (!is.null(AsFactors)) {
    envir.num<- envir
    for (i in AsFactors) {
      envir[, i] <- as.factor(envir[, i])
      envir.num[, i]<- as.numeric(as.factor(envir[, i]))
    }
  }
  
  dist.phylo <- cophenetic(tree)
  org <- organize.pcps(comm = comm, phylodist = dist.phylo, envir = envir.num)
  comm <- org$community
  envir.num <- org$environmental
  dist.phylo <- org$phylodist
  
  ## test for phylogenetic signal in the species distribution
  pred.physig <- as.formula(paste("p.dist", "~", colnames(envir.num[x1]), "+", 
                                 colnames(envir.num[r1]),
                                 collapse = ""))
  phylo.sig.comm <- matrix.p.sig(comm = comm, phylodist = dist.phylo, FUN = FUN.ADONIS2.global, 
                                envir = envir.num, formula = pred.physig, 
                                method.p = "bray", runs = runs)
  
  # test for phylogenetic signal in traits
  if(binary == TRUE){
    bin <- vector()
    for(i in 1:ncol(traits)){
      bin[i] <- is.integer(traits[,i])
    }
    traits.cont <- traits[, which(bin == FALSE)]
    traits.bin <- traits[, which(bin == TRUE)]
    
    # phylogenetic signal for continuous traits (K statistic)
    phy.sig.k <- matrix(NA, nrow = ncol(traits.cont), ncol = 2) # matrix to store the phylogenetic signal for continuous traits
    rownames(phy.sig.k) <- colnames(traits.cont)
    colnames(phy.sig.k) <- c("K", "p")
    
    for (i in 1:ncol(traits.cont)) {
      trait.tmp <- traits.cont[,i]
      names(trait.tmp) <- rownames(traits.cont)
      
      k.trait.tmp <- phylosig(tree = tree, x = trait.tmp, method = "K", test = T)
      phy.sig.k[i,1] <- k.trait.tmp$K
      phy.sig.k[i,2] <- k.trait.tmp$P
    }
    
    # fit the phylogenetic signal for binary traits (Mantel with evolutionary model)
    phy.sig.em <- matrix(NA, nrow = ncol(traits.bin), ncol = 3) # matrix to store the phylogenetic signal for binary traits
    rownames(phy.sig.em) <- colnames(traits.bin)
    colnames(phy.sig.em) <- c("r.Mantel", "p.NULL", "p.BM")
    
    for (i in 1:ncol(traits.bin)) {
      trait.tmp <- traits.bin[,i, drop= FALSE]
      em.trait.tmp <- EM.mantel(tree = tree, traits = trait.tmp, runs = runs)
      phy.sig.em[i,1] <- em.trait.tmp$r.Mantel
      phy.sig.em[i,2] <- em.trait.tmp$p.NULL
      phy.sig.em[i,3] <- em.trait.tmp$p.BM
    }
    
    phylo.sig.trait <- list(k.statistic = phy.sig.k, EM.Mantel = phy.sig.em) # results of phylo tests
    
    # transforming the binary traits in continuous traits, while removing the phylogenetic component
    x <- PVR::PVRdecomp(phy = tree) # Creating the phylogenetic vectors
    
    # creating a matrix to store the vectors
    traits.pvr <- matrix(NA, nrow = nrow(traits), ncol = ncol(traits.bin))
    colnames(traits.pvr) <- colnames(traits[, which(bin == TRUE)])
    
    for (i in 1:ncol(traits.bin)) {
      pvr <- PVR::PVR(x = x, trait = traits.bin[, i], method = "moran")
      vect <- pvr@Selection$Vectors
      glm.pvr <- glm(traits.bin[, i] ~ vect, family = binomial)
      traits.pvr[, i] <- as.matrix(residuals(glm.pvr, type = "response"))
    } # phylogenetic free traits
    
    traits.all <- cbind(traits.cont, traits.pvr) # traits matrix input to CWM
    pgls.option <- c(phy.sig.k[, 1] >= 1, rep(FALSE, ncol(traits.all) - ncol(traits.cont)))
    
    ## run the cwm.sig.glmm for all traits
    cwm.output <- matrix(NA, ncol = ncol(traits.all), nrow = 6)
    rownames(cwm.output) <- c("Intercept.obs", "Slope.obs", "p.site.shuffle.i",
                             "p.site.shuffle.s", "p.trait.shuffle.i", "p.trait.shuffle.s")
    colnames(cwm.output) <- colnames(traits.all)
    
    if (phylo.sig.comm$p.taxa.shuffle[2] < 0.05){ # if assemblages are phylogenetically structured
      for (i in 1:ncol(traits.all)) {
        cwm.model <- cwm.sig.glmm(comm = comm, traits = traits.all, envir = envir, tree = tree,
                                 formula = as.formula(paste(colnames(traits.all[i]), "~", 
                                                            colnames(envir[x1]), "+(1|", 
                                                            colnames(envir[r1]), ")", 
                                                            collapse="")), 
                                 random.term = random.term, PGLS = pgls.option[i], 
                                 runs = runs, AsFactors = AsFactors)
        
        cwm.output[, i] <- round(c(cwm.model$statistic.obs, cwm.model$p.site.shuffle, 
                                  cwm.model$p.trait.shuffle), digits = 4)
      }
      
    } else { # if the assemblages aren't phylogenetically structures
      for (i in 1:ncol(traits.all)) {
        cwm.model <- cwm.sig.glmm(comm = comm, traits = traits.all, envir = envir, tree = tree,
                                 formula = as.formula(paste(colnames(traits.all[i]), "~", 
                                                            colnames(envir[x1]), "+(1|", 
                                                            colnames(envir[r1]), ")", 
                                                            collapse="")), 
                                 random.term = random.term, PGLS = FALSE, 
                                 runs = runs, AsFactors = AsFactors)
        
        cwm.output[, i] <- round(c(cwm.model$statistic.obs, cwm.model$p.site.shuffle, 
                                  cwm.model$p.trait.shuffle), digits = 4)
      }
    }
    
  }
  else { # if there aren't binary traits in the traits matrix
    phy.sig.k <- matrix(NA, nrow = ncol(traits), ncol = 2)
    rownames(phy.sig.k) <- colnames(traits)
    colnames(phy.sig.k) <- c("K", "p")
    
    for (i in 1:ncol(traits)) {
      trait.tmp <- traits[, i]
      names(trait.tmp) <- rownames(traits)
      
      k.trait.tmp<- phylosig(tree = tree, x = trait.tmp, method = "K", test = T)
      phy.sig.k[i, 1] <- k.trait.tmp$K
      phy.sig.k[i, 2] <- k.trait.tmp$P
    }
    
    phylo.sig <- list(k.statistic = phy.sig.k)
    pgls.option <- phy.sig.k[, 1] >= 1
    
    ## run the cwm.sig.glmm for all traits
    cwm.output <- matrix(NA, ncol = ncol(traits.all), nrow = 6)
    rownames(cwm.output) <- c("Intercept.obs", "Slope.obs", "p.site.shuffle.i",
                             "p.site.shuffle.s", "p.trait.shuffle.i", "p.trait.shuffle.s")
    colnames(cwm.output) <- colnames(traits.all)
    
    if (phylo.sig.comm$p.taxa.shuffle[2] < 0.05){
      for (i in 1:ncol(traits.all)) {
        cwm.model <- cwm.sig.glmm(comm = comm, traits = traits.all, envir = envir, tree = tree,
                                 formula = as.formula(paste(colnames(traits.all[i]), "~", 
                                                            colnames(envir[x1]), "+(1|", 
                                                            colnames(envir[r1]), ")", 
                                                            collapse="")), 
                                 random.term = random.term, PGLS = pgls.option[i], 
                                 runs = runs, AsFactors = AsFactors)
        
        cwm.output[, i] <- round(c(cwm.model$statistic.obs, cwm.model$p.site.shuffle, 
                                  cwm.model$p.trait.shuffle), digits = 4)
      }
      
    } else {
      for (i in 1:ncol(traits.all)) {
        cwm.model <- cwm.sig.glmm(comm = comm, traits = traits.all, envir = envir, tree = tree,
                                 formula = as.formula(paste(colnames(traits.all[i]), "~", 
                                                            colnames(envir[x1]), "+(1|", 
                                                            colnames(envir[r1]), ")", 
                                                            collapse="")), 
                                 random.term = random.term, PGLS = FALSE, 
                                 runs = runs, AsFactors = AsFactors)
        
        cwm.output[, i] <- round(c(cwm.model$statistic.obs, cwm.model$p.site.shuffle, 
                                  cwm.model$p.trait.shuffle), digits = 4)
      }
    }
  }
  
  return(list(CWM.output = cwm.output, 
              Trait.phylosig = phylo.sig.trait,
              Comm.phylosig = phylo.sig.comm))
  
}
