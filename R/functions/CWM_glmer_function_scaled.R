# function to perform the CWM analyses employing a GLMM model
# comm= community data
# envir= environmental data
# formula= the lme4 formula with random terms
# PGLS= if the PGLS should be performed
# AsFactors= which environmental traits are factors
# runs= number of permutations
# random= a string with the name of the random term

# formula = as.formula(paste(colnames(traits.all[i]), "~", 
#                            colnames(envir[x1]), "+(1|", 
#                            colnames(envir[r1]), ")", 
#                            collapse=""))

cwm.sig.glmm_scaled <- function(comm, traits, envir, tree, formula, 
                        AsFactors, PGLS, 
                        runs, random.term){
  
  comm <- SYNCSA::organize.syncsa(comm = comm, phylodist = cophenetic(tree))$community
  traits <- SYNCSA::organize.syncsa(comm = comm, traits = traits)$traits
  comm <- sqrt(comm)
  RES <- list(call = match.call())
  envir <- as.data.frame(envir)
  if (!is.null(AsFactors)) {
    for (i in AsFactors) {
      envir[, i] <- as.factor(envir[, i])
    }
  }
  envir.class <- SYNCSA::var.type(envir)
  RES$envir.class <- envir.class
  trait.name <- as.character(stats::as.formula(formula)[[2]]) # trait
  trait <- traits %>% select(all_of(trait.name)) %>% 
    mutate(Species = rownames(.))
  
  RES$formula <- formula
  RES$PGLS <- PGLS
  if(PGLS == TRUE){
    cor.phy <- ape::corBrownian(phy = tree, form = ~Species)
    model.pgls <- as.formula(paste(trait.name, "~1",
                                   collapse = ""))
    fit.pgls <- nlme::gls(model = model.pgls, data = as.data.frame(trait), 
                          correlation = cor.phy)
    RES$PGLS.model <- fit.pgls
    res.fit.pgls <- as.data.frame(residuals(fit.pgls, type =
                                              "normalized"))
    colnames(res.fit.pgls) <- colnames(trait)[1]
    RES$trait.residuals <- res.fit.pgls
    trait.temp <- scale(res.fit.pgls) # scale traits before CWM
  } else {
    trait.temp <- scale(traits) # scale traits before CWM
  }
  
  MT <- SYNCSA::matrix.t(comm, trait.temp, scale =FALSE)$matrix.T
  data.obs <- as.data.frame(cbind(MT, envir))
  mod.obs <- lmerTest::lmer(formula, data = data.obs)
  t.obs <- matrix(NA, 1, 2)
  colnames(t.obs) <- c("Intercept", "Slope")
  t.obs[,1] <- base::summary(mod.obs)$coefficients[7]
  t.obs[,2] <- base::summary(mod.obs)$coefficients[8]
  RES$model <- mod.obs
  RES$statistic.obs <- t.obs
  t.null.site <- matrix(NA, runs, 2)
  t.null.trait <- matrix(NA, runs, 2)
  colnames(t.null.site) <- colnames(t.null.trait) <- colnames(t.obs)
  posit <- which(colnames(envir) == random.term)
  mt <- as.data.frame(cbind(MT, envir %>% select(random.term)))
  
  for (i in 1:runs){
    MT.null.site <- SYNCSA::permut.row.matrix(mt)$permut.matrix
    data.null.site <- as.data.frame(cbind(MT.null.site, envir %>% select(-random.term)))
    
    mod.null.site <- lmerTest::lmer(formula, data =
                                      data.null.site)
    t.null.site[i,1] <- base::summary(mod.null.site)$coefficients[7]    
    t.null.site[i,2] <- base::summary(mod.null.site)$coefficients[8]
    
    MT.null.trait <- SYNCSA::matrix.t(comm,
                                      SYNCSA::permut.row.matrix(trait.temp)$permut.matrix, scale =
                                        FALSE)$matrix.T
    data.null.trait <- as.data.frame(cbind(MT.null.trait,
                                           envir))
    mod.null.trait <- lmerTest::lmer(formula, data =
                                       data.null.trait)
    t.null.trait[i,1] <- base::summary(mod.null.trait)$coefficients[7]
    t.null.trait[i,2] <- base::summary(mod.null.trait)$coefficients[8]
  }
  prob.site.shuffle.intercept <- (sum(ifelse(t.null.site[,1]>=t.obs[,1], 1,
                                             0))+1)/(runs+1)
  prob.site.shuffle.slope <- (sum(ifelse(t.null.site[,2]>=t.obs[,2], 1,
                                         0))+1)/(runs+1)
  
  RES$p.site.shuffle <- as.vector(rbind(prob.site.shuffle.intercept,
                                        prob.site.shuffle.slope))
  
  prob.trait.shuffle.intercept <-(sum(ifelse(t.null.trait[,1]>=t.obs[,1],
                                             1, 0))+1)/(runs+1)
  prob.trait.shuffle.slope <-(sum(ifelse(t.null.trait[,2]>=t.obs[,2],
                                         1, 0))+1)/(runs+1)
  RES$p.trait.shuffle <- as.vector(rbind(prob.trait.shuffle.intercept, 
                                         prob.trait.shuffle.slope))
  
  return(RES)
}
