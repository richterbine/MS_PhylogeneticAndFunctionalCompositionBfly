EM.mantel<-function(tree, traits, runs = 999, euclidean= TRUE, sqrtPhylo=FALSE, checkdata = TRUE, ...){
  phylo.dist<-cophenetic(tree)
  if(sqrtPhylo){
    phylo.dist<-sqrt(phylo.dist)
  }
  if(checkdata){
    if(is.null(tree$tip.label)){
      stop("\n Error in tip labels of tree\n")
    }
    if(is.null(rownames(traits))){
      stop("\n Error in row names of traits\n")
    }
    match.names <- match(rownames(traits),rownames(phylo.dist))
    if(sum(is.na(match.names)) > 0){
      stop("\n There are species from traits data that are not on phylogenetic tree\n")
    }
    phylo.dist <- phylo.dist[match.names, match.names]
  }
  if(length(tree$tip.label) > dim(traits)[1]){
    warning("Tree have more species that species in traits data")
  }
  if(dim(phylo.dist)[1] != dim(traits)[1] & checkdata == FALSE){
    stop("\n Different number of species in tree and in traits data, use checkdata = TRUE\n")
  }
  gow.dist<-gowdis(traits, ...)
  if(euclidean){
    gow.sim<-1-gow.dist
    gow.dist<-sqrt(1-gow.sim)
  }
  traits.attr<-attr(gow.dist, "Types", exact = TRUE)
  res.mantel<-mantel(phylo.dist,gow.dist,permutations=runs)
  res.BM<-matrix(NA,runs,1)
  for(k in 1:runs){
    traits_sim<-matrix(NA,length(tree$tip.label),dim(traits)[2])
    rownames(traits_sim)<-tree$tip.label
    for(i in 1:dim(traits)[2]){
      traits_sim[,i]<-rTraitCont(tree,model="BM")
    }
    traits_sim<-decostand(traits_sim,method="standardize",MARGIN=2)
    traits_sim<-as.data.frame(traits_sim)	
    for(i in 1:dim(traits)[2]){
      if(traits.attr[i] == "B" | traits.attr[i] == "A"){
        probs<-sum(traits[,i])/dim(traits)[1]
        threshold<-quantile(traits_sim[,i],probs=1-probs)
        traits_sim[,i]<-ifelse(traits_sim[,i]>=threshold,1,0)
      }
      if(traits.attr[i] == "N" | traits.attr[i] == "O"){
        n.levels<-length(levels(traits[,i]))
        traits.levels<-levels(traits[,i])
        probs<-cumsum(table(traits[,i]))/sum(table(traits[,i]))
        probs<-probs[1:(n.levels-1)]
        threshold<-quantile(traits_sim[,i],probs=probs)
        threshold<-c(min(traits_sim[,i]),threshold,max(traits_sim[,i]))
        temp<-matrix(NA,length(traits_sim[,i]),1)
        for(j in 1:n.levels){
          if(j < n.levels){
            temp[1:length(traits_sim[,i]),1]<-ifelse(traits_sim[,i]>=threshold[j] & traits_sim[,i]<threshold[j+1], traits.levels[j],temp)
          }
          if(j == n.levels){
            temp[1:length(traits_sim[,i]),1]<-ifelse(traits_sim[,i]>=threshold[j] & traits_sim[,i]<=threshold[j+1], traits.levels[j],temp)
          }
        }
        traits_sim[,i]<-as.factor(temp)
        if(traits.attr[i] == "O"){
          traits_sim[,i]<-ordered(temp,levels=levels(traits[,i]))
        }
      }
    }
    if(checkdata == TRUE){
      match.names <- match(rownames(traits),rownames(traits_sim))
      traits_sim<-traits_sim[match.names,,drop=FALSE]
    }
    gow.dist.BM<-gowdis(traits_sim, ...)
    if(euclidean){
      gow.sim.BM<-1-gow.dist.BM
      gow.dist.BM<-sqrt(1-gow.sim.BM)
    }
    res.mantel.BM<-mantel(phylo.dist,gow.dist.BM,permutations=0)
    res.BM[k,1]<-res.mantel.BM$statistic
  }
  p.BM<-(sum(ifelse(res.BM[,1]>=res.mantel$statistic,1,0))+1)/(runs+1)
  p.NULL<-res.mantel$signif
  r.Mantel<-res.mantel$statistic
  RES<-list(perm.NULL=res.mantel$perm,perm.BM=res.BM[,1],r.Mantel=r.Mantel,p.NULL=p.NULL,p.BM=p.BM)
  return(RES)
}
