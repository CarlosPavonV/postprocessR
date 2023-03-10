sampleBPPtrees <- function(infile,ntrees){
  #Takes the FigTree-readable output of BPP analysis A00 and generates a sample of trees (in a multiPhylo object) based on the 95% HPD of divergence times
  #infile = character string, name of file containing the FigTree-readable output of BPP
  #ntrees = integer, number of trees to return
  require(MCMCtreeR)
  require(BioGeoBEARS)
  require(phytools)
  arb <- readMCMCtree(infile)
  arb.ape <- arb$apePhy
  arb.list <- as.list(1:ntrees)
  for (j in 1:ntrees){
    daughter.table <- matrix(ncol=6,nrow=arb.ape$Nnode)
    colnames(daughter.table) <- c("Node","Desc1","Desc2","Length1","Length2","Parent")
    for (i in 1:arb.ape$Nnode){
      daughter.table[i,1] <- Ntip(arb.ape)+i
      daughter.table[i,2:3] <- get_daughters(Ntip(arb.ape)+i,arb.ape)
      daughter.table[i,6] <- get_parent(Ntip(arb.ape)+i,arb.ape)
    }
    ages <- runif(1,arb$nodeAges[1,2],arb$nodeAges[1,3])#Get root age
    names(ages) <- Ntip(arb.ape)+1
    for (i in 2:arb$apePhy$Nnode){
      max.age <- arb$nodeAges[i,3]
      if(max.age>ages[which(names(ages)==as.character(daughter.table[i,6]))]){
        new.max <- (ages[which(names(ages)==as.character(daughter.table[i,6]))])-((ages[which(names(ages)==as.character(daughter.table[i,6]))])/1000)
      } else{
        new.max <- max.age
      }
      new.age <- runif(1,arb$nodeAges[i,2],new.max)
      ages <- c(ages,new.age)
      names(ages)[i] <- Ntip(arb.ape)+i
    }
    
    for (i in which(daughter.table[,2]%in%daughter.table[,1])){
      daughter.table[i,4] <- ages[which(names(ages)==as.character(daughter.table[i,1]))]-ages[which(names(ages)==as.character(daughter.table[i,2]))]
    }
    for (i in which(daughter.table[,3]%in%daughter.table[,1])){
      daughter.table[i,5] <- ages[which(names(ages)==as.character(daughter.table[i,1]))]-ages[which(names(ages)==as.character(daughter.table[i,3]))]
    }
    for (i in which(daughter.table[,2]%in%(1:length(arb.ape$tip.label)))){
      daughter.table[i,4] <- ages[which(names(ages)==as.character(daughter.table[i,1]))]
    }
    for (i in which(daughter.table[,3]%in%(1:length(arb.ape$tip.label)))){
      daughter.table[i,5] <- ages[which(names(ages)==as.character(daughter.table[i,1]))]
    }
    arbol.tmp <- arb.ape
    for (i in 1:length(arbol.tmp$edge.length)){
      col.query <- which(daughter.table[which(daughter.table[,1]==arbol.tmp$edge[i,1]),]==arbol.tmp$edge[i,2])
      arbol.tmp$edge.length[i] <- daughter.table[which(daughter.table[,1]==arbol.tmp$edge[i,1]),(col.query+2)]
    }
    if(!is.ultrametric(arbol.tmp)){
      arbol.tmp <- force.ultrametric(arbol.tmp,method="nnls")
    }
    arb.list[[j]] <- arbol.tmp
  }
  class(arb.list) <- "multiPhylo"
  return(arb.list)
}
