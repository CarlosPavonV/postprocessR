biogeobears2simmap <- function(phy,ana,clado){
  #Creates a simmap object from the output of stochastic character mapping with BioGeoBears
  #Currently handles a single stochastic map, but can be used for more than one map using a for loop or lapply
  #phy = the tree used for the BioGeoBears analysis in phylo format
  #ana = a data frame corresponding to an ana_events_table from BioGeoBears
  #clado = a data frame corresponding to a clado_events_table from BioGeoBears
  require(phytools)
  require(phangorn)
  '%notin%' <- Negate('%in%')
  phy.sim <- phy
  not.symp <- which(clado$clado_event_type!="sympatry (y)")
  not.symp <- not.symp[not.symp>length(phy$tip.label)]
  for (i in not.symp){
    ance <- paste((strsplit(strsplit(clado$clado_event_txt[i],"->")[[1]][1],"")[[1]]),collapse="+")
    daug1 <- gsub("(?<=.)(?=.)","+",(strsplit(strsplit(clado$clado_event_txt[i],"->")[[1]][2],",")[[1]][1]),perl=TRUE)
    daug2 <- gsub("(?<=.)(?=.)","+",(strsplit(strsplit(clado$clado_event_txt[i],"->")[[1]][2],",")[[1]][2]),perl=TRUE)
    phy.sim <- paintSubTree(phy.sim,
                            clado$daughter_nds[i][[1]][1],
                            anc.state=ance,
                            state=daug1,
                            stem=T)
    phy.sim <- paintSubTree(phy.sim,
                            clado$daughter_nds[i][[1]][2],
                            anc.state=ance,
                            state=daug2,
                            stem=T)
  }
  rownames(ana) <- 1:nrow(ana)
  ana.tip <- which(ana$node.type=="tip")
  for (i in ana.tip){
    daug <- gsub("(?<=.)(?=.)","+",ana$new_rangetxt[i],perl=TRUE)
    phy.sim <-paintSubTree(phy.sim,ana$node[i],state=daug,
                           stem=((phy.sim$edge.length[which(phy.sim$edge[,2]==ana$node[i])]-ana$event_time[which(ana$node==ana$node[i])])
                                 /phy.sim$edge.length[which(phy.sim$edge[,2]==ana$node[i])]))
  }
  ana.int <- which(rownames(ana) %notin% ana.tip)
  for (i in ana.int){
    node.term <- ana$node[i]
    daug <- gsub("(?<=.)(?=.)","+",ana$new_rangetxt[i],perl=TRUE)
    daug.t <- (phy.sim$edge.length[which(phy.sim$edge[,2]==ana$node[i])])-(ana$event_time[which(ana$node==ana$node[i])])
    rowname.list <- strsplit(rownames(phy.sim$mapped.edge),",")
    edge.last <- c()
    for (j in 1:length(rowname.list)){
      if (as.character(node.term)==rowname.list[[j]][2]){
        edge.last <- j
      }
    }
    edge.last <- as.integer(edge.last)
    if(daug %in% colnames(phy.sim$mapped.edge)){
      phy.sim$mapped.edge[edge.last,which(phy.sim$mapped.edge[edge.last,]>0)] <- (phy.sim$mapped.edge[edge.last,which(phy.sim$mapped.edge[edge.last,]>0)])-daug.t
      phy.sim$mapped.edge[edge.last,which(colnames(phy.sim$mapped.edge)==daug)] <- daug.t
    } else{
      phy.sim$mapped.edge <- cbind(phy.sim$mapped.edge,0)
      colnames(phy.sim$mapped.edge) <- c(colnames(phy.sim$mapped.edge)[-ncol(phy.sim$mapped.edge)],daug)
      phy.sim$mapped.edge[edge.last,which(phy.sim$mapped.edge[edge.last,]>0)] <- (phy.sim$mapped.edge[edge.last,which(phy.sim$mapped.edge[edge.last,]>0)])-daug.t
      phy.sim$mapped.edge[edge.last,which(colnames(phy.sim$mapped.edge)==daug)] <- daug.t
    }
    phy.sim$maps[[edge.last]] <- phy.sim$mapped.edge[edge.last,which(phy.sim$mapped.edge[edge.last,]>0)]
  }
  nodes.biogeo <- sapply(1:nrow(clado),function (i) {paste((strsplit(strsplit(clado$clado_event_txt[i],"->")[[1]][1],"")[[1]]),collapse="+")})
  node.states <- data.frame(node=names(getStates(phy.sim,"nodes")),in.phylo=getStates(phy.sim,"nodes"),in.biogeo=nodes.biogeo[(length(phy.sim$tip.label)+1):length(nodes.biogeo)],stringsAsFactors=F)
  rownames(node.states) <- 1:nrow(node.states)
  node.states$node <- as.integer(node.states$node)
  bad.nodes <- c()
  for (i in 1:nrow(node.states)){
    if (node.states$in.phylo[i]!=node.states$in.biogeo[i]){
      bad.nodes <- c(bad.nodes,node.states$node[i])
    }
  }
  for (i in bad.nodes){
    desc1 <- clado$daughter_nds[i][[1]][1]
    desc2 <- clado$daughter_nds[i][[1]][2]
    edge1 <- which(phy.sim$edge[,1]==i&phy.sim$edge[,2]==desc1)
    edge2 <- which(phy.sim$edge[,1]==i&phy.sim$edge[,2]==desc2)
    if(names(phy.sim$maps[[edge1]][which(phy.sim$maps[[edge1]]>0)])==node.states$in.biogeo[which(node.states$node==i)]){
      phy.sim$maps[[edge1]] <- phy.sim$maps[[edge1]][which(phy.sim$maps[[edge1]]>0)]
    } else{
      if (names(phy.sim$maps[[edge2]][which(phy.sim$maps[[edge2]]>0)])==node.states$in.biogeo[which(node.states$node==i)]){
        phy.sim$maps[[edge2]] <- phy.sim$maps[[edge2]][which(phy.sim$maps[[edge2]]>0)]
      } else{
        if(node.states$in.biogeo[which(node.states$node==i)] %notin% names(phy.sim$maps[[edge1]])){
          names.edge1 <- c(node.states$in.biogeo[which(node.states$node==i)],names(phy.sim$maps[[edge1]]))
          phy.sim$maps[[edge1]] <- c(0,phy.sim$maps[[edge1]])
          names(phy.sim$maps[[edge1]]) <- names.edge1
        } else{
          if(node.states$in.biogeo[which(node.states$node==i)] %notin% names(phy.sim$maps[[edge2]])){
            names.edge2 <- c(node.states$in.biogeo[which(node.states$node==i)],names(phy.sim$maps[[edge2]]))
            phy.sim$maps[[edge2]] <- c(0,phy.sim$maps[[edge2]])
            names(phy.sim$maps[[edge2]]) <- names.edge2
          }
        }
      }
    }
  }
  node.states2 <- data.frame(node=names(getStates(phy.sim,"nodes")),in.phylo=getStates(phy.sim,"nodes"),in.biogeo=nodes.biogeo[(length(phy.sim$tip.label)+1):length(nodes.biogeo)],stringsAsFactors=F)
  rownames(node.states2) <- 1:nrow(node.states2)
  node.states2$node <- as.integer(node.states2$node)
  bad.nodes2 <- c()
  for (i in 1:nrow(node.states2)){
    if (node.states2$in.phylo[i]!=node.states2$in.biogeo[i]){
      bad.nodes2 <- c(bad.nodes2,node.states2$node[i])
    }
  }
  if(length(bad.nodes2)>0){
    for (i in bad.nodes2){
      desc1 <- clado$daughter_nds[i][[1]][1]
      desc2 <- clado$daughter_nds[i][[1]][2]
      edge1 <- which(phy.sim$edge[,1]==i&phy.sim$edge[,2]==desc1)
      edge2 <- which(phy.sim$edge[,1]==i&phy.sim$edge[,2]==desc2)
      if ((node.states$in.biogeo[which(node.states$node==i)] %notin% names(phy.sim$maps[[edge1]]))
          && (node.states$in.biogeo[which(node.states$node==i)]!=getStates(phy.sim,"tips")[desc1])){
        names.edge1 <- c(node.states$in.biogeo[which(node.states$node==i)],names(phy.sim$maps[[edge1]]))
        phy.sim$maps[[edge1]] <- c(0,phy.sim$maps[[edge1]])
        names(phy.sim$maps[[edge1]]) <- names.edge1
      } else{
        names.edge2 <- c(node.states$in.biogeo[which(node.states$node==i)],names(phy.sim$maps[[edge2]]))
        phy.sim$maps[[edge2]] <- c(0,phy.sim$maps[[edge2]])
        names(phy.sim$maps[[edge2]]) <- names.edge2
      }
    }
  }
  node.states3 <- data.frame(node=names(getStates(phy.sim,"nodes")),in.phylo=getStates(phy.sim,"nodes"),in.biogeo=nodes.biogeo[(length(phy.sim$tip.label)+1):length(nodes.biogeo)],stringsAsFactors=F)
  rownames(node.states3) <- 1:nrow(node.states3)
  node.states3$node <- as.integer(node.states3$node)
  bad.nodes3 <- c()
  for (i in 1:nrow(node.states3)){
    if (node.states3$in.phylo[i]!=node.states3$in.biogeo[i]){
      bad.nodes3 <- c(bad.nodes3,node.states3$node[i])
    }
  }
  if (length(bad.nodes3)>0){
    print("Sorry, operation failed, there must be a bug in my code :(")
  } else{
    return(phy.sim)
  }
}
