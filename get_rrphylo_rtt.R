get_rrphylo_rtt <- function(rr,nsamples,statistic){
  #Estimates a summary statistic for the absolute rates through time from the RRphylo output
  #Right now it can only handle univariate data and extant-only trees
  #rr = list, RRphylo output
  #nsamples = numeric, the number of points at time where the mean rate is estimated
  #statistic = function, the R function used to calculate the summary statistic
  require(dispRity)
  tree <- rr$tree
  ages <- tree.age(tree, order = "past", fossil = TRUE, digits = 10)
  ages$elements <- 1:nrow(ages)
  ages$ages[1:Ntip(tree)] <- 0
  branch.table <- data.frame(parent=rep(0,length(tree$edge.length)),
                             parent.age=rep(0,length(tree$edge.length)),
                             daughter=rep(0,length(tree$edge.length)),
                             daughter.age=rep(0,length(tree$edge.length)))
  branch.table$parent <- tree$edge[,1]
  branch.table$daughter <- tree$edge[,2]
  branch.table$parent.age <- ages$ages[match(branch.table$parent,ages$elements)]
  branch.table$daughter.age <- ages$ages[match(branch.table$daughter,ages$elements)]
  time.increase <- max(ages$ages)/(nsamples-1)
  times <- seq(max(ages$ages),0,by=-time.increase)
  rates <- rr$rates
  rates[,1] <- abs(rates[,1])
  row.names(rates)[(Nnode(tree)+1):nrow(rates)] <- 1:Ntip(tree)
  rate.list <- list()
  rate.list[[1]] <- rates[1,1]
  for(i in 2:length(times)){
    branches.tmp <- branch.table$daughter[which(branch.table$daughter.age<=times[i] & branch.table$parent.age>times[i])]
    rate.list[[i]] <- rates[which(as.numeric(row.names(rates)) %in% branches.tmp),1]
  }
  summarized.rates <- data.frame(time_before_present=times,
                                 statistic_of_absolute_rate=unlist(lapply(rate.list,statistic)))
  return(summarized.rates)
}