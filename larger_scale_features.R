library(ape)
library(phangorn)
library(seqinr)

#' @param tree a phylogenetic tree
#' @param maxdate - the date of the most recent tip in the tree, in numerical format (or a format that can be added to a number) 
#' @return a 2-column matrix, time and number of splitting events before that time. times are splitting times of the tree
divtt <- function(tree, maxdate=NULL) {
  x = ltt.plot.coords(tree,backward = F)
  # i want to take all the negative increments out and just have the births
  # remove rows of x where N decreases. then the num divs is just 1: (num events) 
  toRemove = 1+which(diff(x[,2])<0)
  x=x[-toRemove,]
  x[,2]=1:nrow(x)
  if (!is.null(maxdate)) {x[,1]= x[,1]+ maxdate - max(node.depth.edgelength(tree))} 
  return(data.frame(time=x[,1], Nsplits=x[,2]))
}

#' @param tree a phylogenetic tree
#' @param numbins number of bins you want to group the times into
#' @return list with group (the bin number of the tips) and times of bin breaks, in tree time (root at 0)
timechop <- function(tree, numbins) {
  tipheights =  node.depth.edgelength(tree)[1:length(tree$tip.label)] 
  bk = quantile(tipheights, seq(from = 0, to=1,length.out = numbins+1))
  group = cut(tipheights, bk, include.lowest = T,labels = 1:numbins)
  return(list(group=group, times=bk))
}

#' @param tree a phylogenetic tree 
#' @param tipbins output of hte function timechop
#' @param aln - output of read.fasta; the sequences for the tips
#' @param gdtree - tree, in units of genetic distance (not timed), though this is not essential
#' @return genetic distance profile over time (a data frame) 
gendistprofiles <- function(tree, tipbins, aln=NULL,scale=TRUE,gdtree=NULL, mode="alignment") { 

  if (mode == "alignment")  sf=ifelse(scale==TRUE, length(aln[[1]]), 1)
  gps=unique(tipbins$group)
  num.bins=length(gps)
  cum.diversity=0*(1:num.bins) 
  now.diversity=0*(1:num.bins)
  num.tips=0*(1:num.bins)

  isDiverseAnyTime = rep(F, length(aln[[1]]))

  for (k in 1:length(gps)) {
    thesetips = tree$tip.label[which(tipbins$group == k)]
   if (mode == "alignment") {
     myind = vapply(thesetips, function(x) grep(x, names(aln)), FUN.VALUE = 1)
    mycharmat=matrix("hi",nrow = length(myind), ncol = length(aln[[myind[1]]]))
    for (n in 1:length(myind)) mycharmat[n,]=aln[[myind[n]]] # store sequence as character matrix
    nposs=vapply(1:ncol(mycharmat), function(x) { 
      myvals=unique(mycharmat[,x])
      return(length(myvals[which(myvals != "-")]))
      }, FUN.VALUE = 1)
    isDiverse = vapply(1:ncol(mycharmat), function(x) { !hasminus[x] & (nposs[x]>1)}, FUN.VALUE = T)
    isDiverseAnyTime= (isDiverse | isDiverseAnyTime)
    cum.diversity[k] = sum(isDiverseAnyTime)/sf 
    now.diversity[k]=sum(isDiverse)/sf
    num.tips[k]=length(myind)
   }
    if (mode == "tree") {
        # sum of path lengths in subtree of just these tips
      thistree = drop.tip(tree, which( ! (tree$tip.label %in% thesetips)))
      now.diversity[k] = sum(thistree$edge.length)
      allthesetips = tree$tip.label[which(as.numeric(tipbins$group) <= k)]
      cum.diversity[k]=sum( drop.tip(tree,which( !(tree$tip.label %in% allthesetips)))$edge.length)
      num.tips[k]=length(thesetips)
    }
  }
  return(data.frame(cum.diversity=cum.diversity, now.diversity=now.diversity,
              time=tipbins$times[2:length(tipbins$times)]))
} 

#' @param timetree a phylogenetic tree (in units of time; the clade of focus) 
#' @param gdtree a phylogenetic tree in units of genetic distance, which contains all tips of the timetree
#' @return lineages through genetic distance space, using the part of gdtree containing the tips in timetree
lttgendist <- function(timetree, gdtree) {
 mytree = drop.tip(gdtree, which( !(gdtree$tip.label  %in% timetree$tip.label))) # tree of only my tips
 return(ltt.plot.coords(mytree,backward = F))
}
  
#' @param tree a phylogenetic tree
#' @param gdtree a phylogenetic tree in units of genetic distance, which contains all tips of the timetree
#' @return a 2-column matrix, genetic dist from root and number of splitting events before that distance. 
divttgendist <- function(timetree, gdtree,maxdate=NULL) {
  mytree = drop.tip(gdtree, which( !(gdtree$tip.label  %in% timetree$tip.label))) # tree of only my tips
  x = ltt.plot.coords(tree,backward = F)
  # take all the negative increments out and just have the births
  toRemove = 1+which(diff(x[,2])<0)
  x=x[-toRemove,]
  x[,2]=1:nrow(x)
  if (!is.null(maxdate)) {x[,1]= x[,1]+ maxdate - max(node.depth.edgelength(tree))} 
  return(data.frame(time=x[,1],Ngdist=x[,2]))
}



