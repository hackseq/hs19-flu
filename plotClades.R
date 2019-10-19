
plotClades <- function(tree, cladelist, show.tip.label=F,trysize="big") {
   # get some colours - get length(myclades$nodes) of them 
  jColors <- c('chartreuse3', 'cornflowerblue', 'darkgoldenrod1', 'peachpuff3',
               'mediumorchid2', 'turquoise3', 'wheat4', 'slategray2')
  mycols=sample(jColors, length(cladelist$nodes), replace=TRUE)
  # define the groups 
  # define cols
  if (trysize=="small") {
    groups=cladelist$clades 
 ecol= edge.color(tree, groups,col=mycols) 
 plot(tree, edge.color = ecol,no.margin = TRUE, edge.width = 2, show.tip.label= show.tip.label)
} else {
    # for bigger trees
    tipcols=rep("grey", length(tree$edge[,1]))
 for (k in 1:length(groups)) { 
   thisgroup=groups[[k]]
  for (n in 1:length(thisgroup)) { 
    thisedge=which(tree$edge[,2] == thisgroup[n])
    tipcols[thisedge] = mycols[k]
  }
 }
  plot(tree, edge.color=tipcols, show.tip.label=F)
}
  }
