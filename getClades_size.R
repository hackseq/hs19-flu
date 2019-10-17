
#' @param rt A tree in phylo format; should be rooted
#' @param minSize The minimum total size of a clade
#' @param maxSize The maximum size for a  clade.
#' @return A list of nodeids, the tips in the node's clade, the clade sizes, and the branch length from which the clade descends
#' @examples
#' myclades <- pickClades(tree)
pickClades <- function(rt, minSize=100,maxSize=200) {
  nTips=length(rt$tip.label)
  myroot=nTips+1
  nnodes=nTips-1
  nodeids=(nTips+1):(nTips+nnodes) # of internal nodes
  dfsall=dfs(graph(rt$edge),root=nTips+1)$order # igraph. perhaps not necessary
  dfsnodes=as.vector(dfsall[dfsall>nTips]) # 
  
  ChildMatrix=t(sapply(nodeids, function(x) rt$edge[rt$edge[,1]==x, 2])) ##has two dimensions but by default, only the second is getting used
  rownames(ChildMatrix)=nodeids # each row lists the two children of a node, names are node ids
  
  NodeDescOf <- function(node) { mydescs=ChildMatrix[as.character(node),]
  return(mydescs[mydescs > nTips])} ##node means always nodeids
  
  allD=allDescendants(rt) ## index with number of descendant . requires phangorn 
  allCladeSizes=sapply(nodeids,function(x) sum(allD[[x]] <= nTips)) 
  
  rejectFlag=-1+0*(1:nnodes) ##list of length nnodes with only 0 in it
  names(rejectFlag)=nodeids
  for (k in 2:length(dfsnodes)){ 
    ii=dfsnodes[k]
    if (rejectFlag[ii-nTips] == -1) {
      # if clade is too big, reject and move to descendant, which will be next on the list as we did a dfs
      if (allCladeSizes[ii-nTips] > maxSize) { rejectFlag[ii-nTips]=1 }
      
      # if clade is too small, reject and reject all descendants
      if (allCladeSizes[ii-nTips] < minSize) { rejectFlag[ii-nTips]=1; iiDescs=allD[[ii]][allD[[ii]]>nTips]; 
      rejectFlag[as.character(iiDescs)]=1 
      }
      # if a clade is just right, accept but reject all descendants
      if (allCladeSizes[ii-nTips] >= minSize & allCladeSizes[ii-nTips] <= maxSize) {
        rejectFlag[ii-nTips]=0 # accept
        iiDescs=allD[[ii]][allD[[ii]]>nTips]; # but reject the descendants
        rejectFlag[as.character(iiDescs)]=1  
      }
    }
  } # end main loop
  rejectFlag[as.character(myroot)]=1; 
  goodnodes = nTips+ which(rejectFlag == 0) # should be the NODE index since I added nTips
  goodsizes = allCladeSizes[which(rejectFlag==0)]
  goodclades = sapply(goodnodes,function(x) allD[[x]][allD[[x]] <=nTips])
  
  # want branch lengths from which these nodes descend
  brlens=sapply(goodnodes,function(x) { rt$edge.length[which(rt$edge[,2]==x)]})
  return(list(nodes=goodnodes,clades=goodclades, sizes=goodsizes,brlens=brlens))
}

