# needs phangorn and igraph

#' @param rt A tree in phylo format; should be rooted and timed (or with lengths in units corresponding to the TimeFrame argument)
#' @param MinTotalSize The minimum total size of a clade
#' @param MinTrimSize The minimum size for a trimmed clade. These will be in the inputs for the learning
#' @param OverlapCutoff The maximum fraction of overlap two clades can have and both be accepted
#' @param TimeFrame The maximum time between a node and the last tip in the node's trimmed clade
#' @return A list of nodeids, the node's relevant parent (RP), the node's clade size (untrimmed), the node's 
#' trimmed size (trimsize), whether the node was rejected (rejected), and the tips in the node's trimmed clade (trimclades, a list)
#' @examples
#' myclades <- getclades(tree)
  getClades2 <- function(rt,MinTotalSize=10, MinTrimSize=10, OverlapCutoff=0.8, TimeFrame=0.5,truncate=NA) { 
    # set up
    nTips=length(rt$tip.label)
    myroot=nTips+1
    nnodes=nTips-1
    nodeids=(nTips+1):(nTips+nnodes)
    dfsall=dfs(graph(rt$edge),root=nTips+1)$order # igraph
    dfsnodes=as.vector(dfsall[dfsall>nTips]) # they are already in order but I don't know that this will be true for any tree
    # consider removing that
    RP=NA+(1:nnodes); names(RP)=nodeids
    ChildMatrix=t(sapply(nodeids, function(x) rt$edge[rt$edge[,1]==x, 2]))
    rownames(ChildMatrix)=nodeids # each row lists the two children of a node, names are node ids
    NodeDescOf <- function(node) { mydescs=ChildMatrix[as.character(node),]
    return(mydescs[mydescs > nTips])}
    
    RP[as.character(NodeDescOf(myroot))]=myroot # alternatively could subtract nTips from the node ids toget the row numbers
    # which might be much faster for the whole tree ; doing as.char now for clarity... hmm.
    
    # extra initialisation that solved problem of the rare node not getting an RP, in the staph data: 
    RP=unlist(sapply(nodeids, function(x) rt$edge[which(rt$edge[,2]==x),1]))
    RP <- c(myroot, RP)
    names(RP)=nodeids
    
    # compute heights 
    allHeights=node.depth.edgelength(rt); # heights from root to node, in units of branch length in the tree. assumes timed tree

    # need to define something like allDates, at least for the tips, giving the dates from the metadata, since these can't be 
    # computed from heights of nodes without a timed tree
    
    # compute descendants and clade sizes 
    allD=allDescendants(rt) 
    allCladeSizes=sapply(nodeids,function(x) sum(allD[[x]] <= nTips)) 
    
    # tips within the TimeFrame for each node
    allTrimmedClades = sapply(nodeids, function(x) {  myTipDes=allD[[x]][allD[[x]]<=nTips]
        myTipTimes=allHeights[myTipDes] # here, would need something like allDates[myTipDesc]
        return(myTipDes[myTipTimes <= allHeights[x]+TimeFrame]) }) # here, replace with [myTipTimes == 2006] 
    
    if (!is.na(truncate)) {
    allFullTruncClades=sapply(nodeids,function(x) { myTipDes=allD[[x]][allD[[x]] <=nTips]
          myTipTimes=allHeights[myTipDes]
          return(myTipDes[myTipTimes<= allHeights[x]+truncate]) })
    fulltruncs=sapply(nodeids,function(x) length(allFullTruncClades[[x-nTips]]))
    }
    
    # sizes of trimmed clades 
    allTrimmedSizes = sapply(nodeids, function(x) length(allTrimmedClades[[x-nTips]]))
    
    
    rejectFlag=0*(1:nnodes)
    names(rejectFlag)=nodeids

    
    # main loop
    for (k in 2:length(dfsnodes)) {
      ii = dfsnodes[k] 
      
      if (rejectFlag[as.character(ii)] != 1) {
        # does ii have an RP? 
        rpii=RP[as.character(ii)]
        
        # if size is small, reject the node and all its descendants
        if (allCladeSizes[ii-nTips] < MinTotalSize) { 
          rejectFlag[as.character(ii)]= 1 ; iiDescs=allD[[ii]][allD[[ii]]> nTips]; 
          rejectFlag[as.character(iiDescs)]=1 
        }
        # -- if size of Ci(T) too small but Ci is big enough, set flag for non-use of i,
        #  can still use i's descendants
        if (allTrimmedSizes[ii-nTips] < MinTrimSize & allCladeSizes[ii-nTips] >= MinTotalSize) { 
          rejectFlag[as.character(ii)]=1 # as.char or ii-nTips; same effect
          RP[as.character(NodeDescOf(ii))]= ii
        }
        
        # -- if size is big enough, check intersection of clade with relevant parent's clade 
        if (allTrimmedSizes[ii-nTips] >= MinTrimSize & allCladeSizes[ii-nTips] >= MinTotalSize) { 
          # check intersection 
          myintersect = intersect(allTrimmedClades[[rpii-nTips]], allTrimmedClades[[ii-nTips]]) 
          
          #     -- if overlap is "big", set flag for non-use and set relevant parent (RP) of i's children to RP of i. 
          # overlap is the portion of ii's trimmed clade that is contained in the parent's trimmed clade
          if (length(myintersect) > OverlapCutoff*allTrimmedSizes[[ii-nTips]] ) {
            rejectFlag[as.character(ii)]=1
            RP[as.character(NodeDescOf(ii))]=rpii
          } else {   # If overlap small - keep i, set RP of i's children to i; do not reject i
            RP[as.character(NodeDescOf(ii))]=ii
          } # end if - else on the intersection
        } # end if size is big enough
      } # end if not reject flag
    } # end main loop
    rejectFlag[as.character(myroot)]=1
return(list(nodes=nodeids,RP=RP, sizes= allCladeSizes,trimsize= allTrimmedSizes,
            rejected= rejectFlag,trimclades=allTrimmedClades,
            fullclades=ifelse(is.na(truncate),NA, allFullTruncClades), 
            fulltruncs=fulltruncs))
  }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        