df=read.csv("Data",sep= ",",header=T,stringsAsFactors=FALSE)
names(df)=c("Clade","numberTipsClade","numberTipsTrimmed","sackin",
            "colless","Variance","I2","B1","B2","avgLadder","ILnumber","pitchforks",
            "maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries","BS","descinm","getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist",
            "diameter", "WienerIndex", "betweenness", "closeness", "eigenvector","MeadianEp","MaxEp","MeanEp")

allHeights=node.depth.edgelength(tree)
allD=allDescendants(tree) 
TimeFrame=3.4
#nodeids are the ids if the clades (the first column of the data)
nodeids=df$Clade
nTips=length(tree$tip.label)
# # tips within the TimeFrame for each node
allTrimmedClades = sapply(nodeids, function(x) {  myTipDes=allD[[x]][allD[[x]]<=nTips]
myTipTimes=allHeights[myTipDes] # here, would need something like allDates[myTipDesc]
return(myTipDes[myTipTimes <= allHeights[x]+TimeFrame]) })
# # sizes of trimmed clades 
allTrimmedSizes_3.4 = sapply(allTrimmedClades, function(x) length(allTrimmedClades[x]))
Labels=allTrimmedSizes_3.4/df$numberTipsTrimmed