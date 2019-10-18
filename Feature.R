source("~/Tree_Statistics.R")
library(ape)
library(phangorn)
library(apTreeshape)
library(igraph)
library(phyloTop)
library(moments)
library(devtools)
#devtools::install_github('Leonardini/treeCentrality')
library(treeCentrality)
library(seqinr)

#==========================================================================================
#read the tree
load("flutree.Rdata")
plot(flutree,show.tip.label = FALSE)
#read the protein data
Pdata=read.fasta("~/PH3N2.fa",seqtype="AA")

L=sapply(names(Pdata), function(x) strsplit(x,"_"))
Lab=numeric()
for(i in L){
  Lab=c(Lab,i[1])
}
names(Pdata)=Lab
#read the aux table with the RNA and protein labels
MappingData=read.csv("RNAP.csv",sep= ",",header=T,stringsAsFactors=FALSE)
MappingData=MappingData[,2:ncol(MappingData)]
treeinds=match(flutree$tip.label, MappingData[,1])
all(flutree$tip.label == MappingData[treeinds,1])
MappingData=MappingData[treeinds,]
#extract the clades
myclades=getClades2(flutree)
Final_trimmedClades=myclades$trimclades[myclades$rejected==0]
Final_trimmedClades_root=as.numeric(names(which(myclades$rejected==0)))
allHeights=node.depth.edgelength(flutree); max(allHeights)  #
hdata=data.frame(tiplab=flutree$tip.label, height=allHeights[1:length(flutree$tip.label)])
#remove the clades that do not have enough time to grows
res=numeric()
for(i in 1:length(Final_trimmedClades_root)){
  print(i)
  if( max(allHeights)-allHeights[Final_trimmedClades_root[i]]<=1.4){res=c(res,i)}
}

Final_root=Final_trimmedClades_root[-res] 

# Ntip1=length(tr1$tip.label)
# Ntip2=length(Final_trimmedClades[[348]])
# tr=drop.tip(tr1,setdiff(tr1$tip.label,Auxdata[Final_trimmedClades[[1]],2]), trim.internal = TRUE)
Clade_featurs=function(flutree,tr1,tr,root,MappingData,hdata, Pdata){
  Ntip1=length(tr1$tip.label)
  Ntip2=length(tr$tip.label)
  features=numeric(34)
  Epitop=getEpitopeDist(match(tr$tip.label,hdata$tiplab), MappingData,hdata, Pdata, pastperiod=5, D0=14)
  features[1]=root
  features[2]=Ntip1
  features[3]=Ntip2
  features[4]=sackin(as.treeshape(tr),"pda")
  features[5]=colless(as.treeshape(tr),"pda")
  features[6]=var(node.depth.edgelength(tr)[1:length(tr$tip.label)])
  features[7]=computeI2(tr)
  features[8]=computeB1(tr)
  features[9]=computeB2(tr)
  features[10]=avgLadder(tr, normalise = TRUE)
  features[11]=ILnumber(tr, normalise = TRUE)
  features[12]=pitchforks(tr, normalise = TRUE)
  features[13]=maxHeight(tr, normalise = TRUE)
  features[14]=computeMaxWidth(tr)
  features[15]=computeDelW(tr)
  features[16]=computeStairs1(tr)
  features[17]=computeStairs2(tr)
  features[19]=computeCherries(tr, DOUBLE = FALSE)
  features[20]=BS(tr)
  features[21]=descinm(tr,m=2)$W
  features[22]=getstattest(tr)$W
  features[23]=skewness(i_bl(tr))
  features[24]=kurtosis(i_bl(tr))
  features[25]=tips_pairwise_distance(tr)
  features[26]=tips_pairwise_distance_max(tr)
  features[27:31]=computeNetworkStats (tr, weight = FALSE, meanpath = FALSE, maxOnly = TRUE)
  features[32]=median(Epitop)
  features[33]=max(Epitop)
  features[34]=mean(Epitop)
  return(features)
}
#=================extract the features=================================================
Data_Clade=matrix(0,length(Final_trimmedClades),34)

for(i in match(Final_root,Final_trimmedClades_root)){
  print(i)
  tr1=extract.clade(flutree,Final_trimmedClades_root[i])
  tr=drop.tip(tr1,setdiff(tr1$tip.label,hdata$tiplab[Final_trimmedClades[[i]]]), trim.internal = TRUE)
  root=Final_trimmedClades_root[i]
  Data_Clade[i,]=Clade_featurs(flutree,tr1,tr,root,MappingData,hdata, Pdata)
}

write.csv(Data_Clade,"~/Data.csv")
#=====================================================================================
