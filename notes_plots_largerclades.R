library(seqinr)

source("getClades_size.R")
source("plotClades.R")
source("larger_scale_features.R")

timeinfo = read.csv("Data/H3N2_labels_clades_dates.csv", header = F); head(timeinfo)
aln=read.fasta("Data/H3N2.fa")
gdh3n2 = read.tree("rootedRAxML_H3N2.tree")

#############################################
# notes and testing my larger-scale features
#############################################

flu=read.tree("timetreeH3N2copy.nwk")
gdh3n2=read.tree("rootedRAxML_H3N2.tree")

clades2019= pickClades(flu, minSize = 300, maxSize = 500)

plotClades(flu, clades2019, show.tip.label = F, trysize = "big")
# my functions are timechop (sets time bins), divtt, ltt, gendistprofiles (diversity through time),
# lttgendist,  and divttgendist 
myclade = extract.clade(flu, node = clades2019$nodes[6])
plot(myclade,show.tip.label=F)


## compute things 

# things that are NOT cumulative 
lttinfo=as.data.frame(ltt.plot.coords(myclade,backward = F))
glttd=lttgendist(timetree = myclade, gdtree = gdh3n2) # ltt through gen dist 

# a couple functions make things both cum and non-cum 
gg=gendistprofiles(myclade, timechop(myclade, 30),gdtree = gdh3n2,mode = "tree")
ggt = gendistprofiles(myclade, timechop(myclade,30), aln=aln, mode="alignment")
# things that ARE cumulative 

# number of tips through time 
d=divtt(myclade) ;
dg=divttgendist(timetree=myclade,gdtree=gdh3n2,maxdate=NULL) # tips through gen dist 




p1=ggtree(myclade) 

## plot things that are NOT cumulative 
# lineages through time
n1=ggplot(data=lttinfo, aes(x=time,y=N))+geom_line()+ggtitle("Lineages through time") + 
  theme(plot.title = element_text(size = 8))
n2= ggplot(data.frame(glttd), aes(x=time, y=N))+geom_line()+
  ggtitle("Lineages through gen dist pseudo-time") + 
  theme(plot.title = element_text(size = 8))
n3=ggplot(data=gg, aes(x=time, y=now.diversity))+geom_line()+
  ggtitle("Gen. diversity in gd tree") + 
  theme(plot.title = element_text(size = 8))
n4=ggplot(data=ggt, aes(x=time, y=now.diversity))+geom_point()+
  ggtitle("Gen. div in alignment")  + 
  theme(plot.title = element_text(size = 8))

# things that are cumulative 
c1=ggplot(data=d, aes(x=time, y=Nsplits))+geom_line()+ ggitle("Num tips through time") + 
  theme(plot.title = element_text(size = 8)) # tips through time 
c2=ggplot(dg, aes(x=time,y=Ngdist))+geom_line()+
  ggtitle("Lineages through genetic distance pseudo-time") + 
  theme(plot.title = element_text(size = 8))
c3=ggplot(data=gg, aes(x=time, y=cum.diversity))+geom_line()+
  ggtitle("Cumulative gen diversity in gd tree") + 
  theme(plot.title = element_text(size = 8))
c4=ggplot(data=ggt, aes(x=time, y=cum.diversity))+geom_point()+
  ggtitle("Cumulative gen diversity in alignment") + 
  theme(plot.title = element_text(size = 8))

# source("multiplot.R") # available in scater package which my computer says isn't avail for R 3.4.4
# so instead I have read it fro mthe disk as I have used it before 
dev.set(2)
multiplot(p1,n1,n2,n3,n4, ncol=1)
dev.set(5)
multiplot(p1,c1,c2,c3,c4, ncol=1)


## have made all that into a function getManyPlots.R


flu2014 <- truncate_tree_by_date(flu, alignH3N2_labels_clades_dates, as.Date("2014-02-28"))
comp2014 <- truncate_tree_by_date(flu, alignH3N2_labels_clades_dates, as.Date("2015-05-30"))

flu2015 <- truncate_tree_by_date(flu, alignH3N2_labels_clades_dates, as.Date("2015-02-28"))
comp2015 <- truncate_tree_by_date(flu, alignH3N2_labels_clades_dates, as.Date("2016-05-30"))

flu2016 <- truncate_tree_by_date(flu, alignH3N2_labels_clades_dates, as.Date("2016-02-28"))
comp2016 <- truncate_tree_by_date(flu, alignH3N2_labels_clades_dates, as.Date("2017-05-30"))

flu2017 <- truncate_tree_by_date(flu, alignH3N2_labels_clades_dates, as.Date("2017-02-28"))
comp2017 <- truncate_tree_by_date(flu, alignH3N2_labels_clades_dates, as.Date("2018-05-30"))

growth2014
growth2015
growth2016
growth2017


getManyPlots(extract.clade(flu2016, c2016$nodes[15]), aln, gdh3n2)






