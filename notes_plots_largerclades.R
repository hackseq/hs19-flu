
timeinfo = read.csv("Data/alignH3N2_labels_clades_dates.csv", header = F); head(timeinfo)

library(seqinr)
aln=read.fasta("Data/alignH3N2.fa")
gdh3n2 = read.tree("rootedRAxML_H3N2.tree")

#############################################
# notes and testing my larger-scale features
#############################################

flu=read.tree("flutreeH3N2_HA.tree")
gdh3n2=read.tree("rootedRAxML_H3N2.tree")

clades2019= pickClades(flu, minSize = 300, maxSize = 500)

plotClades(flu, clades2019, show.tip.label = F)
# my functions are timechop (sets time bins), divtt, ltt, gendistprofiles (diversity through time),
# lttgendist,  and divttgendist 

# number of tips through time 
myclade = extract.clade(flu, node = clades2019$nodes[6])
plot(myclade,show.tip.label=F)

d=divtt(myclade) ;
ggplot(data=d, aes(x=time, y=N))+geom_line()

# lineages through time
lttinfo=as.data.frame(ltt.plot.coords(myclade,backward = F))
ggplot(data=lttinfo, aes(x=time,y=N))+geom_line()

# genetic diversity over time as measured in the UNTIMED phylo tree 
# requires the rooted raxml (but untimed) phylogenetic tree
# eg gdh3n2 = read.tree("rootedRAxML_H3N2.tree")
gg=gendistprofiles(myclade, timechop(myclade, 30),gdtree = gdh3n2,mode = "tree")
plot(gg$time, gg$cum.diversity,'l')

# genetic diversity over time as measured in the *alignment*- requires alignment file
# eg: aln = read.fasta("Data/alignH3N2.fa" which in turn requires the seqinr package
ggt = gendistprofiles(myclade, timechop(myclade,30), aln=aln, mode="alignment")
ggplot(data=ggt, aes(x=time, y=cum.diversity))+geom_point()
ggplot(data=ggt, aes(x=time, y=now.diversity))+geom_point()

# like lineages through time, but lineages through genetic distance 
dg=divttgendist(timetree=myclade,gdtree=gdh3n2,maxdate=NULL)
ggplot(dg, aes(x=time,y=N))+geom_line()













