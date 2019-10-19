 getManyPlots <- function(myclade, alignment, gdtree) {
  ## compute things 
  
  # things that are NOT cumulative 
  lttinfo=as.data.frame(ltt.plot.coords(myclade,backward = F))
  glttd=lttgendist(timetree = myclade, gdtree = gdtree) # ltt through gen dist 
  
  # a couple functions make things both cumulative and non-cumulative
  gg=gendistprofiles(myclade, timechop(myclade, 30),gdtree = gdtree,mode = "tree")
  ggt = gendistprofiles(myclade, timechop(myclade,30), aln=alignment, mode="alignment")

 # things that ARE cumulative 
  d=divtt(myclade) ; # tips through time 
  dg=divttgendist(timetree=myclade,gdtree=gdtree,maxdate=NULL) # tips through gen dist  space 
  
  ## plot things
  # plot the tree itself - requires ggtree 
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
  
  dev.new()
  multiplot(p1,n1,n2,n3,n4, ncol=1)
  dev.new()
  multiplot(p1,c1,c2,c3,c4, ncol=1)
  
 }
 