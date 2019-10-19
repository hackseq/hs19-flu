#need ape package for makeNodeLabel function
library(ape)

#load Rdata files
flutree2018_5<-load.Rdata2("../../flutree2018-5.Rdata",path=getwd())
flutreeB_HA<-load.Rdata2("../../flutreeB_HA.Rdata",path=getwd())
flutreeB_NA<-load.Rdata2("../../flutreeB_NA.Rdata",path=getwd())
flutreeH1N1_HA<-load.Rdata2("../../flutreeH1N1_HA.Rdata",path=getwd())
flutreeH3N2_HA<-load.Rdata2("../../flutreeH3N2_HA.Rdata",path=getwd())
flutreeH3N2_NA<-load.Rdata2("../../lutreeH3N2_NA.Rdata",path=getwd())

#apply makeNodeLabels
makeNodeLabel(flutree2018_5)
makeNodeLabel(flutreeB_HA)
makeNodeLabel(flutreeB_NA)
makeNodeLabel(flutreeH1N1_HA)
makeNodeLabel(flutreeH3N2_HA)
makeNodeLabel(flutreeH3N2_NA)

#write newick files
write.tree(flutree2018_5, file="clade_assignments/trees/flutree2018_5.nwk")
write.tree(flutreeB_HA, file="clade_assignments/trees/flutreeB_HA.nwk")
write.tree(flutreeB_NA, file="clade_assignments/trees/flutreeB_NA.nwk")
write.tree(flutreeH1N1_HA, file="clade_assignments/trees/flutreeH1N1_HA.nwk")
write.tree(flutreeH3N2_HA, file="clade_assignments/trees/flutreeH3N2_HA.nwk")
write.tree(flutreeH3N2_NA, file="clade_assignments/trees/flutreeH3N2_NA.nwk")


