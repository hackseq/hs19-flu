# This script is for going from the Rdata files with the tree predictions to nwk files

library("ape")
load("flutreeB_NA.Rdata")
labeled <- makeNodeLabel(flutree)
write.tree(labeled, file="flutreeB_NA.nwk")

load("flutreeB_HA.Rdata")
labeled <- makeNodeLabel(flutree)
write.tree(labeled, file="flutreeB_HA.nwk")

load("flutreeB_NA.Rdata")
labeled <- makeNodeLabel(flutree)
write.tree(labeled, file="flutreeH1N1_HA.nwk")

load("flutreeB_NA.Rdata")
labeled <- makeNodeLabel(flutree)
write.tree(labeled, file="flutreeH3N2_HA.nwk")

load("flutreeB_NA.Rdata")
labeled <- makeNodeLabel(flutree)
write.tree(labeled, file="flutreeH3N2_NA.nwk")

load("flutree2018-5.Rdata")
labeled <- makeNodeLabel(flutree)
write.tree(labeled, file="flutreeH3N2_NA.nwk")