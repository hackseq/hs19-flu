#!/usr/bin/env Rscript

library(igraph)
library(phangorn)
library(ape)
library(readr)
library(dplyr)
library(tibbletime)
library(treeio)
library(ggplot2)
library(phytools)

source("getClades_size.R")
source("larger_scale_features.R")



# alignB_labels_clades_dates <- load_labels_clades_dates("Data/alignB_labels_clades_dates.csv")
# alignBNA_labels_clades_dates <- load_labels_clades_dates("Data/alignBNA_labels_clades_dates.csv")
# alignH1N1_labels_clades_dates <- load_labels_clades_dates("Data/alignH1N1_labels_clades_dates.csv")
# alignH1N1NA_labels_clades_dates <- load_labels_clades_dates("Data/alignH1N1NA_labels_clades_dates.csv")
alignH3N2_labels_clades_dates <- load_labels_clades_dates("Data/alignH3N2_labels_clades_dates.csv")
# alignH3N2NA_labels_clades_dates <- load_labels_clades_dates("Data/alignH3N2NA_labels_clades_dates.csv")
# PH3N2_labels_clades_dates <- load_labels_clades_dates("Data/PH3N2_labels_clades_dates.csv")

alignH3N2_labels_clades_dates  %>% pull(date) %>% max()

flutreeH3N2_HA <- read.tree("flutreeH3N2_HA.tree")


flutreeH3N2_HA_2014 <- truncate_tree_by_date(flutreeH3N2_HA, alignH3N2_labels_clades_dates, as.Date("2014-02-28"))
flutreeH3N2_HA_2015 <- truncate_tree_by_date(flutreeH3N2_HA, alignH3N2_labels_clades_dates, as.Date("2015-02-28"))
flutreeH3N2_HA_2016 <- truncate_tree_by_date(flutreeH3N2_HA, alignH3N2_labels_clades_dates, as.Date("2016-02-28"))
flutreeH3N2_HA_2017 <- truncate_tree_by_date(flutreeH3N2_HA, alignH3N2_labels_clades_dates, as.Date("2017-02-28"))


#trees <- c(flutreeH3N2_HA_2014, flutreeH3N2_HA_2015, flutreeH3N2_HA_2016, flutreeH3N2_HA_2017)
#class(trees) <- "multiPhylo"
#ggtree(trees) + facet_wrap(~.id) + theme_tree2()

clades_2016 = pickClades(flutreeH3N2_HA_2016, 200, 400)
clades_2017 = pickClades(flutreeH3N2_HA_2017, 200, 400)

tree_2016_plot <-  ggtree(flutreeH3N2_HA_2016, mrsd="2016-02-28")
for (i in seq_along(clades_2016$nodes)) {
  tree_2016_plot <- tree_2016_plot + geom_cladelabel(node=clades_2016$nodes[i], label=clades_2016$nodes[i], align=T) 
}
tree_2016_plot <- tree_2016_plot + theme_tree2()
tree_2016_plot

tree_2017_plot <- ggtree(flutreeH3N2_HA_2017, mrsd="2017-02-28")
for (i in seq_along(clades_2017$nodes)) {
  tree_2017_plot <- tree_2017_plot + geom_cladelabel(node=clades_2017$nodes[i], label=clades_2017$nodes[i], align=T) 
}
tree_2017_plot <- tree_2017_plot + theme_tree2()
tree_2017_plot


growth_ratios_for_clades(flutreeH3N2_HA_2016, flutreeH3N2_HA_2017, 200, 400)



