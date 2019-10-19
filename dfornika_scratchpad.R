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


#' @param labels_clades_dates_path path to three-column csv file containing labels, clade names and dates
load_labels_clades_dates <- function(labels_clades_dates_path) {
  labels_clades_dates <- read_csv(
    labels_clades_dates_path,
           col_names = FALSE, 
           col_types = cols(X1=col_character(), X2=col_character(), X3=col_date("%Y-%m-%d"))
    )
  names(labels_clades_dates) <- c('label','clade_name','date')
  labels_clades_dates <- as_tbl_time(labels_clades_dates, index=date)
}


#' @param tree phylogenetic tree
#' @param labels_dates data frame with fields 'label' and 'date'
#' @param truncation_date date after which all tips will be dropped from the tree
truncate_tree_by_date <- function(tree, labels_dates, truncation_date) {
  labels_to_exclude <- labels_dates %>%
    filter(date > truncation_date) %>%
    pull(label)
  ladderize(drop.tip(tree, labels_to_exclude))
}


# alignB_labels_clades_dates <- load_labels_clades_dates("Data/alignB_labels_clades_dates.csv")
# alignBNA_labels_clades_dates <- load_labels_clades_dates("Data/alignBNA_labels_clades_dates.csv")
# alignH1N1_labels_clades_dates <- load_labels_clades_dates("Data/alignH1N1_labels_clades_dates.csv")
# alignH1N1NA_labels_clades_dates <- load_labels_clades_dates("Data/alignH1N1NA_labels_clades_dates.csv")
alignH3N2_labels_clades_dates <- load_labels_clades_dates("Data/alignH3N2_labels_clades_dates.csv")
# alignH3N2NA_labels_clades_dates <- load_labels_clades_dates("Data/alignH3N2NA_labels_clades_dates.csv")
# PH3N2_labels_clades_dates <- load_labels_clades_dates("Data/PH3N2_labels_clades_dates.csv")

alignH3N2_labels_clades_dates  %>% pull(date) %>% max()

flutreeH3N2_HA <- read.tree("flutreeH3N2_HA.tree")


flutreeH3N2_HA_2014 <- truncate_tree_by_date(flutreeH3N2_HA, alignH3N2_labels_clades_dates, as.Date("2014-07-31"))
flutreeH3N2_HA_2015 <- truncate_tree_by_date(flutreeH3N2_HA, alignH3N2_labels_clades_dates, as.Date("2015-07-31"))
flutreeH3N2_HA_2016 <- truncate_tree_by_date(flutreeH3N2_HA, alignH3N2_labels_clades_dates, as.Date("2016-07-31"))
flutreeH3N2_HA_2017 <- truncate_tree_by_date(flutreeH3N2_HA, alignH3N2_labels_clades_dates, as.Date("2017-07-31"))


#trees <- c(flutreeH3N2_HA_2014, flutreeH3N2_HA_2015, flutreeH3N2_HA_2016, flutreeH3N2_HA_2017)
#class(trees) <- "multiPhylo"
#ggtree(trees) + facet_wrap(~.id) + theme_tree2()

clades_2016 = pickClades(flutreeH3N2_HA_2016, 200, 400)
clades_2017 = pickClades(flutreeH3N2_HA_2017, 200, 400)

tree_2016_plot <-  ggtree(flutreeH3N2_HA_2016, mrsd="2016-07-31")
for (i in seq_along(clades_2016$nodes)) {
  tree_2016_plot <- tree_2016_plot + geom_cladelabel(node=clades_2016$nodes[i], label=clades_2016$nodes[i], align=T) 
}
tree_2016_plot <- tree_2016_plot + theme_tree2()
tree_2016_plot

tree_2017_plot <- ggtree(flutreeH3N2_HA_2017, mrsd="2017-07-31")
for (i in seq_along(clades_2017$nodes)) {
  tree_2017_plot <- tree_2017_plot + geom_cladelabel(node=clades_2017$nodes[i], label=clades_2017$nodes[i], align=T) 
}
tree_2017_plot <- tree_2017_plot + theme_tree2()
tree_2017_plot

ratios = c()

for (i in seq_along(clades_2016$clades)) {
  clade_2016_labels <- flutreeH3N2_HA_2016$tip.label[clades_2016$clades[[i]]]
  clade_2016_mrca <- findMRCA(flutreeH3N2_HA_2016, clade_2016_labels)
  clade_2017_mrca <- findMRCA(flutreeH3N2_HA_2017, clade_2016_labels)

  clade_2016_mrca_descendants <- getDescendants(flutreeH3N2_HA_2016, clade_2016_mrca)
  clade_2017_mrca_descendants <- getDescendants(flutreeH3N2_HA_2017, clade_2017_mrca)
  ratio <- length(clade_2017_mrca_descendants) / length(clade_2016_mrca_descendants)
  ratios <- append(ratios, ratio)
}

