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
source("")


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

#' @param tree1 phylogenetic tree, truncated at an earlier time than tree2
#' @param tree2 phylogenetic tree, truncated at an later time than tree1
#' @param min_clade_size minimum clade size
#' @param max_clade_size max clade size
growth_ratios_for_clades <- function(tree1, tree2, min_clade_size, max_clade_size){
  tree1_clades <- pickClades(tree1, min_clade_size, max_clade_size)
  ratios <- c()
  for (i in seq_along(tree1_clades$clades)) {
    tree1_clade_tip_labels <- tree1$tip.label[tree1_clades$clades[[i]]]
    tree1_clade_mrca <- findMRCA(tree1, tree1_clade_tip_labels)
    tree2_clade_mrca <- findMRCA(tree2, tree1_clade_tip_labels)
    tree1_clade_mrca_descendants <- getDescendants(tree1, tree1_clade_mrca)
    tree2_clade_mrca_descendants <- getDescendants(tree2, tree2_clade_mrca)
    ratio <- length(tree2_clade_mrca_descendants) / length(tree1_clade_mrca_descendants)
    ratios <- append(ratios, ratio)
  }
  data.frame("clade" = tree1_clades$nodes, "growth_ratio" = ratios)
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


growth_ratios_for_clades(flutreeH3N2_HA_2016, flutreeH3N2_HA_2017, 200, 400)



