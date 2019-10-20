# Script to generate features from NA sequences  
# Pak Yu and Quang Nguyen
# Updated 10/20

library(ape)
library(phangorn)
library(seqinr)
library(stringr)
library(e1071)

#Utility functions  
source("./NA_mining_utils.R")


# Load files  
HA_tree <- read.tree(file = "./flutreeH3N2_HA.tree") # tree file 
HA_dat <- readRDS(file = "./processed_dat/H3N2_HA_processed_nonepitope.rds") # processed file 
alignment_file <- read.fasta(file = "./Data/alignH3N2NA.fa") # NA aligned sequences 

alignment_key <- read.csv(file = "./Data/alignH3N2NA_labels_clades_dates.csv",
                          stringsAsFactors = F, header = F)

map_file <- read.table(file = "./HA_NA_map.txt", sep = ',')
HA_subtrees <- HA_dat$subtrees  


# getting all phylogenetic trees from subtrees
tree_clades <- tree_trim(HA_tree, HA_subtrees)

feat_mat <- t(sapply(1:length(tree_clades), function(x){
  print(x)
  get_subtree_NA_feat(tree_clades[[x]], alignment_file = alignment_file,
                      alignment_key = alignment_key, map_file = map_file)
}))
colnames(feat_mat) <- c("min_hamming_distance", "max_hamming_distance", "median_hamming_distance", 
                        "min_gc_content", "max_gc_content", "mean_gc_content")
feat_mat <- as.data.frame(feat_mat)
feat_mat <- cbind(as.numeric(names(tree_clades)), feat_mat)
colnames(feat_mat)[1] <- "root_node_id"

write.csv(feat_mat, file = "./processed_dat/NA_hamming_gc.csv")



