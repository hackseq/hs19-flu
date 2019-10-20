# Script to extract tree features (similar to Maryam's script) from NA trees for all HA subtrees 
# Quang Nguyen
# Updated 10/20/19

library(ape)
library(phangorn)
library(apTreeshape)
library(igraph)
library(phyloTop)
library(moments)
library(treeCentrality)

source('./NA_mining_utils.R')
source('./Tree_Statistics.R')

NA_tree <- read.tree(file = "./clade_assignments/trees/flutreeH3N2_NA.nwk")
HA_subtrees <- readRDS(file = "./processed_dat/H3N2_HA_processed.rds")$subtree
HA_tree <- read.tree(file = "./clade_assignments/trees/flutreeH3N2_HA.nwk")
map_file <- read.table(file = "./HA_NA_map.txt", sep = ",")
tree_clades <- tree_trim(HA_tree, HA_subtrees)


NA_matched_tree_features <- c()
for (i in 1:length(tree_clades)){
  if (i %% 5 == 0){
    print(i)
  }
  tree <- tree_clades[[i]]
  tips <- as.vector(na.omit(sapply(tree$tip.label, HA_to_NA, map_file = map_file)))
  if (length(tips) <= 1){
    NA_matched_tree_features <- rbind(NA_matched_tree_features,rep(NA,27)) # if there are no tips or that there are only 1 match, no tree to be have 
  } else {
    tips <- tips[which(tips %in% NA_tree$tip.label)] # for now remove those matches that somehow is not in the list of tips from files
    NA_subtree <- keep.tip(NA_tree, tip = tips)
    NA_matched_tree_features <- rbind(NA_matched_tree_features, Clade_features_modified(tree))
  }
}

colnames(NA_matched_tree_features) <- c("sackin","colless","Variance","I2","B1","B2","avgLadder","ILnumber","pitchforks",
                                        "maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries","BS","descinm","getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist",
                                        "diameter", "WienerIndex", "betweenness", "closeness", "eigenvector")

NA_matched_tree_features <- as.data.frame(NA_matched_tree_features)
NA_matched_tree_features$Nodeid <- as.numeric(names(tree_clades))
NA_matched_tree_features <- NA_matched_tree_features[,c(ncol(NA_matched_tree_features),1:ncol(NA_matched_tree_features) - 1)]

write.csv(NA_matched_tree_features, file = "./processed_dat/subtree_NA_tree_features.csv")
