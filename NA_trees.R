# Script to generate features from NA sequences  
# Pak Yu and Quang Nguyen
# Updated 10/20

library(ape)
library(phangorn)
library(seqinr)
library(stringr)
library(e1071)

# Load files  
HA_tree <- read.tree(file = "./flutreeH3N2_HA.tree") # tree file 
HA_dat <- readRDS(file = "./processed_dat/H3N2_HA_processed_nonepitope.rds") # processed file 
alignment_file <- read.fasta(file = "./Data/alignH3N2NA.fa") # NA aligned sequences 

alignment_key <- read.csv(file = "./Data/alignH3N2NA_labels_clades_dates.csv",
                          stringsAsFactors = F, header = F)

map_file <- read.table(file = "./HA_NA_map.txt", sep = ',')
HA_subtrees <- HA_dat$subtrees  

# Extracting subtrees in tree form with pruning using Maryam's code
#' @param flutree Tree full tree file  
#' @param myclades Output from getClades2() function
#' @return List of phylogenetic trees with names as the root node 
tree_trim <- function(flutree, myclades){
  # getting clades without rejection
  Final_trimmedClades=myclades$trimclades[myclades$rejected==0]
  # get root ids from final_trimmed_clades
  Final_trimmedClades_root=as.numeric(names(which(myclades$rejected==0)))
  # Get branch lengths  
  allHeights=node.depth.edgelength(flutree); max(allHeights)
  # grab all tips and their lengths 
  hdata=data.frame(tiplab=flutree$tip.label, height=allHeights[1:length(flutree$tip.label)])
  # remove subtrees without proper growth 
  print(max(allHeights)-allHeights[Final_trimmedClades_root[314]] <= 1.4) 
  res=numeric()
  for(i in 1:length(Final_trimmedClades_root)){
    if((max(allHeights)-allHeights[Final_trimmedClades_root[i]]) <= 1.4){
      res=c(res,i)
    }
  }
  Final_root=Final_trimmedClades_root[-res] 
  result <- list()
  for (i in match(Final_root, Final_trimmedClades_root)){
    if (i %% 5 == 0){print(i)}
    tr1 <- extract.clade(flutree,Final_trimmedClades_root[i]) # get all from root 
    tr <- drop.tip(tr1,setdiff(tr1$tip.label,hdata$tiplab[Final_trimmedClades[[i]]]), trim.internal = TRUE)
    root <- Final_trimmedClades_root[i]
    result[[as.character(root)]] <- tr
  }
  return(result)
}

#' Function to map HA sample names to NA sample names  
#' @param HA HA index, as a string (can get from)
HA_to_NA <- function(HA, map_file){
  idx <- which(map_file[,2] == paste0(">",HA,"_A"))
  output <- stringr::str_match(map_file[idx,1], regex(">(.*?)_A"))[2]
  return(output)
}

#' Extracting features from NA sequences using subtrees  
#' @param subtree The subtree as a phylo object
#' @param alignment_file The alignment file in fasta format and loaded as read.fasta from the seqinr package
#' @param alignment_key The clade key csv file for HA alignments with sample names as the first column. 
#'                      Sequence of samples should be the same as that of the HA alignment fasta file  
#' @param map_file The mapping file of samples between HA and NA, generated courtesy of Pak. HA second column, NA first column
#' @return A vector of length 6 with min, max and median/mean of hamming distance and GC content respectively
get_subtree_NA_feat<- function(subtree, alignment_file, alignment_key, map_file){
  labs <- subtree$tip.label
  # get index matching to alignment_file 
  match_idx <- numeric(length(labs))
  for (i in 1:length(labs)){
    NA_lab <- HA_to_NA(HA = labs[i], map_file = map_file)
    if (!is.na(NA_lab)){
      match_idx[i] <- which(alignment_key[,1] == HA_to_NA(HA = labs[i],map_file = map_file))
    }
  }
  if (length(which(match_idx != 0)) == 0){ # means that there is no match between HA and NA
    return(rep(NA,6))
  } else {
    match_idx <- match_idx[match_idx != 0] # only grab non-zero 
    # generating a matrix of sequence where each letter is it's own element in a vector 
    seq_mat <- c()
    for (i in 1:length(match_idx)){
      seq_mat <- rbind(seq_mat, getSequence(alignment_file[match_idx[i]])[[1]])
    }
    if (nrow(seq_mat) == 1){ # there is only one match 
      gc_content <- rep(seqinr::GC(seq_mat[1,]),3)
      distance <- c(0,0,0)
      output <- c(distance, gc_content)
      return(output)
    } else {
      gc_content <- t(apply(seq_mat, 1, function(x){
        seqinr::GC(x)
      }))
      distance <- as.vector(e1071::hamming.distance(seq_mat))
      distance <- distance[distance != 0] # remove 0s since they are identical sequences non informative
      output <- c(min(distance), max(distance), median(distance), min(gc_content), max(gc_content), mean(gc_content))
      return(output)
    }
  }
}


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



