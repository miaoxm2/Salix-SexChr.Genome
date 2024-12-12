### cluster gene trees

library(ape)
# set path to the tree folder
tree_folder <- "/16-genes/tree"  


# extract all tree files 
tree_files <- list.files(path = tree_folder, pattern = "*.treefile", full.names = TRUE)

# read all tree files
trees <- lapply(tree_files, read.tree)

#modify tip label, remove gene name, only remain species name and chr
modify_tree_tips <- function(tree) {
  tree$tip.label <- sapply(tree$tip.label, function(label) {
    modified_label <- gsub("\\G.*", "", label)
    modified_label1 <- gsub("\\_chr15.*", "", modified_label)
    return(modified_label1)
  })
  return(tree)
}
modified_trees <- lapply(trees, modify_tree_tips)

# set rooted tip
root_tip <- "Potri.015" 

# reroot all trees
modified_trees <- lapply(modified_trees, function(tree) root(tree, outgroup = root_tip, resolve.root = TRUE))


# create an empty matrix
n_trees <- length(modified_trees)
distance_matrix <- matrix(0, nrow = n_trees, ncol = n_trees)

# calculate topo distance between trees
for (i in 1:(n_trees-1)) {
  for (j in (i + 1):n_trees) {
    # No.i -- No.j
    distance <- dist.topo(modified_trees[[i]], modified_trees[[j]])
    # store distance
    distance_matrix[i, j] <- distance
    distance_matrix[j, i] <- distance  
  }
}


# hierarchy cluster and plot 
hc <- hclust(as.dist(distance_matrix), method = "average")
plot(hc, main = "Hierarchical Clustering of Trees")


library(ggplot2)
# choose cluster level
k <- 4
clusters <- cutree(hc, k)  

# create dataframe of tree name and number of cluster
cluster_results <- data.frame(Tree = tree_files, Cluster = clusters)
# check cluster=N
cluster_results[which(cluster_results$Cluster==3),]
# output the dataframe

# plot consensus tree for each cluster level in one page
consensus_trees <- list()


par(mfrow = c(3, 2))  
for (cluster_num in unique(cluster_results$Cluster)) {
  
  cluster_trees <- modified_trees[cluster_results$Cluster == cluster_num]
  
  consensus_tree <- root(consensus(cluster_trees, p = 0.8), outgroup = root_tip, resolve.root = TRUE)


  plot(consensus_tree,cex = 2, main = paste0(length(cluster_trees),"genes - Consensus Tree for Cluster", cluster_num))
# nodelabels(round(consensus_tree$node.label,2))  # add support number. should be numeric. round()cannot work after reroot
}


pdf("sco_trees.4sp6cluster.pdf", width = 10, height = 10)


dev.off()
write.csv(cluster_results, file = "tree_clusters2.csv", row.names = FALSE)


