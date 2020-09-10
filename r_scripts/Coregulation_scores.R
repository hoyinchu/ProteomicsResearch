## This script creates co-regulation scores for protein pairs in ProteomeHD
## Hyperparameters of treeClust and TOM's sigmoid funcion were previously optimised

# Load required libraries
library(data.table); library(treeClust); library(WGCNA)

# Set seed to make results reproducible
set.seed(42)

# Make sure the working directory is set to the folder downloaded from Github. If it isn't, set it with setwd()
getwd()

## Load ProteomeHD
pQTL_protein <- fread("pQTL/pQTL_protein_converted.csv")                             # Read in fast using data.table
pQTL_protein <- data.frame(pQTL_protein, row.names = "uniprotswissprot")    # Convert to data.frame; use protein IDs as rownames
pQTL_protein <- pQTL_protein[, grep("GM", colnames(pQTL_protein)) ]                 # Keep only columns with SILAC ratios

## Obtain treeClust distances
tc_distances <- treeClust.dist( pQTL_protein,
                                d.num = 2,
                                verbose = TRUE,
                                rcontrol = rpart.control(cp = 0.105),
                                control = treeClust.control(serule = 1.8) )

## Turn the distance matrix into a similarity matrix
tc_sim_symm <- 1-as.matrix(tc_distances)

## Calculate the adjacency matrix using the sigmoid function
adj_mat <- sigmoidAdjacencyFunction(tc_sim_symm, mu = 0.91, alpha = 37)

## Get the Topological Overlap Matrix
tom_sim <- TOMsimilarity( adj_mat, TOMDenom = "mean" )
colnames(tom_sim) <- colnames(tc_sim_symm)
row.names(tom_sim) <- colnames(tc_sim_symm)


## Test if network is scale-free
connectivity <- colSums(adj_mat, na.rm = TRUE) - 1
png("ScaleFreeNess.png")
scaleFreePlot(connectivity)
dev.off()


## Turn similarity matrices into long-format, remove duplicates and merge into final table
tc_sim_dt <- as.data.table( reshape2::melt( tc_sim_symm )) 
tc_sim_dt <- tc_sim_dt[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tc_sim = value ) ]
tc_sim_dt <- tc_sim_dt[ Protein_1 > Protein_2 ]                  # Removes duplicate pairs (incl. self-comparisons)

adj_mat_dt <- as.data.table( reshape2::melt( adj_mat )) 
adj_mat_dt <- adj_mat_dt[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tc_adj = value ) ]
adj_mat_dt <- adj_mat_dt[ Protein_1 > Protein_2 ]                # Removes duplicate pairs (incl. self-comparisons)

tc_tom_dt <- as.data.table( reshape2::melt( tom_sim )) 
tc_tom_dt <- tc_tom_dt[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tc_tom = value ) ]
tc_tom_dt <- tc_tom_dt[ Protein_1 > Protein_2 ]                  # Removes duplicate pairs (incl. self-comparisons)

tc_dt <- merge( tc_sim_dt, adj_mat_dt, by = c("Protein_1", "Protein_2"))
tc_dt <- merge( tc_dt, tc_tom_dt,      by = c("Protein_1", "Protein_2"))

## Write out the combined result file
fwrite(tc_dt, "pQTL_treeClust_similarities.csv")

## Write out simplified result file containing only the final scores
tc_dt_final <- tc_dt[, .(Protein_1, Protein_2, coregulation_score = tc_tom)]
fwrite(tc_dt_final, "pQTL_coregulation_scores.csv")



