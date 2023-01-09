# load required packages
# install.packages("ape")
library(ape)
# install.packages("tidyverse")
library(tidyverse)
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("ggtree")
library(ggtree)
# if (!require("remotes", quietly = TRUE))
#   install.packages("remotes")
# remotes::install_github("statdivlab/groves")
library(groves)
# sessionInfo()
# 09/29/22 run uses:
# groves_0.0.0.9000, ggtree_3.4.2, tidyverse_1.3.2, ape_5.6-2
# data was generated with GTDB R202

# read in gene names 
gene_names <- stringr::str_sub(list.files("data/prevotella/gene_trees/"), 1, -5)

# start by computing the log map visualization for all the gene trees and the phylogenomic
# tree created from concatenating all 63 genes

# make a vector of paths to .txt files of all trees to go into the plot 
paths <- c(paste0("data/prevotella/gene_trees/", list.files("data/prevotella/gene_trees")),
           "data/prevotella/phylogenomic_trees/concat_tree.txt")
# compute log map coordinates for all trees 
lm_res <- compute_logmap(tree_paths = paths,
                         tree_names = c(gene_names, "phylogenomic"))
# we can see that the base tree was chosen to be the phylogenomic tree
# this means that this tree has the lowest mean distance from all of the 
# other trees out of the full set 
print(lm_res$base_lab)
plot_res <- plot_logmap(vectors = lm_res$vectors, phylogenomic = 64,
            title = "Prevotella Trees", tree_names = c(gene_names, "phylogenomic"),
            phylogenomic_name = "$\\bar{T}_p^{full}$") 
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz.jpeg")

# compute proportion of variance explained by each component 
pca_res <- stats::prcomp(lm_res$vectors)
pca_res$sdev^2/sum(pca_res$sdev^2)
 
# add the concatenated tree from the reduced set (without the three outliers)
new_vectors <- add_vector(new_tree_path = "data/prevotella/phylogenomic_trees/red_concat.txt",
                          base_path = "data/prevotella/phylogenomic_trees/concat_tree.txt",
                          vectors = lm_res$vectors,
                          new_name = "reduced phylogenomic")
plot_res <- plot_logmap(vectors = new_vectors, phylogenomic = 64, other_tree = 65, 
            phylogenomic_name = "$\\bar{T}_p^{full}$",
            other_name = "$\\bar{T}_p^{reduced}$", ignore_in_pca = 65,
            title = "Prevotella Trees", 
            tree_names = c(gene_names, "phylogenomic", "reduced phylogenomic"))
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_two_phy.jpeg")

# next, look at individual gene trees 
# we can see 3 clear outliers in the log map visualization 
# if we look at the output of plot_logmap, we can see that the two outliers in PC1
# are DMRL_synthase and GTP_cyclohydroI, and the outlier in PC2 is BacA 
# let's add these labels to our plot 
plot_res <- plot_logmap(vectors = new_vectors, phylogenomic = 64, other_tree = 65, 
            phylogenomic_name = "$\\bar{T}_p^{full}$",
            other_name = "$\\bar{T}_p^{reduced}$", ignore_in_pca = 65,
            title = "Prevotella Trees", 
            tree_names = c(gene_names, "phylogenomic", "reduced phylogenomic"),
            trees_to_label = c("BacA", "DMRL_synthase", "GTP_cyclohydroI"))
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_labeled.jpeg", width = 7, height = 3.5, dpi = 300)

# plot trees 
# phylogenomic tree 
phylogenomic <- read.tree("data/prevotella/phylogenomic_trees/concat_tree.txt")
# rename labels from accession numbers to numbers 
match_df <- data.frame(oldlab = phylogenomic$tip.label,
                       newlab = 1:length(phylogenomic$tip.label))
phylogenomic_num <- rename_labs(match_df, phylogenomic)
phylogenomic_num_mid <- phytools::midpoint.root(phylogenomic_num)
ggtree(phylogenomic_num_mid) + geom_tiplab(size = 2) + 
  xlim(c(0, 0.55)) + ggtitle("Phylogenomic Tree") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/phylogenomic_tree.jpeg")
# BacA gene tree 
BacA <- read.tree("data/prevotella/gene_trees/BacA.txt")
BacA_num <- rename_labs(match_df, BacA)
BacA_num_mid <- phytools::midpoint.root(BacA_num)
ggtree(BacA_num_mid) + geom_tiplab(size = 2) + 
  xlim(c(0, 3.3)) + ggtitle("BacA Tree") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/BacA_tree.jpeg")
# DMRL_synthase
DMRL_synthase <- read.tree("data/prevotella/gene_trees/DMRL_synthase.txt")
DMRL_synthase_num <- rename_labs(match_df, DMRL_synthase)
DMRL_synthase_num_mid <- phytools::midpoint.root(DMRL_synthase_num)
ggtree(DMRL_synthase_num_mid) + geom_tiplab(size = 2) + 
  xlim(c(0, 6)) + ggtitle("DMRL_synthase Tree") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/DMRL_synthase.jpeg")
# GTP_cyclohydroI 
GTP_cyclohydroI <- read.tree("data/prevotella/gene_trees/GTP_cyclohydroI.txt")
GTP_cyclohydroI_num <- rename_labs(match_df, GTP_cyclohydroI)
GTP_cyclohydroI_num_mid <- phytools::midpoint.root(GTP_cyclohydroI_num)
ggtree(GTP_cyclohydroI_num_mid) + geom_tiplab(size = 2) + 
  xlim(c(0, 2.8)) + ggtitle("GTP_cyclohydroI Tree") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/GTP_cyclohydroI.jpeg")

# Plot gene tree support values 
merge_txt(paths[1:63], "data/prevotella/gene_trees.txt")
gene_trees <- read.tree("data/prevotella/gene_trees.txt")
support <- check_gene_support(main_tree = phylogenomic,
                  trees = gene_trees,
                  rooted = FALSE) 
plot_support(phylogenomic_num, support, color_branch = TRUE,
             title = "Phylogenomic Tree Gene Support", xlim_max = 1.1)
ggsave("figures/prevotella/gene_tree_support.jpeg")

# try with midpoint rooted tree 
gene_trees_mid <- phangorn::midpoint(rename_labs(match_df, gene_trees))
support_mid <- check_gene_support(main_tree = phylogenomic_num_mid,
                                  trees = gene_trees_mid, 
                                  rooted = TRUE)
plot_support(phylogenomic_num_mid, support_mid, color_branch = TRUE,
             title = "Phylogenomic Tree Gene Support", xlim_max = 0.54)
ggsave("figures/prevotella/gene_tree_support.jpeg")

# Plot bootstrap support values 
phylo_bs <- phangorn::midpoint(rename_labs(match_df, 
                                           ape::read.tree("data/prevotella/phylogenomic_trees/concat_IQTREEE.treefile")))
boot <- as.numeric(phylo_bs$node.label)/100
boot[1] <- 1
plot_support(phylo_bs, boot, color_branch = TRUE, 
             title = "Phylogenomic Bootstrap Support", xlim_max = 0.54,
             support_type = "boot")
ggsave("figures/prevotella/bootstrap_support.jpeg")

# make table to match tip number to accession number to taxonomic data 
taxonomic_info <- read.table("data/prevotella/NCBI_genomes_summary_info.tsv", 
                             sep = '\t', header = TRUE)
tax_index <- which(taxonomic_info$input_accession %in% match_df$old)
tax_table <- select(taxonomic_info[tax_index, ], 
                    input_accession, organism_name, infraspecific_name)
tax_table <- arrange(tax_table, input_accession)
num_acc <- arrange(match_df, oldlab)
tax_table$number <- num_acc$newlab  
tax_table <- arrange(tax_table, number)
tax_table <- select(tax_table, -number)
# install.packages("xtable")
library(xtable)
xtable(tax_table, type = "latex")

# check for trees that are not binary
binary_trees <- check_binary(gene_trees)
gene_names[!binary_trees]

# compute Frechet mean tree using Megan Owen's SturmMean.jar implementation 
fm_output <- system2("java", args = c("-jar", "scripts/jar/SturmMean.jar", "-o", 
                                      "data/prevotella/frechet_mean_tree.txt", 
                                      "data/prevotella/gene_trees.txt"), stdout = T)
frechet_mean <- ape::read.tree(text = fm_output[12])
ape::is.binary(frechet_mean)

# compare two phylogenomic trees 
# compare BHV distances between all trees 
merge_txt(c(paths, "data/prevotella/phylogenomic_trees/red_concat.txt"), 
          "data/prevotella/all_trees.txt")
all_trees <- read.tree("data/prevotella/all_trees.txt")
bhv_dists <- compute_geodesic("data/prevotella/all_trees.txt")
bhv_dist_vec <- bhv_dists[lower.tri(bhv_dists)]
phy_bhv <- bhv_dists[nrow(bhv_dists), nrow(bhv_dists) - 1]
mean(bhv_dist_vec >= phy_bhv)
# 100% of BHV distances between trees are larger than the distance between the two 
# phylogenomic trees
# compare RF distances between all trees
rf_dists <- as.matrix(phangorn::RF.dist(all_trees))
rf_dist_vec <- rf_dists[lower.tri(rf_dists)]
phy_rf <- rf_dists[nrow(rf_dists), nrow(rf_dists) - 1]
mean(rf_dist_vec >= phy_rf)
# 99.9% of RF distances between trees are larger than the distance between the
# two phylogenomic trees
# compare branch length sums between all trees
branch_sums <- lapply(all_trees, function(tree) {sum(tree$edge.length)})
branch_sum_dists <- as.matrix(stats::dist(data.frame(x = unlist(branch_sums)), 
                                          method = "euclidean"))
branch_sum_dist_vec <- branch_sum_dists[lower.tri(branch_sum_dists)]
phy_branch_sum <- branch_sum_dists[nrow(branch_sum_dists), nrow(branch_sum_dists) - 1]
mean(branch_sum_dist_vec >= phy_branch_sum)
# 99% of squared branch sum differences are larger than the squared branch sum
# difference between the two phylogenomic trees

# check relative weight of branch to tip 51 in PC1 
max_branch <- max(DMRL_synthase$edge.length)
lm_comp <- which(lm_res$vectors == max_branch, arr.ind = TRUE)[2]
pca_res <- stats::prcomp(lm_res$vectors, rank. = 2)
weight_51 <- pca_res$rotation[lm_comp, 1]
weight_51/sum(abs(pca_res$rotation[, 1]))
