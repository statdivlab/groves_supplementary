# load required packages
# install.packages("ape")
library(ape)
# install.packages("phangorn")
library(phangorn)
# install.packages("tidyverse")
library(tidyverse)
# install.packages("latex2exp")
library(latex2exp)
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("ggtree")
library(ggtree)
# if (!require("remotes", quietly = TRUE))
#   install.packages("remotes")
# remotes::install_github("statdivlab/groves")
library(groves)
# sessionInfo()
# 09/22/22 run uses:
# groves_0.0.0.9000, ggtree_3.4.2, tidyverse_1.3.2, ape_5.6-2, phangorn_2.10.0, latex2exp_0.9.5
# data was generated with GTDB R202

# read in gene names 
gene_names <- stringr::str_sub(list.files("data/strep/gene_trees/"), 1, -5)

# start by computing the log map visualization for all the gene trees and the phylogenomic
# tree created from concatenating all 196 genes

# make a vector of paths to .txt files of all trees to go into the plot 
paths <- c(paste0("data/strep/gene_trees/", list.files("data/strep/gene_trees")),
           "data/strep/phylogenomic_trees/concat_tree.txt")
# compute log map coordinates for all trees 
lm_res <- compute_logmap(tree_paths = paths,
                         tree_names = c(gene_names, "phylogenomic"))
# we can see that the base tree was chosen to be the phylogenomic tree
print(lm_res$base_lab)
plot_res <- plot_logmap(vectors = lm_res$vectors, phylogenomic = length(paths),
            title = "Streptococcus Trees", tree_names = c(gene_names, "phylogenomic"),
            phylogenomic_name = "$\\bar{T}_p^{full}$")
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz.jpeg")

# compute proportion of variance explained by each component 
pca_res <- stats::prcomp(lm_res$vectors)
pca_res$sdev^2/sum(pca_res$sdev^2)

# let's plot the ribosomal trees in the plot
ribs <- grep("Ribosom", gene_names)
ribs <- ifelse(1:length(paths) %in% ribs, "ribosomal", "other")
plot_res <- plot_logmap(vectors = lm_res$vectors, phylogenomic = length(paths), group = ribs, alpha = 0.8,
            title = "Streptococcus Trees", tree_names = c(gene_names, "phylogenomic"),
            phylogenomic_name = "$\\bar{T}_p^{full}$")
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz_rib.jpeg")

# add the phylogenomic tree from the ribosomal set 
paths <- c(paths, "data/strep/phylogenomic_trees/concat_tree_rib.txt")
new_vectors <- add_vector(new_tree_path = paths[length(paths)],
                          base_path = "data/strep/gene_trees/Cpn60_TCP1.txt",
                          vectors = lm_res$vectors,
                          new_name = "ribosomal phylogenomic")
plot_res <- plot_logmap(vectors = new_vectors, phylogenomic = nrow(new_vectors) - 1, 
            other_tree = nrow(new_vectors), 
            phylogenomic_name = "$\\bar{T}_p^{full}$", 
            other_name = "$\\bar{T}_p^{rib}$", ignore_in_pca = nrow(new_vectors),
            title = "Streptococcus Trees", group = c(ribs, "other"),
            tree_names = c(gene_names, "phylogenomic", "ribosomal phylogenomic"),
            alpha = 0.8)
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz_full.jpeg")
# add names for outlying trees in PC1 and PC2 
plot_res <- plot_logmap(vectors = new_vectors, phylogenomic = nrow(new_vectors) - 1, 
            other_tree = nrow(new_vectors), 
            phylogenomic_name = "$\\bar{T}_p^{full}$", 
            other_name = "$\\bar{T}_p^{rib}$", ignore_in_pca = nrow(new_vectors),
            title = "Streptococcus Trees", group = c(ribs, "other"),
            tree_names = c(gene_names, "phylogenomic", "ribosomal phylogenomic"),
            alpha = 0.8,
            trees_to_label = c("DUF3270", "DUF1934"))
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz_full_labeled.jpeg")

# compare the two phylogenomic trees
# compare BHV distances between all trees 
merge_txt(paths, "data/strep/all_trees.txt")
all_trees <- read.tree("data/strep/all_trees.txt")
bhv_dists <- compute_geodesic("data/strep/all_trees.txt")
bhv_dist_vec <- bhv_dists[lower.tri(bhv_dists)]
phy_bhv <- bhv_dists[nrow(bhv_dists), nrow(bhv_dists) - 1]
mean(bhv_dist_vec >= phy_bhv)
# 75% of BHV distances between trees are larger than the distance between the two 
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
# 38% of squared branch sum differences are larger than the squared branch sum
# difference between the two phylogenomic trees

# color the log map plot by sum of branch lengths
plot_res <- plot_logmap(vectors = new_vectors, phylogenomic = nrow(new_vectors) - 1, 
            other_tree = nrow(new_vectors), 
            phylogenomic_name = "$\\bar{T}_p^{full}$", 
            other_name = "$\\bar{T}_p^{rib}$", ignore_in_pca = nrow(new_vectors),
            title = "Streptococcus Trees", group = unlist(branch_sums),
            tree_names = c(gene_names, "phylogenomic", "ribosomal phylogenomic"),
            alpha = 0.8)
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz_sum_branches.jpeg")
# principal component 1 seems to be governed by branch lengths

# re-scale gene trees to divide each branch by the sum of branch lengths
rescaled_trees <- standardize_branches(all_trees)
new_paths <- write_trees_txt(rescaled_trees, "data/strep/rescaled_trees", 
                             c(gene_names, "phylo", "ribo_phylo"))
rescaled_lm_res <- compute_logmap(tree_paths = new_paths,
                         tree_names = c(gene_names, "phylogenomic", "ribosomal phylogenomic"))
rescaled_lm_res$base_lab
# the new base tree is the phylogenomic tree
plot_res <- plot_logmap(vectors = rescaled_lm_res$vectors, phylogenomic = length(gene_names) + 1, 
            other_tree = length(gene_names) + 2, 
            phylogenomic_name = "$\\bar{T}_p^{full}$", 
            other_name = "$\\bar{T}_p^{rib}$",
            title = "Standardized Streptococcus Trees", group = c(ribs, "other"),
            tree_names = c(gene_names, "phylogenomic", "ribosomal phylogenomic"),
            alpha = 0.8)
# add names of trees that have large magnitudes of PC1 or PC2 
plot_res <- plot_logmap(vectors = rescaled_lm_res$vectors, phylogenomic = length(gene_names) + 1, 
            other_tree = length(gene_names) + 2, 
            phylogenomic_name = "$\\bar{T}_p^{full}$", 
            other_name = "$\\bar{T}_p^{rib}$",
            title = "Standardized Streptococcus Trees", group = c(ribs, "other"),
            tree_names = c(gene_names, "phylogenomic", "ribosomal phylogenomic"),
            alpha = 0.8, 
            trees_to_label = c("eIF-1a", "Ribosomal_S17", "Ribosomal_S21"))
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/standardized_lm_plot.jpeg")

# check for trees that are not binary
binary_trees <- check_binary(all_trees)
c(gene_names, "concat", "ribosomal_concat")[!binary_trees]

# compute Frechet mean tree using Megan Owen's SturmMean.jar implementation 
merge_txt(paths[1:(length(paths) - 2)], "data/strep/gene_trees.txt")
fm_output <- system2("java", args = c("-jar", "scripts/jar/SturmMean.jar", "-o", 
                        "data/strep/frechet_mean_tree.txt", 
                        "data/strep/gene_trees.txt"), stdout = T)
frechet_mean <- ape::read.tree(text = fm_output[12])
ape::is.binary(frechet_mean)

# compute summary statistics of phylogenomic tree branch lengths
full_phylo <- ape::read.tree("data/strep/phylogenomic_trees/concat_tree.txt")
rib_phylo <- ape::read.tree("data/strep/phylogenomic_trees/concat_tree_rib.txt")
summary(full_phylo$edge.length)
summary(rib_phylo$edge.length)
