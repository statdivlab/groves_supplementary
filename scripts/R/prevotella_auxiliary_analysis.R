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
                        title = "Phylogenomic base tree", tree_names = c(gene_names, "phylogenomic"),
                        phylogenomic_name = "$\\bar{T}_p^{full}$") 
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_phylogenomic_base.png")

# compare this to the visualization with the log map with different base trees 
# rank trees by mean squared distance from other trees 
merge_txt(paths, "data/prevotella/all_trees.txt")
all_trees <- read.tree("data/prevotella/all_trees.txt")
bhv_dists <- compute_geodesic("data/prevotella/all_trees.txt")
sq_dists <- bhv_dists^2
diag(sq_dists) <- NA
mean_dists <- rowMeans(sq_dists, na.rm = TRUE)
rank_df <- data.frame(rank = rank(mean_dists), mean_dist = mean_dists, 
                      name = c(gene_names, "phylogenomic")) %>%
  arrange(rank)
head(rank_df)
# we can see that the phylogenomic tree has the lowest mean squared distance from 
# the other trees, followed by the trees for genes LYTB, Voltage_CLC, YicC_N,
# Acyltransf_2, and DUF4924 

# log map visualization with LYTB as a base tree 
LYTB_lm_res <- compute_logmap(tree_paths = paths,
                         tree_names = c(gene_names, "phylogenomic"),
                         base_lab = "LYTB")
LYTB_plot_res <- plot_logmap(vectors = LYTB_lm_res$vectors, phylogenomic = 64,
                        title = "LYTB base tree", tree_names = c(gene_names, "phylogenomic"),
                        phylogenomic_name = "$\\bar{T}_p^{full}$",
                        other_tree = which(gene_names == "LYTB"),
                        other_name = "LYTB") 
LYTB_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_LYTB_base.png")

# log map visualization with Voltage_CLC as a base tree 
Voltage_CLC_lm_res <- compute_logmap(tree_paths = paths,
                              tree_names = c(gene_names, "phylogenomic"),
                              base_lab = "Voltage_CLC")
Voltage_CLC_plot_res <- plot_logmap(vectors = Voltage_CLC_lm_res$vectors, phylogenomic = 64,
                             title = "Voltage_CLC base tree", tree_names = c(gene_names, "phylogenomic"),
                             phylogenomic_name = "$\\bar{T}_p^{full}$",
                             other_tree = which(gene_names == "Voltage_CLC"),
                             other_name = "Voltage_CLC") 
Voltage_CLC_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_Voltage_CLC_base.png")

# log map visualization with YicC_N as a base tree 
YicC_N_lm_res <- compute_logmap(tree_paths = paths,
                                     tree_names = c(gene_names, "phylogenomic"),
                                     base_lab = "YicC_N")
YicC_N_plot_res <- plot_logmap(vectors = YicC_N_lm_res$vectors, phylogenomic = 64,
                                    title = "YicC_N base tree", tree_names = c(gene_names, "phylogenomic"),
                                    phylogenomic_name = "$\\bar{T}_p^{full}$",
                                    other_tree = which(gene_names == "YicC_N"),
                                    other_name = "YicC_N") 
YicC_N_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_YicC_N_base.png")

# log map visualization with Acyltransf_2  as a base tree 
Acyltransf_2_lm_res <- compute_logmap(tree_paths = paths,
                                tree_names = c(gene_names, "phylogenomic"),
                                base_lab = "Acyltransf_2")
Acyltransf_2_plot_res <- plot_logmap(vectors = Acyltransf_2_lm_res$vectors, phylogenomic = 64,
                               title = "Acyltransf_2 base tree", tree_names = c(gene_names, "phylogenomic"),
                               phylogenomic_name = "$\\bar{T}_p^{full}$",
                               other_tree = which(gene_names == "Acyltransf_2"),
                               other_name = "Acyltransf_2") 
Acyltransf_2_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_Acyltransf_2_base.png")

# log map visualization with DUF4924  as a base tree 
DUF4924_lm_res <- compute_logmap(tree_paths = paths,
                                      tree_names = c(gene_names, "phylogenomic"),
                                      base_lab = "DUF4924")
DUF4924_plot_res <- plot_logmap(vectors = DUF4924_lm_res$vectors, phylogenomic = 64,
                                     title = "DUF4924 base tree", tree_names = c(gene_names, "phylogenomic"),
                                     phylogenomic_name = "$\\bar{T}_p^{full}$",
                                     other_tree = which(gene_names == "DUF4924"),
                                     other_name = "DUF4924") 
DUF4924_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_DUF4924_base.png")

# look at results with an outlier as the base tree
# log map visualization with BacA as base tree 
BacA_lm_res <- compute_logmap(tree_paths = paths,
                                 tree_names = c(gene_names, "phylogenomic"),
                                 base_lab = "BacA")
BacA_plot_res <- plot_logmap(vectors = BacA_lm_res$vectors, phylogenomic = 64,
                                title = "BacA base tree", tree_names = c(gene_names, "phylogenomic"),
                                phylogenomic_name = "$\\bar{T}_p^{full}$",
                                other_tree = which(gene_names == "BacA"),
                                other_name = "BacA") 
BacA_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_BacA_base.png")

# log map visualization with DMRL_synthase as base tree 
DMRL_synthase_lm_res <- compute_logmap(tree_paths = paths,
                              tree_names = c(gene_names, "phylogenomic"),
                              base_lab = "DMRL_synthase")
DMRL_synthase_plot_res <- plot_logmap(vectors = DMRL_synthase_lm_res$vectors, phylogenomic = 64,
                             title = "DMRL_synthase base tree", tree_names = c(gene_names, "phylogenomic"),
                             phylogenomic_name = "$\\bar{T}_p^{full}$",
                             other_tree = which(gene_names == "DMRL_synthase"),
                             other_name = "DMRL_synthase") 
DMRL_synthase_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_DMRL_synthase_base.png")

# log map visualization with GTP_cyclohydroI as base tree 
GTP_cyclohydroI_lm_res <- compute_logmap(tree_paths = paths,
                                       tree_names = c(gene_names, "phylogenomic"),
                                       base_lab = "GTP_cyclohydroI")
GTP_cyclohydroI_plot_res <- plot_logmap(vectors = GTP_cyclohydroI_lm_res$vectors, phylogenomic = 64,
                                      title = "GTP_cyclohydroI base tree", tree_names = c(gene_names, "phylogenomic"),
                                      phylogenomic_name = "$\\bar{T}_p^{full}$",
                                      other_tree = which(gene_names == "GTP_cyclohydroI"),
                                      other_name = "GTP_cyclohydroI") 
GTP_cyclohydroI_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_GTP_cyclohydroI_base.png")

# compute RF distances (measure of topology only distance)
RF_dists <- as.matrix(phangorn::RF.dist(all_trees))
diag(RF_dists) <- NA
mean_dists <- rowMeans(RF_dists, na.rm = TRUE)
RF_rank_df <- data.frame(rank = rank(mean_dists), mean_dist = mean_dists, 
                      name = c(gene_names, "phylogenomic")) %>%
  arrange(rank)
head(RF_rank_df)
tail(RF_rank_df)

# look at base tree as tree with largest mean RF distance from other trees
# log map visualization with dUTPase as base tree 
dUTPase_lm_res <- compute_logmap(tree_paths = paths,
                                         tree_names = c(gene_names, "phylogenomic"),
                                         base_lab = "dUTPase")
dUTPase_plot_res <- plot_logmap(vectors = dUTPase_lm_res$vectors, phylogenomic = 64,
                                        title = "dUTPase base tree", tree_names = c(gene_names, "phylogenomic"),
                                        phylogenomic_name = "$\\bar{T}_p^{full}$",
                                        other_tree = which(gene_names == "dUTPase"),
                                        other_name = "dUTPase") 
dUTPase_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_dUTPase_base.png")

# log map visualization with Tyr_Deacylase as base tree 
Tyr_Deacylase_lm_res <- compute_logmap(tree_paths = paths,
                                 tree_names = c(gene_names, "phylogenomic"),
                                 base_lab = "Tyr_Deacylase")
Tyr_Deacylase_plot_res <- plot_logmap(vectors = Tyr_Deacylase_lm_res$vectors, phylogenomic = 64,
                                title = "Tyr_Deacylase base tree", tree_names = c(gene_names, "phylogenomic"),
                                phylogenomic_name = "$\\bar{T}_p^{full}$",
                                other_tree = which(gene_names == "Tyr_Deacylase"),
                                other_name = "Tyr_Deacylase") 
Tyr_Deacylase_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_Tyr_Deacylase_base.png")

# log map visualization with DHquinase_II as base tree 
DHquinase_II_lm_res <- compute_logmap(tree_paths = paths,
                                       tree_names = c(gene_names, "phylogenomic"),
                                       base_lab = "DHquinase_II")
DHquinase_II_plot_res <- plot_logmap(vectors = DHquinase_II_lm_res$vectors, phylogenomic = 64,
                                      title = "DHquinase_II base tree", tree_names = c(gene_names, "phylogenomic"),
                                      phylogenomic_name = "$\\bar{T}_p^{full}$",
                                      other_tree = which(gene_names == "DHquinase_II"),
                                      other_name = "DHquinase_II") 
DHquinase_II_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_DHquinase_II_base.png")

# also investigate Frechet mean tree
fm_output <- system2("java", args = c("-jar", "scripts/jar/SturmMean.jar", "-o", 
                                      "data/prevotella/frechet_mean_tree.txt", 
                                      "data/prevotella/gene_trees.txt"), stdout = T)
frechet_mean <- ape::read.tree(text = fm_output[12])
ape::is.binary(frechet_mean)
# distance between phylogenomic tree and frechet mean tree
merge_txt(c(paths, "data/prevotella/frechet_mean_tree.txt"), 
          "data/prevotella/phylo_and_frechet.txt")
phylo_frechet_dist <- compute_geodesic("data/prevotella/phylo_and_frechet.txt")
dist_vec <- phylo_frechet_dist[lower.tri(phylo_frechet_dist)]
phy_bhv <- phylo_frechet_dist[64, 65]
mean(dist_vec < phy_bhv)

# run MDS on the same data
MDS_res <- compute_MDS(dist_matrix = bhv_dists,
                         tree_names = c(gene_names, "phylogenomic"))
MDS_plot <- plot_MDS(df = MDS_res$df, phylogenomic = 64,
  title = "MDS of Prevotella trees", tree_names = c(gene_names, "phylogenomic"),
  phylogenomic_name = "$\\bar{T}_p^{full}$") 
MDS_plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/MDS_plot.png")

# use a randomly generated tree as the base tree
# seed of 1
set.seed(1)
rand_tree1 <- rtree(length(all_trees[[1]]$tip.label), rooted = FALSE)
rand_tree1$tip.label <- all_trees[[1]]$tip.label
write.tree(rand_tree1, "data/prevotella/rand_tree1.txt")
rand_tree1_lm_res <- compute_logmap(tree_paths = c(paths, "data/prevotella/rand_tree1.txt"),
                              tree_names = c(gene_names, "phylogenomic", "random"),
                              base_lab = "random")
rand_tree1_plot_res <- plot_logmap(vectors = rand_tree1_lm_res$vectors,
                                   phylogenomic = 64,
                                   phylogenomic_name = "$\\bar{T}_p^{full}$",
                                   title = "Random base tree", 
                                   tree_names = c(gene_names, "phylogenomic", "random"),
                                   other_tree = 65,
                                   other_name = "random") 
rand_tree1_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_rand_base1.png")
RF_dists <- as.matrix(phangorn::RF.dist(c(all_trees, rand_tree1)))
diag(RF_dists) <- NA
rowMeans(RF_dists, na.rm = TRUE)
# seed of 2
set.seed(2)
rand_tree2 <- rtree(length(all_trees[[1]]$tip.label), rooted = FALSE)
rand_tree2$tip.label <- all_trees[[1]]$tip.label
write.tree(rand_tree2, "data/prevotella/rand_tree2.txt")
rand_tree2_lm_res <- compute_logmap(tree_paths = c(paths, "data/prevotella/rand_tree2.txt"),
                                    tree_names = c(gene_names, "phylogenomic", "random"),
                                    base_lab = "random")
rand_tree2_plot_res <- plot_logmap(vectors = rand_tree2_lm_res$vectors,
                                   phylogenomic = 64,
                                   phylogenomic_name = "$\\bar{T}_p^{full}$",
                                   title = "Random base tree", 
                                   tree_names = c(gene_names, "phylogenomic", "random"),
                                   other_tree = 65,
                                   other_name = "random") 
rand_tree2_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_rand_base2.png")
RF_dists <- as.matrix(phangorn::RF.dist(c(all_trees, rand_tree2)))
diag(RF_dists) <- NA
rowMeans(RF_dists, na.rm = TRUE)
# seed of 3
set.seed(3)
rand_tree3 <- rtree(length(all_trees[[1]]$tip.label), rooted = FALSE)
rand_tree3$tip.label <- all_trees[[1]]$tip.label
write.tree(rand_tree3, "data/prevotella/rand_tree3.txt")
rand_tree3_lm_res <- compute_logmap(tree_paths = c(paths, "data/prevotella/rand_tree3.txt"),
                                    tree_names = c(gene_names, "phylogenomic", "random"),
                                    base_lab = "random")
rand_tree3_plot_res <- plot_logmap(vectors = rand_tree3_lm_res$vectors,
                                   phylogenomic = 64,
                                   phylogenomic_name = "$\\bar{T}_p^{full}$",
                                   title = "Random base tree", 
                                   tree_names = c(gene_names, "phylogenomic", "random"),
                                   other_tree = 65,
                                   other_name = "random") 
rand_tree3_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_rand_base3.png")
RF_dists <- as.matrix(phangorn::RF.dist(c(all_trees, rand_tree3)))
diag(RF_dists) <- NA
rowMeans(RF_dists, na.rm = TRUE)
# seed of 4
set.seed(4)
rand_tree4 <- rtree(length(all_trees[[1]]$tip.label), rooted = FALSE)
rand_tree4$tip.label <- all_trees[[1]]$tip.label
write.tree(rand_tree4, "data/prevotella/rand_tree4.txt")
rand_tree4_lm_res <- compute_logmap(tree_paths = c(paths, "data/prevotella/rand_tree4.txt"),
                                    tree_names = c(gene_names, "phylogenomic", "random"),
                                    base_lab = "random")
rand_tree4_plot_res <- plot_logmap(vectors = rand_tree4_lm_res$vectors,
                                   phylogenomic = 64,
                                   phylogenomic_name = "$\\bar{T}_p^{full}$",
                                   title = "Random base tree", 
                                   tree_names = c(gene_names, "phylogenomic", "random"),
                                   other_tree = 65,
                                   other_name = "random") 
rand_tree4_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_rand_base4.png")
RF_dists <- as.matrix(phangorn::RF.dist(c(all_trees, rand_tree4)))
diag(RF_dists) <- NA
rowMeans(RF_dists, na.rm = TRUE)
