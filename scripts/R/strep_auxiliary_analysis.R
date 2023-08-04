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
lm_res_strep <- compute_logmap(tree_paths = paths,
                               tree_names = c(gene_names, "phylogenomic"))
# we can see that the base tree was chosen to be the phylogenomic tree
print(lm_res_strep$base_lab)
plot_res <- plot_logmap(vectors = lm_res_strep$vectors, phylogenomic = length(paths),
                        title = "Phylogenomic base tree", tree_names = c(gene_names, "phylogenomic"),
                        phylogenomic_name = "$\\bar{T}_p^{full}$")
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz_phylo_base.jpeg")

# compare this to the visualization with the log map with different base trees 
# rank trees by mean squared distance from other trees 
merge_txt(paths, "data/strep/all_trees.txt")
all_trees <- read.tree("data/strep/all_trees.txt")
bhv_dists <- compute_geodesic("data/strep/all_trees.txt")
sq_dists <- bhv_dists^2
diag(sq_dists) <- NA
mean_dists <- rowMeans(sq_dists, na.rm = TRUE)
rank_df <- data.frame(rank = rank(mean_dists), mean_dist = mean_dists, 
                      name = c(gene_names, "phylogenomic")) %>%
  arrange(rank)
head(rank_df)
# we can see that the phylogenomic tree has the lowest mean squared distance from 
# the other trees, followed by the trees for genes MnmE_helical, NFACT_N, 
# Ribosomal_L9_C, NAPRTase, Ribosomal_S30AE

# log map visualization with MnmE_helical as a base tree 
MnmE_helical_lm_res <- compute_logmap(tree_paths = paths,
                              tree_names = c(gene_names, "phylogenomic"),
                              base_lab = "MnmE_helical")
MnmE_helical_plot_res <- plot_logmap(vectors = MnmE_helical_lm_res$vectors, phylogenomic = 64,
                             title = "MnmE_helical base tree", tree_names = c(gene_names, "phylogenomic"),
                             phylogenomic_name = "$\\bar{T}_p^{full}$",
                             other_tree = which(gene_names == "MnmE_helical"),
                             other_name = "MnmE_helical") 
MnmE_helical_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz_MnmE_helical_base.png")

# log map visualization with NFACT_N as a base tree 
NFACT_N_lm_res <- compute_logmap(tree_paths = paths,
                                      tree_names = c(gene_names, "phylogenomic"),
                                      base_lab = "NFACT_N")
NFACT_N_plot_res <- plot_logmap(vectors = NFACT_N_lm_res$vectors, phylogenomic = 64,
                                     title = "NFACT_N base tree", tree_names = c(gene_names, "phylogenomic"),
                                     phylogenomic_name = "$\\bar{T}_p^{full}$",
                                     other_tree = which(gene_names == "NFACT_N"),
                                     other_name = "NFACT_N") 
NFACT_N_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz_NFACT_N_base.png")

# log map visualization with Ribosomal_L9_C as a base tree 
Ribosomal_L9_C_lm_res <- compute_logmap(tree_paths = paths,
                                 tree_names = c(gene_names, "phylogenomic"),
                                 base_lab = "Ribosomal_L9_C")
# Ribosomal_L9_C is not resolved so we can't use it as a base tree 

# log map visualization with NAPRTase as a base tree 
NAPRTase_lm_res <- compute_logmap(tree_paths = paths,
                                 tree_names = c(gene_names, "phylogenomic"),
                                 base_lab = "NAPRTase")
NAPRTase_plot_res <- plot_logmap(vectors = NAPRTase_lm_res$vectors, phylogenomic = 64,
                                title = "NAPRTase base tree", tree_names = c(gene_names, "phylogenomic"),
                                phylogenomic_name = "$\\bar{T}_p^{full}$",
                                other_tree = which(gene_names == "NAPRTase"),
                                other_name = "NAPRTase") 
NAPRTase_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz_NAPRTase_base.png")

# log map visualization with Ribosomal_S30AE as the base tree 
Ribosomal_S30AE_lm_res <- compute_logmap(tree_paths = paths,
                                  tree_names = c(gene_names, "phylogenomic"),
                                  base_lab = "Ribosomal_S30AE")
# Ribonsomal_S30AE is not resolved so we can't use it as a base tree

# look at effect of using outlier trees as base trees
# log map visualization with DUF3270 as a base tree 
DUF3270_lm_res <- compute_logmap(tree_paths = paths,
                                  tree_names = c(gene_names, "phylogenomic"),
                                  base_lab = "DUF3270")
# DUF3270 is not resolved so we can't use it as a base tree 

# log map visualization with EcsB as a base tree 
EcsB_lm_res <- compute_logmap(tree_paths = paths,
                                  tree_names = c(gene_names, "phylogenomic"),
                                  base_lab = "EcsB")
EcsB_plot_res <- plot_logmap(vectors = EcsB_lm_res$vectors, phylogenomic = 64,
                                 title = "EcsB base tree", tree_names = c(gene_names, "phylogenomic"),
                                 phylogenomic_name = "$\\bar{T}_p^{full}$",
                                 other_tree = which(gene_names == "EcsB"),
                                 other_name = "EcsB") 
EcsB_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz_EcsB_base.png")

# log map visualization with DUF1934 as a base tree 
DUF1934_lm_res <- compute_logmap(tree_paths = paths,
                              tree_names = c(gene_names, "phylogenomic"),
                              base_lab = "DUF1934")
DUF1934_plot_res <- plot_logmap(vectors = DUF1934_lm_res$vectors, phylogenomic = 64,
                             title = "DUF1934 base tree", tree_names = c(gene_names, "phylogenomic"),
                             phylogenomic_name = "$\\bar{T}_p^{full}$",
                             other_tree = which(gene_names == "DUF1934"),
                             other_name = "DUF1934") 
DUF1934_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz_DUF1934_base.png")

# also investigate Frechet mean tree
fm_output <- system2("java", args = c("-jar", "scripts/jar/SturmMean.jar", "-o", 
                                      "data/strep/frechet_mean_tree.txt", 
                                      "data/strep/gene_trees.txt"), stdout = T)
frechet_mean <- ape::read.tree(text = fm_output[12])
ape::is.binary(frechet_mean)
# distance between phylogenomic tree and frechet mean tree
merge_txt(c(paths, "data/strep/frechet_mean_tree.txt"), 
          "data/strep/phylo_and_frechet.txt")
phylo_frechet_dist <- compute_geodesic("data/strep/phylo_and_frechet.txt")
dist_vec <- phylo_frechet_dist[lower.tri(phylo_frechet_dist)]
phy_bhv <- phylo_frechet_dist[197, 198]
mean(dist_vec < phy_bhv)

# run MDS on the same data
MDS_res <- compute_MDS(dist_matrix = bhv_dists,
                       tree_names = c(gene_names, "phylogenomic"))
MDS_plot <- plot_MDS(df = MDS_res$df, phylogenomic = 64,
                     title = "MDS of Streptococcus trees", tree_names = c(gene_names, "phylogenomic"),
                     phylogenomic_name = "$\\bar{T}_p^{full}$") 
MDS_plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/MDS_plot.png")

# run Fisher's discriminant analysis on log map vectors between ribosomal and 
# non-ribomsomal genes 
ribs <- grep("Ribosom", gene_names)
ribs <- ifelse(1:length(paths) %in% ribs, "ribosomal", "other")
lm_df <- as.data.frame(lm_res_strep$vectors)
lm_df$rib <- ribs
# run LDA on 75% of trees 
set.seed(1)
n <- nrow(lm_df)
training_ind <- sample(x = 1:n, size = round(n*.75))
testing_ind <- (1:n)[-training_ind]
lda_res <- MASS::lda(rib ~ ., data = lm_df[training_ind, ])
# predict label on remaining 25% of trees
lda_pred <- predict(lda_res, newdata = lm_df[testing_ind, ])
preds <- lda_pred$class  
correct <- mean(preds == lm_df$rib[testing_ind])  
# 63% of the test set is classified correctly 
rank_df <- data.frame(coord = 1:nrow(lda_res$scaling), 
                      scale = lda_res$scaling[,1],
                      rank = rank(abs(lda_res$scaling[,1]))) %>%
  arrange(desc(rank))
head(rank_df)
# we can see that the log map vectors with the largest scalings used to 
# separate the two classes are vectors 78, 17, and 86. We could map these back
# to tree space and see which internal edges of the base tree map to edges 
# 78, 17, and 86
phylo <- read.tree("data/strep/phylogenomic_trees/concat_tree.txt")
phylo_mid <- phytools::midpoint.root(phylo)
phylo_mid$tip.label <- 1:length(phylo_mid$tip.label)
ggtree(phylo_mid) + geom_tiplab(size = 2) + 
  xlim(c(0, .4)) + ggtitle("Phylogenomic tree") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/phylogenomic_tree.jpeg")
# color those internal edges
order_df <- data.frame(phylo_rank = rank(phylo$edge.length),
                       vector_rank = rank(lm_res_strep$vectors[197,]))
tree_ind <- order_df$phylo_rank[which(order_df$vector_rank == 78)]
ggtree(phylo_mid) + geom_tiplab(size = 2) + 
  xlim(c(0, .4)) + ggtitle("Phylogenomic tree") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_highlight(node = phylo$edge[tree_ind, 2], fill = "red")
ggsave("figures/strep/phylo_with_lab.jpeg")

# use a randomly generated tree as the base tree
# seed of 1
set.seed(1)
rand_tree1 <- rtree(length(all_trees[[1]]$tip.label), rooted = FALSE)
rand_tree1$tip.label <- all_trees[[1]]$tip.label
write.tree(rand_tree1, "data/strep/rand_tree1.txt")
rand_tree1_lm_res <- compute_logmap(tree_paths = c(paths, "data/strep/rand_tree1.txt"),
                                    tree_names = c(gene_names, "phylogenomic", "random"),
                                    base_lab = "random")
rand_tree1_plot_res <- plot_logmap(vectors = rand_tree1_lm_res$vectors,
                                   phylogenomic = 197,
                                   phylogenomic_name = "$\\bar{T}_p^{full}$",
                                   title = "Random base tree", 
                                   tree_names = c(gene_names, "phylogenomic", "random"),
                                   other_tree = 198,
                                   other_name = "random") 
rand_tree1_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/strep/lm_viz_rand_base1.png")
RF_dists <- as.matrix(phangorn::RF.dist(c(all_trees, rand_tree1)))
diag(RF_dists) <- NA
rowMeans(RF_dists, na.rm = TRUE)
# seed of 2
set.seed(2)
rand_tree2 <- rtree(length(all_trees[[1]]$tip.label), rooted = FALSE)
rand_tree2$tip.label <- all_trees[[1]]$tip.label
write.tree(rand_tree2, "data/strep/rand_tree2.txt")
rand_tree2_lm_res <- compute_logmap(tree_paths = c(paths, "data/strep/rand_tree2.txt"),
                                    tree_names = c(gene_names, "phylogenomic", "random"),
                                    base_lab = "random")
rand_tree2_plot_res <- plot_logmap(vectors = rand_tree2_lm_res$vectors,
                                   phylogenomic = 197,
                                   phylogenomic_name = "$\\bar{T}_p^{full}$",
                                   title = "Random base tree", 
                                   tree_names = c(gene_names, "phylogenomic", "random"),
                                   other_tree = 198,
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
write.tree(rand_tree3, "data/strep/rand_tree3.txt")
rand_tree3_lm_res <- compute_logmap(tree_paths = c(paths, "data/strep/rand_tree3.txt"),
                                    tree_names = c(gene_names, "phylogenomic", "random"),
                                    base_lab = "random")
rand_tree3_plot_res <- plot_logmap(vectors = rand_tree3_lm_res$vectors,
                                   phylogenomic = 197,
                                   phylogenomic_name = "$\\bar{T}_p^{full}$",
                                   title = "Random base tree", 
                                   tree_names = c(gene_names, "phylogenomic", "random"),
                                   other_tree = 198,
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
write.tree(rand_tree4, "data/strep/rand_tree4.txt")
rand_tree4_lm_res <- compute_logmap(tree_paths = c(paths, "data/strep/rand_tree4.txt"),
                                    tree_names = c(gene_names, "phylogenomic", "random"),
                                    base_lab = "random")
rand_tree4_plot_res <- plot_logmap(vectors = rand_tree4_lm_res$vectors,
                                   phylogenomic = 197,
                                   phylogenomic_name = "$\\bar{T}_p^{full}$",
                                   title = "Random base tree", 
                                   tree_names = c(gene_names, "phylogenomic", "random"),
                                   other_tree = 198,
                                   other_name = "random") 
rand_tree4_plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/prevotella/lm_viz_rand_base4.png")
RF_dists <- as.matrix(phangorn::RF.dist(c(all_trees, rand_tree4)))
diag(RF_dists) <- NA
rowMeans(RF_dists, na.rm = TRUE)
