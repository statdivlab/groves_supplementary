# This code produces the figures in the section 5 of [..... update here]
# Note that this script will not run because the data is not included in this repository. 
# This is because unlike the two data analyses, this is not publically available data. 

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
# 09/22/22 run uses:
# groves_0.0.0.9000, ggtree_3.4.2, tidyverse_1.3.2, ape_5.6-2
# data was generated with GTDB R202

# read in gene names 
gene_names <- stringr::str_sub(list.files("../TreeViz/data/prevotella/analysis_data/gene_trees/"),
                               1, -5)

# start by computing the log map visualization for all the gene trees and the phylogenomic
# tree created from concatenating all 65 genes

# we will do this twice for two versions of the phylogenomic tree created with two different
# runs of IQTREE. These different trees are referred to as 'concat_tree_v1' and 'concat_tree_v2'.

# version 1 
# make a vector of paths to .txt files of all trees to go into the plot 
paths_v1 <- c(paste0("../TreeViz/data/prevotella/analysis_data/gene_trees/", 
                  list.files("../TreeViz/data/prevotella/analysis_data/gene_trees/")),
           "../TreeViz/data/prevotella/analysis_data/phylogenomic_trees/concat_tree_v1.txt")
# compute log map coordinates for all trees 
lm_res_v1 <- compute_logmap(tree_paths = paths_v1,
                         tree_names = c(gene_names, "phylogenomic"))
# we can see that the base tree is the phylogenomic tree
print(lm_res_v1$base_lab)
plot_res <- plot_logmap(vectors = lm_res_v1$vectors, phylogenomic = length(paths_v1),
            title = "PCA of Euclidean representation of Prevotella trees", tree_names = c(gene_names, "phylogenomic"),
            phylogenomic_name = "$\\bar{T}_p$, run 2",
            trees_to_label = c("BacA"))
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/other_prevotella/lm_viz_run2.jpeg", width = 7, height = 3.5, dpi = 300)

# version 2 
# make a vector of paths to .txt files of all trees to go into the plot 
paths_v2 <- c(paste0("../TreeViz/data/prevotella/analysis_data/gene_trees/", 
                     list.files("../TreeViz/data/prevotella/analysis_data/gene_trees/")),
              "../TreeViz/data/prevotella/analysis_data/phylogenomic_trees/concat_tree_v2.txt")
# compute log map coordinates for all trees with the phylogenomic tree as the base tree
lm_res_v2 <- compute_logmap(tree_paths = paths_v2,
                         tree_names = c(gene_names, "phylogenomic"),
                         base_lab = "phylogenomic")
print(lm_res_v2$base_lab)
plot_res <- plot_logmap(vectors = lm_res_v2$vectors, phylogenomic = length(paths_v2),
            title = "PCA of Euclidean representation of Prevotella trees", tree_names = c(gene_names, "phylogenomic"),
            phylogenomic_name = "$\\bar{T}_p$, run 1",
            trees_to_label = "BacA")
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/other_prevotella/lm_viz_run1.jpeg", width = 7, height = 3.5, dpi = 300)
# compute log map coordinates for all trees with minimum mean BHV tree as base tree 
lm_res_v2.2 <- compute_logmap(tree_paths = paths_v2,
                              tree_names = c(gene_names, "phylogenomic"))
# we can see the base tree is HSP90
print(lm_res_v2.2$base_lab)
plot_res <- plot_logmap(vectors = lm_res_v2.2$vectors, phylogenomic = length(paths_v2),
            title = "", tree_names = c(gene_names, "phylogenomic"),
            phylogenomic_name = "$\\bar{T}_p^{full}$, run 1",
            trees_to_label = "BacA")
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/other_prevotella/lm_viz_v2.2.png")

# look at gene tree support 
merge_txt(paths_v1[1:65], "../TreeViz/data/prevotella/analysis_data/gene_trees.txt")
gene_trees <- read.tree("../TreeViz/data/prevotella/analysis_data/gene_trees.txt")
concat_v1 <- ape::read.tree("../TreeViz/data/prevotella/analysis_data/phylogenomic_trees/concat_tree_v1.txt")
concat_v2 <- ape::read.tree("../TreeViz/data/prevotella/analysis_data/phylogenomic_trees/concat_tree_v2.txt")
# v1
support_v1 <- check_gene_support(main_tree = concat_v1,
                              trees = gene_trees,
                              rooted = FALSE)
median(support_v1)
plot_support(concat_v1, support_v1, color_branch = TRUE,
             title = "Phylogenomic Tree Gene Support", xlim_max = 1.1)
ggsave("figures/other_prevotella/gene_tree_support_v1.png")
# v2
support_v2 <- check_gene_support(main_tree = concat_v2,
                                 trees = gene_trees,
                                 rooted = FALSE)
median(support_v2)
plot_support(concat_v2, support_v2, color_branch = TRUE,
             title = "Phylogenomic Tree Gene Support", xlim_max = 1.1)
ggsave("figures/other_prevotella/gene_tree_support_v2.png")
