# simulation to see differences between using Frechet mean tree as base versus 
# using the tree in the set with the minimum mean squared BHV distance from all 
# other trees in the set
set.seed(1)
trees <- ape::rmtree(50, 10)
paths <- rep(NA, 50)
for (i in 1:50) {
  path <- paste0("data/frechet_sim/tree", i, ".txt")
  write.tree(trees[[i]], path)
  paths[i] <- path
}
ape::write.tree(trees, "data/frechet_sim/trees.txt")
fm_output <- system2("java", args = c("-jar", "scripts/jar/SturmMean.jar", "-o", 
                                      "data/frechet_sim/frechet_mean_tree.txt", 
                                      "data/frechet_sim/trees.txt"), stdout = T)
frechet_mean <- ape::read.tree(text = fm_output[12])
write.tree(frechet_mean, "data/frechet_sim/frechet_tree.txt")
ape::is.binary(frechet_mean)

merge_txt(c(paths, "data/frechet_sim/frechet_tree.txt"),
          "data/frechet_sim/all_trees.txt")
bhv_dists <- compute_geodesic("data/frechet_sim/all_trees.txt")
sq_dists <- bhv_dists^2
diag(sq_dists) <- NA
mean_dists <- rowMeans(sq_dists, na.rm = TRUE)

lm_res_frechet <- compute_logmap(tree_paths = c(paths, "data/frechet_sim/frechet_tree.txt"))
# we can see that the base tree was chosen to be the Frechet mean tree
print(lm_res_frechet$base_lab)
plot_res <- plot_logmap(vectors = lm_res_frechet$vectors,
                        title = "Simulated trees, Frechet base tree", phylogenomic = 51,
                        phylogenomic_name = "Frechet mean tree")
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Name), nudge_x = .2)
ggsave("figures/sim/frechet_base.png")
# now find the gene tree that has the lowest mean squared BHV distance from other trees
base <- which.min(mean_dists[1:50])
lm_res_base <- compute_logmap(tree_paths = c(paths, "data/frechet_sim/frechet_tree.txt"),
                                 base_lab = base)
# we can see that the base tree is tree 24
print(lm_res_base$base_lab)
plot_res <- plot_logmap(vectors = lm_res_base$vectors,
                        title = "Simulated trees, min mean sq distance base", phylogenomic = 51,
                        phylogenomic_name = "Frechet mean tree",
                        other_tree = base,
                        other_name = "base tree")
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Name), nudge_x = .2)
ggsave("figures/sim/min_mean_dist.png")
# what if we make the base tree number the tree with the largest mean squared BHV distance from other trees
base <- which.max(mean_dists[1:50])
lm_res_base_max <- compute_logmap(tree_paths = c(paths, "data/frechet_sim/frechet_tree.txt"),
                              base_lab = base)
# we can see that the base tree is tree 40
print(lm_res_base_max$base_lab)
plot_res <- plot_logmap(vectors = lm_res_base_max$vectors,
                        title = "Simulated trees, max mean sq distance base", phylogenomic = 51,
                        phylogenomic_name = "Frechet mean tree",
                        other_tree = base,
                        other_name = "base tree")
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Name), nudge_x = .2)
ggsave("figures/sim/max_mean_dist.png")
# what if we make the base tree the outlier from the viz with the Frechet mean tree as base
base <- 20
lm_res_base_out <- compute_logmap(tree_paths = c(paths, "data/frechet_sim/frechet_tree.txt"),
                                  base_lab = base)
plot_res <- plot_logmap(vectors = lm_res_base_out$vectors,
                        title = "Simulated trees, outlier base", phylogenomic = 51,
                        phylogenomic_name = "Frechet mean tree",
                        other_tree = base,
                        other_name = "base tree")
plot_res$plot + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = Name), nudge_x = .2)
ggsave("figures/sim/outlier_base.png")
