library(tidyverse)

# --- Data ---
load_gene_info = function() read_csv('data/gene-info.csv')

load_correlations = function () read_csv('data/corrs.csv')
load_corr_matrix = function () load_correlations() |> as.matrix()

# --- Network construction ---
connectivity = function(m) rowSums(m, na.rm=TRUE)

adjacency = function(corr_matrix, beta) {
  a = abs(corr_matrix)
  0 -> a[is.na(a)]
  0 -> diag(a)

  a^beta
}

tom = function(corr_matrix, beta) {
  a = adjacency(corr_matrix, beta)

  k = connectivity(a)
  k_min = outer(k, k, pmin)

  (a %*% a + a) / (k_min + 1 - a)
}

# --- Scale-free topology ---
k_distribution = function(a) {
  k = connectivity(a)
  binned_k = hist(k, plot = FALSE, breaks=50)
  p_k = binned_k$counts / sum(binned_k$counts)

  nonzero = p_k > 0
  list(p_k = log10(p_k[nonzero]), use_k = log10(binned_k$mids[nonzero]))
}

scale_free_r2 = function(a) {
  k_dist = k_distribution(a)
  cor(lk_dist$p_k, k_dist$use_k)^2
}

tabulate_betas = function(corr_matrix, betas = 1:20) {
  map_dfr(betas, function(beta) {
    a = adjacency(corr_matrix, beta)
    tibble(beta = beta, r2 = scale_free_r2(a), mean_k = mean(connectivity(a)))
  })
}

plot_log_log = function(corr_matrix, beta) {
  a = adjacency(corr_matrix, beta)
  k_dist = k_distribution(a)

  x = k_dist$use_k
  y = k_dist$p_k

  plot(x, y)
  abline(lm(y ~ x), col = "red")
  k_dist
}

# --- Clustering ---
cluster_genes = function(tom) {
  distances = as.dist(1 - tom)
  hclust(distances, method="average")
}