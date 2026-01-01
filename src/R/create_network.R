library(tidyverse)

load_gene_info = function() {
  read_csv('data/gene-info.csv')
}

load_correlations = function () {
  read_csv('data/corrs.csv') 
}

connectivity = function(m) {
  rowSums(m, na.rm=T)
}

adjacency = function(correlations, beta) {
  a = abs(correlations)
  0 -> a[is.na(a)]
  0 -> diag(a)

  a^beta
}

r_squared = function(m) {
  k = connectivity(m)
  binned_k = hist(k, plot=FALSE)
  p_k = binned_k$counts / sum(binned_k$counts)

  # Filter out empty bins
  nonzero = p_k > 0
  p_k[nonzero] -> p_k
  use_k = binned_k$mids[nonzero]

  cor(log10(p_k), log10(use_k))^2
}

tabulate_betas = function(correlations) {
  1:20 |>
    map_dfr(function(beta) {
      a = adjacency(correlations, beta)
      r2 = r_squared(a)
      k = connectivity(a)
      tibble(beta = beta, r2 = r2, mean_k = mean(k))
    })
}

minima = function(m) {
  k = connectivity(m)
  outer(k, k, pmin)
}

tom = function(correlations, beta) {
  a = adjacency(correlations, beta)
  (a %*% a + a) / (minima(a) + 1 - a)
}

cluster_tom = function(tom) {
  distances = as.dist(1 - tom)
  hclust(distances, method="average")
}