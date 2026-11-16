# Generic UniFrac scenarios (report Figure 3) with square-root-transformed unweighted UniFrac,
# the distance used by Dubois et al. (2017). Same settings as the UniFrac rows in
# sim_repeated_rda.R: 24 individuals x 4 samples, 5 diet variables, no diet effect,
# cor(diet | ID) = cor(microbiome | ID) = 0, 0.5, 0.8; 500 replicates.
# Appends/replaces rows with dist = "unifrac_sqrt" in sim_repeated_rda_results.tsv.
# Run: OPENBLAS_NUM_THREADS=1 Rscript --vanilla sim_unifrac_sqrt.R  (from this folder)
source("sim_repeated_rda_lib.R")
RNGkind("L'Ecuyer-CMRG"); set.seed(20172)
NREP <- 500; NPERM <- 499; NCORE <- min(80, detectCores())

one_rep_sqrt <- function(sc, tree) {
  s <- sim_data(sc$nind, sc$k, sc$p, sc$m, sc$rho_d, sc$rho_y, sc$R2true)
  n <- sc$nind * sc$k
  G <- gower(sqrt(z_to_uf(s$z, tree, sc$load, sc$mu)))
  A <- model.matrix(~ 0 + factor(s$id)) / sc$k
  Gm <- t(A) %*% G %*% A; Gm <- Gm - outer(rowMeans(Gm), colMeans(Gm), "+") + mean(Gm)
  rbind(naive = unlist(test_marg(s$d, G, perm_free(n, NPERM))),
        block = unlist(test_marg(s$d, G, perm_block(sc$nind, sc$k, NPERM))),
        means = unlist(test_marg(t(A) %*% s$d, Gm, perm_free(sc$nind, NPERM))))
}

out <- do.call(rbind, lapply(c(0, 0.5, 0.8), function(rr) {
  sc <- list(nind = 24, k = 4, p = 5, m = 20, rho_d = rr, rho_y = rr, R2true = 0, ntax = 300)
  tree <- make_tree(sc$ntax)
  sc$load <- matrix(rnorm(sc$m * sc$ntax, sd = 0.8), sc$m); sc$mu <- rnorm(sc$ntax, sd = 2)
  arr <- simplify2array(mclapply(seq_len(NREP), function(i) one_rep_sqrt(sc, tree),
                                 mc.cores = NCORE, mc.set.seed = TRUE))
  data.frame(dist = "unifrac_sqrt", nind = 24, k = 4, p = 5, m = 20, rho_d = rr, rho_y = rr, R2true = 0,
             method = dimnames(arr)[[1]], reject = rowMeans(arr[, "p", ] < 0.05),
             R2 = rowMeans(arr[, "R2", ]), adjR2 = rowMeans(arr[, "adjR2", ]),
             adjR2_sd = apply(arr[, "adjR2", ], 1, sd))
}))
out[] <- lapply(out, function(x) if (is.numeric(x)) round(x, 4) else x)
old <- read.delim("sim_repeated_rda_results.tsv")
old <- old[old$dist != "unifrac_sqrt", ]
write.table(rbind(old, out), "sim_repeated_rda_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
