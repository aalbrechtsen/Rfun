# Report Figure 3: generic UniFrac scenarios with the paper's sample sizes and number of diet variables.
#  - 24 individuals, 128 samples: 15 with 4 samples, 5 with 8, 4 with 7 (as assumed for Dubois et al.)
#  - 16 diet variables (Dubois et al. Fig. 3b)
#  - no geography, no diet effect; cor(diet | ID) = cor(microbiome | ID) = 0, 0.5, 0.8
#  - unweighted UniFrac and square-root unweighted UniFrac (the paper's distance)
# Methods: naive RDA (free permutation, as in the paper), RDA with ID permutation (whole individuals
# swapped among individuals with the same number of samples), individual means (n = 24) computed from
# 4 randomly chosen samples per person, so every mean has the same precision and individuals are exchangeable.
# Appends/replaces rows with dist "unifrac_paper" / "unifrac_sqrt_paper" in sim_repeated_rda_results.tsv
# (k is stored as 5.33, the mean number of samples per person).
# Run: OPENBLAS_NUM_THREADS=1 Rscript --vanilla sim_unifrac_paper.R  (from this folder)
source("sim_repeated_rda_lib.R")
RNGkind("L'Ecuyer-CMRG"); set.seed(20173)
NREP <- 500; NPERM <- 499; NCORE <- min(80, detectCores())
KVEC <- c(rep(4, 15), rep(8, 5), rep(7, 4)); P <- 16; M <- 20

perm_units <- function(id, strata_ind, B) {
  rows <- split(seq_along(id), id)
  replicate(B, {
    o <- seq_along(rows)
    for (s in unique(strata_ind)) { w <- which(strata_ind == s); o[w] <- w[sample.int(length(w))] }
    idx <- integer(length(id))
    for (i in seq_along(rows)) idx[rows[[i]]] <- rows[[o[i]]]
    idx
  })
}

one_rep <- function(rr, sqrt_d, tree, load, mu) {
  nind <- length(KVEC); id <- rep(seq_len(nind), KVEC); n <- length(id)
  d <- sqrt(rr) * matrix(rnorm(nind * P), nind)[id, ] + sqrt(1 - rr) * matrix(rnorm(n * P), n)
  z <- sqrt(rr) * matrix(rnorm(nind * M), nind)[id, ] + sqrt(1 - rr) * matrix(rnorm(n * M), n)
  D <- z_to_uf(z, tree, load, mu); if (sqrt_d) D <- sqrt(D)
  G <- gower(D)
  # individual means from 4 random samples per person
  sub <- unlist(lapply(split(seq_len(n), id), function(r) r[sample.int(length(r), 4)]))
  A <- model.matrix(~ 0 + factor(id[sub])) / 4
  Gm <- t(A) %*% G[sub, sub] %*% A; Gm <- Gm - outer(rowMeans(Gm), colMeans(Gm), "+") + mean(Gm)
  rbind(naive = unlist(test_marg(d, G, perm_free(n, NPERM))),
        block = unlist(test_marg(d, G, perm_units(id, KVEC, NPERM))),
        means = unlist(test_marg(t(A) %*% d[sub, ], Gm, perm_free(nind, NPERM))))
}

out <- do.call(rbind, lapply(c(FALSE, TRUE), function(sqrt_d) do.call(rbind, lapply(c(0, 0.5, 0.8), function(rr) {
  tree <- make_tree(300); load <- matrix(rnorm(M * 300, sd = 0.8), M); mu <- rnorm(300, sd = 2)
  arr <- simplify2array(mclapply(seq_len(NREP), function(i) one_rep(rr, sqrt_d, tree, load, mu),
                                 mc.cores = NCORE, mc.set.seed = TRUE))
  data.frame(dist = if (sqrt_d) "unifrac_sqrt_paper" else "unifrac_paper", nind = 24,
             k = round(mean(KVEC), 2), p = P, m = M, rho_d = rr, rho_y = rr, R2true = 0,
             method = dimnames(arr)[[1]], reject = rowMeans(arr[, "p", ] < 0.05),
             R2 = rowMeans(arr[, "R2", ]), adjR2 = rowMeans(arr[, "adjR2", ]),
             adjR2_sd = apply(arr[, "adjR2", ], 1, sd))
}))))
out[] <- lapply(out, function(x) if (is.numeric(x)) round(x, 4) else x)
old <- read.delim("sim_repeated_rda_results.tsv")
old <- old[!old$dist %in% c("unifrac_paper", "unifrac_sqrt_paper"), ]
write.table(rbind(old, out), "sim_repeated_rda_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out[, c("dist", "rho_d", "method", "reject", "adjR2", "adjR2_sd")], row.names = FALSE)
