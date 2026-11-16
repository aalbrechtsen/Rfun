# One sample per person (k = 1) for the samples-per-person comparison (report Figure 2).
# With k = 1 there are no repeated samples: permuting individuals equals free permutation,
# individual means equal the samples, and a within-individual test is impossible.
# Appends/replaces the k = 1 rows in sim_repeated_rda_results.tsv.
# Run: OPENBLAS_NUM_THREADS=1 Rscript --vanilla sim_k1.R  (from this folder)
source("sim_repeated_rda_lib.R")
RNGkind("L'Ecuyer-CMRG"); set.seed(20171)
NREP <- 1000; NPERM <- 499; NCORE <- min(80, detectCores())
sc <- list(nind = 24, k = 1, p = 5, m = 20, rho_d = 0.8, rho_y = 0.8, R2true = 0)

res <- mclapply(seq_len(NREP), function(r) {
  s <- sim_data(sc$nind, sc$k, sc$p, sc$m, sc$rho_d, sc$rho_y, sc$R2true)
  G <- gram_euclid(s$z)
  rbind(naive = unlist(test_marg(s$d, G, perm_free(sc$nind, NPERM))),
        block = unlist(test_marg(s$d, G, perm_block(sc$nind, 1, NPERM))),
        means = unlist(test_marg(s$d, G, perm_free(sc$nind, NPERM))))
}, mc.cores = NCORE, mc.set.seed = TRUE)
arr <- simplify2array(res)
out <- data.frame(dist = "euclid", nind = sc$nind, k = 1, p = sc$p, m = sc$m, rho_d = sc$rho_d,
                  rho_y = sc$rho_y, R2true = 0, method = dimnames(arr)[[1]],
                  reject = rowMeans(arr[, "p", ] < 0.05), R2 = rowMeans(arr[, "R2", ]),
                  adjR2 = rowMeans(arr[, "adjR2", ]), adjR2_sd = apply(arr[, "adjR2", ], 1, sd))
out[] <- lapply(out, function(x) if (is.numeric(x)) round(x, 4) else x)
old <- read.delim("sim_repeated_rda_results.tsv")
old <- old[old$k != 1, ]
write.table(rbind(old, out), "sim_repeated_rda_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
