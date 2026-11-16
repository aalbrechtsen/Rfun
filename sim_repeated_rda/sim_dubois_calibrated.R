# Simulation matched to Dubois et al. (2017, Microbiome 5:151), checked against the full text:
#  - 128 toilet-paper samples from 24 people: 15 Nunavut, 9 Montreal (mean 5.3 per person).
#    Per-person counts are not reported; Table S1 lists 61 Nunavut samples, so we use
#    4 per Nunavut participant (60) and 7-8 per Montreal participant (68).
#  - 16 food categories as RDA predictors (Fig. 3; Table S3 lists 15 + cooked game meat),
#    mean frequency per day, centred and scaled.
#  - Response: square-root-transformed unweighted UniFrac distances; RDA in vegan,
#    adjusted R2 (RsquareAdj); permutation tests with free permutation (no strata stated).
#  - Analyses: diet (adj R2 0.11), diet | geography (partial RDA, 0.064),
#    diet within Nunavut only (0.169).
#  - Individual identity explains 45% of unweighted UniFrac variation (Fig. 2b),
#    geography 4% (Fig. 2a).
# Diet has NO effect on the microbiome unless R2true > 0.
#
# Methods compared:
#   naive : all samples, rows permuted freely (as in the paper)
#   block : same RDA, whole individuals permuted (within geography and number of samples)
#   means : individual centroids on mean diet (n = 24); not possible for Nunavut only (15 people, 16 predictors)
#
# Run: OPENBLAS_NUM_THREADS=1 Rscript --vanilla sim_dubois_calibrated.R  (from this folder)

source("sim_repeated_rda_lib.R")
set.seed(2017151)
RNGkind("L'Ecuyer-CMRG")
NREP <- 1000; NPERM <- 499; NCORE <- min(80, detectCores())

# permute whole individuals among individuals in the same stratum (same size, same region)
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

sim_dubois <- function(sc) {
  kvec <- c(rep(sc$kN, sc$nN), sc$kM); nind <- length(kvec); p <- sc$p; m <- sc$m
  id <- rep(seq_len(nind), kvec); n <- length(id)
  g_ind <- rep(c(1, 0), c(sc$nN, length(sc$kM))); g <- g_ind[id]
  gc <- g - mean(g)
  # diet: geography shift + individual habit + time-to-time variation
  u <- rnorm(p); u <- u / sqrt(mean(u^2))
  a <- matrix(rnorm(nind * p), nind)[id, , drop = FALSE]
  d <- outer(gc, u) * sc$geo_diet +
       sqrt(sc$rho_d) * a + sqrt(1 - sc$rho_d) * matrix(rnorm(n * p), n)
  d <- scale(d)
  # microbiome latent: geography + individual + time-to-time
  v <- rnorm(m); v <- v / sqrt(sum(v^2))
  b <- matrix(rnorm(nind * m), nind)[id, , drop = FALSE]
  z <- sqrt(sc$rho_y) * b + sqrt(1 - sc$rho_y) * matrix(rnorm(n * m), n)
  z <- z + outer(gc, v) * sqrt(sc$geo_mb * m / (1 - sc$geo_mb) / var(gc))
  if (sc$R2true > 0) {
    Gam <- matrix(rnorm(p * m), p)
    Gam <- Gam * sqrt(sc$R2true / (1 - sc$R2true) * sum(apply(z, 2, var)) / sum(Gam^2))
    z <- z + d %*% Gam
  }
  list(id = id, g = g, g_ind = g_ind, kvec = kvec, d = d, z = z)
}

centroids <- function(G, id) {
  A <- model.matrix(~ 0 + factor(id)); A <- sweep(A, 2, colSums(A), "/")
  Gm <- t(A) %*% G %*% A
  list(A = A, Gm = Gm - outer(rowMeans(Gm), colMeans(Gm), "+") + mean(Gm))
}

one_rep_dub <- function(sc, tree) {
  s <- sim_dubois(sc); n <- length(s$id); nind <- length(s$kvec); B <- NPERM
  G <- if (sc$dist == "euclid") gram_euclid(s$z) else
    gower(sqrt(z_to_uf(s$z, tree, sc$load, sc$mu)))          # sqrt unweighted UniFrac, as in the paper
  cm <- centroids(G, s$id); dm <- t(cm$A) %*% s$d
  Zid <- model.matrix(~ factor(s$id))[, -1, drop = FALSE]
  geo <- matrix(s$g); geo_m <- matrix(s$g_ind)
  strata <- paste(s$g_ind, s$kvec)
  pb <- perm_units(s$id, strata, B); pf <- perm_free(n, B)
  # Nunavut only
  nu <- which(s$g == 1); Gn <- G[nu, nu]; Gn <- Gn - outer(rowMeans(Gn), colMeans(Gn), "+") + mean(Gn)
  dn <- scale(s$d[nu, ]); idn <- s$id[nu]
  r <- list(
    diet_naive    = test_marg(s$d, G, pf),
    diet_block    = test_marg(s$d, G, pb),
    diet_means    = test_marg(dm, cm$Gm, perm_free(nind, B)),
    dietGeo_naive = test_partial(s$d, geo, G, pf),
    dietGeo_block = test_partial(s$d, geo, G, pb),
    dietGeo_means = test_partial(dm, geo_m, cm$Gm, perm_units(seq_len(nind), strata, B)),
    nunavut_naive = test_marg(dn, Gn, perm_free(length(nu), B)),
    nunavut_block = test_marg(dn, Gn, perm_units(idn, rep(1, sc$nN), B)))
  out <- t(sapply(r, unlist))
  cbind(out, ind_R2 = test_marg(Zid, G, perm_free(n, 1))$R2, geo_R2 = test_marg(geo, G, perm_free(n, 1))$R2)
}

run_dub <- function(sc, nrep) {
  tree <- NULL
  if (sc$dist == "unifrac") {
    tree <- make_tree(sc$ntax)
    sc$load <- matrix(rnorm(sc$m * sc$ntax, sd = 0.8), sc$m); sc$mu <- rnorm(sc$ntax, sd = 2)
  }
  arr <- simplify2array(mclapply(seq_len(nrep), function(i) one_rep_dub(sc, tree),
                                 mc.cores = NCORE, mc.set.seed = TRUE))
  data.frame(dist = sc$dist, rho_d = sc$rho_d, rho_y = sc$rho_y, R2true = sc$R2true,
             method = dimnames(arr)[[1]],
             reject = rowMeans(arr[, "p", ] < 0.05),
             R2 = rowMeans(arr[, "R2", ]), adjR2 = rowMeans(arr[, "adjR2", ]),
             adjR2_q05 = apply(arr[, "adjR2", ], 1, quantile, 0.05),
             adjR2_q95 = apply(arr[, "adjR2", ], 1, quantile, 0.95),
             ind_R2 = mean(arr[1, "ind_R2", ]), geo_R2 = mean(arr[1, "geo_R2", ]), row.names = NULL)
}

base <- list(nN = 15, kN = 4, kM = c(8, 8, 8, 8, 8, 7, 7, 7, 7), p = 16, m = 20,
             rho_d = 0.6, rho_y = 0.55, geo_diet = 1, geo_mb = 0.04, R2true = 0, ntax = 300,
             dist = "euclid")
# microbiome settings per distance, calibrated so that individual identity explains ~45% and
# geography ~4% of variation, as for unweighted UniFrac in the paper (Fig. 2)
CAL <- list(euclid = list(rho_y = 0.35, geo_mb = 0.04), unifrac = list(rho_y = 0.8, geo_mb = 0.09))
mod <- function(...) { sc <- modifyList(base, list(...)); modifyList(sc, CAL[[sc$dist]]) }
scen <- list()
for (dist in c("euclid", "unifrac")) for (rd in c(0, 0.3, 0.6, 0.8))
  scen[[length(scen) + 1]] <- mod(dist = dist, rho_d = rd)
for (dist in c("euclid", "unifrac"))
  scen[[length(scen) + 1]] <- mod(dist = dist, R2true = 0.05)

res <- do.call(rbind, lapply(scen, function(sc) run_dub(sc, if (sc$dist == "unifrac") 500 else NREP)))
res[] <- lapply(res, function(x) if (is.numeric(x)) round(x, 4) else x)
write.table(res, "sim_dubois_calibrated_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
