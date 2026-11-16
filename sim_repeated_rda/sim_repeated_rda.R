# Simulation: distance-based RDA with repeated samples per individual
# (motivated by Dubois et al. 2017: 24 individuals, repeated samples, unweighted UniFrac).
#
# Diet (p variables) and microbiome are both correlated within individuals
# (intraclass correlations rho_d and rho_y). We compare
#   naive    : dbRDA on all samples, rows permuted freely (what ignores individual identity)
#   block    : same statistic, whole individuals permuted (restricted permutation)
#   means    : dbRDA on per-individual centroids / mean diet (n = number of individuals)
#   within   : dbRDA conditioned on individual (Condition(id)), permutations within individual
# and report rejection rate (alpha = 0.05), R2 and adjusted R2 against the true population R2.
#
# Run: OPENBLAS_NUM_THREADS=1 Rscript --vanilla sim_repeated_rda.R   (from this folder; uses many cores)

suppressPackageStartupMessages({ library(parallel); library(ape) })
set.seed(20171)
NREP  <- 1000
NPERM <- 499
NCORE <- min(80, detectCores())

## ---------- core dbRDA machinery (McArdle & Anderson 2001) ----------
gower <- function(D) {                         # D: distance matrix -> Gower-centred G
  A <- -0.5 * as.matrix(D)^2
  n <- nrow(A); J <- diag(n) - 1/n
  J %*% A %*% J
}
gram_euclid <- function(Y) { Yc <- scale(Y, scale = FALSE); tcrossprod(Yc) }
proj_tr <- function(X, G) {                    # tr(H_X G) for column-centred X
  Q <- qr.Q(qr(X)); sum(Q * (G %*% Q))
}
cent <- function(X) scale(X, scale = FALSE)
resid_on <- function(X, Z) qr.resid(qr(Z), X)

# marginal test of X on G; perms is a matrix of row permutations (one per column)
test_marg <- function(X, G, perms) {
  n <- nrow(G); p <- ncol(X); tot <- sum(diag(G))
  X <- cent(X)
  ss <- proj_tr(X, G)
  Fobs <- (ss / p) / ((tot - ss) / (n - p - 1))
  Fp <- apply(perms, 2, function(o) { s <- proj_tr(X[o, , drop = FALSE], G)
    (s / p) / ((tot - s) / (n - p - 1)) })
  R2 <- ss / tot
  list(R2 = R2, adjR2 = 1 - (1 - R2) * (n - 1) / (n - p - 1),
       p = (sum(Fp >= Fobs - 1e-12) + 1) / (ncol(perms) + 1))
}

# partial test X | Z (Z = individual dummies); permute X-residuals within individuals
test_partial <- function(X, Z, G, perms) {
  n <- nrow(G); p <- ncol(X); q <- ncol(Z); tot <- sum(diag(G))
  QZ <- qr.Q(qr(cbind(1, Z))); M <- diag(n) - tcrossprod(QZ)
  Gr <- M %*% G %*% M; totr <- sum(diag(Gr))
  Xr <- M %*% X
  ss <- proj_tr(Xr, Gr)
  df_res <- n - p - q
  Fobs <- (ss / p) / ((totr - ss) / df_res)
  Fp <- apply(perms, 2, function(o) { s <- proj_tr(Xr[o, , drop = FALSE], Gr)
    (s / p) / ((totr - s) / df_res) })
  R2 <- ss / tot                                  # fraction of total variation
  # adjusted, as vegan's varpart unique fraction: adj(X+Z) - adj(Z)
  R2z <- (tot - totr) / tot
  adj <- function(r, k) 1 - (1 - r) * (n - 1) / (n - k - 1)
  list(R2 = R2, adjR2 = adj(R2 + R2z, p + q) - adj(R2z, q),
       p = (sum(Fp >= Fobs - 1e-12) + 1) / (ncol(perms) + 1))
}

perm_free   <- function(n, B) replicate(B, sample.int(n))
perm_block  <- function(nind, k, B) replicate(B, {   # swap whole individuals, keep order
  o <- sample.int(nind); as.vector(t(outer((o - 1) * k, 1:k, "+"))) })
perm_within <- function(nind, k, B) replicate(B,
  as.vector(sapply(0:(nind - 1), function(i) i * k + sample.int(k))))

## ---------- data generation ----------
# diet: d_ij = sqrt(rho_d) a_i + sqrt(1-rho_d) e_ij  (each variable variance 1)
# latent microbiome (m dims): z_ij = d_ij Gamma + sqrt(rho_y) b_i + sqrt(1-rho_y) w_ij
# Gamma scaled so that the population R2 of diet on z equals R2true.
sim_data <- function(nind, k, p, m, rho_d, rho_y, R2true) {
  id <- rep(seq_len(nind), each = k); n <- nind * k
  a <- matrix(rnorm(nind * p), nind)[id, , drop = FALSE]
  d <- sqrt(rho_d) * a + sqrt(1 - rho_d) * matrix(rnorm(n * p), n)
  b <- matrix(rnorm(nind * m), nind)[id, , drop = FALSE]
  noise <- sqrt(rho_y) * b + sqrt(1 - rho_y) * matrix(rnorm(n * m), n)
  z <- noise
  if (R2true > 0) {
    Gam <- matrix(rnorm(p * m), p)
    Gam <- Gam * sqrt(R2true / (1 - R2true) * m / sum(Gam^2))
    z <- z + d %*% Gam
  }
  list(id = id, d = d, z = z)
}

# turn latent z into presence/absence of taxa on a tree -> unweighted UniFrac
make_tree <- function(ntax) {
  tr <- rtree(ntax)
  # edge x tip descendant indicator
  desc <- prop.part(tr)                      # internal nodes
  ntip <- Ntip(tr)
  E <- matrix(0, nrow(tr$edge), ntip)
  for (e in seq_len(nrow(tr$edge))) {
    ch <- tr$edge[e, 2]
    if (ch <= ntip) E[e, ch] <- 1 else E[e, desc[[ch - ntip]]] <- 1
  }
  list(E = E, bl = tr$edge.length, tr = tr)
}
unifrac_uw <- function(P, tree) {            # P: samples x taxa (0/1)
  B <- (P %*% t(tree$E)) > 0                 # samples x edges presence
  B <- B * 1; bl <- tree$bl
  sh <- B %*% (t(B) * bl)                    # shared branch length
  tot <- as.vector(B %*% bl)
  un <- outer(tot, tot, "+") - sh            # union
  1 - sh / un
}
z_to_uf <- function(z, tree, load, mu, depth = 3000) {
  eta <- sweep(z %*% load, 2, mu, "+")
  pr <- exp(eta); pr <- pr / rowSums(pr)
  P <- t(apply(pr, 1, function(q) rmultinom(1, depth, q) > 0)) * 1
  unifrac_uw(P, tree)
}

## ---------- one replicate ----------
one_rep <- function(sc, dist = c("euclid", "unifrac"), tree = NULL) {
  dist <- match.arg(dist)
  s <- sim_data(sc$nind, sc$k, sc$p, sc$m, sc$rho_d, sc$rho_y, sc$R2true)
  n <- sc$nind * sc$k
  G <- if (dist == "euclid") gram_euclid(s$z) else
    gower(z_to_uf(s$z, tree, sc$load, sc$mu))
  Z <- model.matrix(~ factor(s$id))[, -1, drop = FALSE]
  A <- model.matrix(~ 0 + factor(s$id)); A <- A / sc$k     # n x nind averaging
  Gm <- t(A) %*% G %*% A                                  # Gram of individual centroids
  Gm <- Gm - outer(rowMeans(Gm), colMeans(Gm), "+") + mean(Gm)
  dm <- t(A) %*% s$d
  B <- NPERM
  r_naive  <- test_marg(s$d, G, perm_free(n, B))
  r_block  <- test_marg(s$d, G, perm_block(sc$nind, sc$k, B))
  r_means  <- test_marg(dm, Gm, perm_free(sc$nind, B))
  r_within <- test_partial(s$d, Z, G, perm_within(sc$nind, sc$k, B))
  rbind(naive = unlist(r_naive), block = unlist(r_block),
        means = unlist(r_means), within = unlist(r_within))
}

run_scen <- function(sc, dist = "euclid", nrep = NREP) {
  tree <- NULL
  if (dist == "unifrac") {
    tree <- make_tree(sc$ntax)
    sc$load <- matrix(rnorm(sc$m * sc$ntax, sd = 0.8), sc$m)
    sc$mu <- rnorm(sc$ntax, sd = 2)
  }
  res <- mclapply(seq_len(nrep), function(r) one_rep(sc, dist, tree),
                  mc.cores = NCORE, mc.set.seed = TRUE)
  arr <- simplify2array(res)                               # method x stat x rep
  out <- data.frame(method = dimnames(arr)[[1]],
    reject = apply(arr[, "p", ] < 0.05, 1, mean),
    R2 = apply(arr[, "R2", ], 1, mean),
    adjR2 = apply(arr[, "adjR2", ], 1, mean),
    adjR2_sd = apply(arr[, "adjR2", ], 1, sd))
  cbind(dist = dist, as.data.frame(sc[c("nind", "k", "p", "m", "rho_d", "rho_y", "R2true")]),
        out, row.names = NULL)
}

## ---------- scenarios ----------
RNGkind("L'Ecuyer-CMRG")
base <- list(nind = 24, k = 4, p = 5, m = 20, rho_d = 0.5, rho_y = 0.5, R2true = 0, ntax = 300)
mod <- function(...) modifyList(base, list(...))

scen <- list()
# 1. null, grid of within-individual correlations
for (rd in c(0, 0.5, 0.8)) for (ry in c(0, 0.5, 0.8))
  scen[[length(scen) + 1]] <- mod(rho_d = rd, rho_y = ry)
# 2. null, number of samples per individual and number of diet variables
for (k in c(2, 8)) scen[[length(scen) + 1]] <- mod(k = k, rho_d = 0.8, rho_y = 0.8)
for (p in c(2, 10)) scen[[length(scen) + 1]] <- mod(p = p, rho_d = 0.8, rho_y = 0.8)
# 3. true diet effect (same within and between individuals)
for (R2 in c(0.05, 0.11)) for (rr in c(0, 0.5, 0.8))
  scen[[length(scen) + 1]] <- mod(R2true = R2, rho_d = rr, rho_y = rr)

res <- do.call(rbind, lapply(scen, run_scen, dist = "euclid"))
# 4. unweighted UniFrac on simulated taxa presence/absence (null only)
for (rr in c(0, 0.5, 0.8))
  res <- rbind(res, run_scen(mod(rho_d = rr, rho_y = rr), dist = "unifrac", nrep = 500))

res[] <- lapply(res, function(x) if (is.numeric(x)) round(x, 4) else x)
write.table(res, "sim_repeated_rda_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
