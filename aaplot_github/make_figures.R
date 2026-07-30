#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
if (length(file_arg) == 0) {
  script_path <- file.path("aaplot_github", "make_figures.R")
} else {
  script_path <- sub("^--file=", "", file_arg[1])
}

script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
fig_dir <- file.path(script_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(script_dir, "aaplot_functions.R"))
source(file.path(repo_root, "evalAdmix", "visFuns.R"))

data_prefix <- file.path(repo_root, "evalAdmix", "data", "admixTjeck2")
cor_file <- file.path(fig_dir, "admixTjeck2.corres.txt")

if (!file.exists(cor_file)) {
  evaladmix_bin <- file.path(repo_root, "evalAdmix", "evalAdmix")
  if (!file.exists(evaladmix_bin)) {
    stop("evalAdmix binary not found. Run `make` in evalAdmix first.", call. = FALSE)
  }
  cmd <- sprintf(
    "%s -plink %s -fname %s -qname %s -P 2 -o %s",
    shQuote(evaladmix_bin),
    shQuote(data_prefix),
    shQuote(paste0(data_prefix, ".3.P")),
    shQuote(paste0(data_prefix, ".3.Q")),
    shQuote(cor_file)
  )
  status <- system(cmd)
  if (status != 0) {
    stop("evalAdmix failed while generating the residual correlation matrix.", call. = FALSE)
  }
}

pop <- read.table(paste0(data_prefix, ".fam"), stringsAsFactors = FALSE)
q <- read.table(paste0(data_prefix, ".3.Q"))
cor_mat <- as.matrix(read.table(cor_file))

pop_id <- pop[, 2]
ord <- orderInds(pop = as.vector(pop_id), q = q)
q <- q[ord, ]
pop_id <- pop_id[ord]
cor_mat <- cor_mat[ord, ord]

Q <- t(q)
N <- length(pop_id)
pop_names <- unique(pop_id)
Npop <- length(pop_names)
pop_id_int <- rep(seq_len(Npop), table(pop_id)[pop_names])
pop_sep <- c(0, cumsum(table(pop_id_int)))
mean_pop_x <- pop_sep[-1] - table(pop_id_int) / 2

colorpal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

png(file.path(fig_dir, "aaplot_correlations_above.png"), width = 2000, height = 2000, res = 200)
par(mar = c(5.1, 4.1, 10.1, 2.1))
x <- barplot(Q, col = colorpal, space = 0, border = NA, cex.axis = 1.2,
             cex.lab = 1.8, axisnames = FALSE, ylab = "Admixture proportions",
             xlab = "", main = "", cex.main = 1.5, xpd = NA)
text(mean_pop_x, rep(-0.05, length(mean_pop_x)), unique(pop_id), xpd = TRUE, font = 2, cex = 1.5)
abline(v = pop_sep)
text(mean(x), 1.35, "Residual correlations above admixture proportions", font = 2, cex = 1.8, xpd = TRUE)
addKey(from = 1, to = 1.3, N = ncol(Q), maxCor = 0.1)
addCor(cor_mat, from = 1, to = 1.3, maxCor = 0.1, lines = 0.2, popID = pop_id)
dev.off()

png(file.path(fig_dir, "aaplot_correlations_below.png"), width = 2000, height = 2000, res = 200)
par(mar = c(11.1, 4.1, 4.1, 2.1))
x <- barplot(Q, col = colorpal, space = 0, border = NA, cex.axis = 1.2,
             cex.lab = 1.8, axisnames = FALSE, ylab = "Admixture proportions",
             xlab = "", main = "", cex.main = 1.5, xpd = NA)
text(mean_pop_x, rep(1.05, length(mean_pop_x)), unique(pop_id), xpd = TRUE, font = 2, cex = 1.5)
abline(v = pop_sep)
addKey(from = 0, to = -0.3, N = ncol(Q), maxCor = 0.1)
addCor(cor_mat, from = 0, to = -0.3, maxCor = 0.1, lines = 0.2, popID = pop_id)
dev.off()

png(file.path(fig_dir, "aaplot_individual_and_mean_correlations.png"), width = 2000, height = 2000, res = 200)
par(mar = c(9.1, 4.1, 8.1, 2.1))
x <- barplot(Q, col = colorpal, space = 0, border = NA, cex.axis = 1.2,
             cex.lab = 1.8, axisnames = FALSE, ylab = "Admixture proportions",
             xlab = "", main = "", cex.main = 1.5, xpd = NA)
abline(v = pop_sep)
addKey(from = 1, to = 1.3, N = ncol(Q), maxCor = 0.1)
addCor(cor_mat, from = 1, to = 1.3, maxCor = 0.1, lines = 0.2, popID = pop_id)
addCor(cor_mat, from = -0.1, to = -0.4, maxCor = 0.1, lines = 0.2, popID = pop_id, meanCor = TRUE)
text(0, -0.2, "Mean\ncorrelation", adj = 0, xpd = TRUE)
text(mean_pop_x, rep(-0.05, length(mean_pop_x)), unique(pop_id), xpd = TRUE, font = 2, cex = 1.5)
dev.off()

message("Wrote figures to: ", fig_dir)
