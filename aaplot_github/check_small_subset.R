#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
script_path <- if (length(file_arg) == 0) {
  file.path("aaplot_github", "check_small_subset.R")
} else {
  sub("^--file=", "", file_arg[1])
}

script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
fig_dir <- file.path(script_dir, "figures")

source(file.path(script_dir, "aaplot_functions.R"))
source(file.path(repo_root, "evalAdmix", "visFuns.R"))

data_prefix <- file.path(repo_root, "evalAdmix", "data", "admixTjeck2")
cor_file <- file.path(fig_dir, "admixTjeck2.corres.txt")

pop <- read.table(paste0(data_prefix, ".fam"), stringsAsFactors = FALSE)
q <- read.table(paste0(data_prefix, ".3.Q"))
cor_mat <- as.matrix(read.table(cor_file))

pop_id <- pop[, 2]
ord <- orderInds(pop = as.vector(pop_id), q = q)
q <- q[ord, ]
pop_id <- pop_id[ord]
cor_mat <- cor_mat[ord, ord]

# Keep only a few individuals per population so separator problems are visible.
n_per_pop <- 5
keep <- unlist(lapply(unique(pop_id), function(p) which(pop_id == p)[seq_len(n_per_pop)]))
q <- q[keep, ]
pop_id <- pop_id[keep]
cor_mat <- cor_mat[keep, keep]

Q <- t(q)
pop_names <- unique(pop_id)
pop_id_int <- rep(seq_along(pop_names), table(pop_id)[pop_names])
pop_sep <- c(0, cumsum(table(pop_id_int)))
mean_pop_x <- pop_sep[-1] - table(pop_id_int) / 2
colorpal <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999"
)

png(file.path(fig_dir, "aaplot_small_subset_line_check.png"), width = 1600, height = 1200, res = 160)
par(mar = c(7, 4, 8, 2))
x <- barplot(Q, col = colorpal, space = 0, border = NA, axisnames = FALSE,
             ylab = "Admixture proportions", xlab = "", main = "",
             xpd = NA, cex.axis = 1.1, cex.lab = 1.4)
text(mean_pop_x, rep(-0.05, length(mean_pop_x)), unique(pop_id), xpd = TRUE, font = 2, cex = 1.2)
abline(v = pop_sep)
addCor(cor_mat, from = 1, to = 1.3, maxCor = 0.1, lines = 1.4, popID = pop_id)
dev.off()

message("Wrote small-subset line check to: ", file.path(fig_dir, "aaplot_small_subset_line_check.png"))
