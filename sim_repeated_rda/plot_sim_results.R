# Plots of the repeated-sample RDA simulations (white background PNG + PDF)
# Run from sim_repeated_rda/: Rscript --vanilla plot_sim_results.R
G <- read.delim("sim_repeated_rda_results.tsv")
C <- read.delim("sim_dubois_calibrated_results.tsv")
dir.create("plots", showWarnings = FALSE)

COL <- c(naive = "#2a78d6", block = "#eb6834", means = "#1baf7a")
LAB <- c(naive = "RDA (Dubois et al.)", block = "RDA (ID permutation)",
         means = "Individual means")
INK <- "#1b2228"; MUTED <- "#5b6670"; GRID <- "#e3e7ea"

# groups: list of list(label, sub, truth, pts = data.frame(m, y, lo, hi))
dotplot <- function(groups, ylim, ylab, refs = NULL, reflab = NULL, ticks = pretty(ylim),
                    fmt = function(x) format(x), legend_methods = NULL, main = NULL, legend_ncol = NULL, left_mar = 4.6) {
  ng <- length(groups)
  if (is.null(legend_methods)) legend_methods <- names(COL)[names(COL) %in% unlist(lapply(groups, function(g) g$pts$m))]
  nleg_rows <- if (is.null(legend_ncol)) 1 else ceiling(length(legend_methods) / legend_ncol)
  par(mar = c(4.2, left_mar, (if (is.null(main)) 3.2 else 4.4) + 1.1 * (nleg_rows - 1), 5.2), mgp = c(3, .6, 0), las = 1,
      col.axis = MUTED, fg = MUTED, family = "sans")
  plot(NA, xlim = c(.5, ng + .5), ylim = ylim, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       bty = "n", xaxs = "i")
  for (i in seq(2, ng, by = 2)) rect(i - .5, ylim[1], i + .5, ylim[2], col = "#f4f6f7", border = NA)
  abline(h = ticks, col = GRID, lwd = 1)
  axis(2, at = ticks, labels = fmt(ticks), lwd = 0, cex.axis = .9)
  mtext(ylab, side = 2, line = 3.3, las = 0, col = INK, cex = .95)
  for (k in seq_along(refs)) {
    abline(h = refs[k], lty = 2, col = "#8a949d", lwd = 1.3)
    if (!is.null(reflab) && nzchar(reflab[k]))
      text(ng + .58, refs[k], reflab[k], adj = c(0, .5), cex = .78, col = MUTED, xpd = NA)
  }
  for (i in seq_len(ng)) {
    g <- groups[[i]]; n <- nrow(g$pts); off <- (seq_len(n) - (n + 1) / 2) * min(.16, .8 / n)
    if (!is.null(g$truth)) segments(i - .4, g$truth, i + .4, g$truth, lwd = 3, col = INK)
    for (j in seq_len(n)) {
      p <- g$pts[j, ]
      if (!is.na(p$lo)) segments(i + off[j], max(ylim[1], p$lo), i + off[j], min(ylim[2], p$hi),
                                 col = COL[p$m], lwd = 2)
      points(i + off[j], min(max(p$y, ylim[1]), ylim[2]), pch = 21, bg = COL[p$m], col = "white",
             cex = 1.5, lwd = 1.2)
    }
    mtext(g$label, side = 1, at = i, line = .6, col = INK, cex = .85)
    if (!is.null(g$sub)) mtext(g$sub, side = 1, at = i, line = 1.7, col = MUTED, cex = if (is.numeric(g$sub)) .85 else if (ng > 6) .62 else .75)
  }
  if (!is.null(main)) mtext(main, side = 3, line = 2.6 + 1.1 * (nleg_rows - 1), adj = 0, col = INK, font = 2, cex = 1)
  par(xpd = NA)
  usr <- par("usr"); ytop <- usr[4] + diff(usr[3:4]) * 0.03
  legend(x = mean(usr[1:2]), y = ytop, xjust = 0.5, yjust = 0, bty = "n", pch = 21,
         pt.bg = COL[legend_methods], col = "white", pt.cex = 1.4, legend = LAB[legend_methods],
         text.col = INK, cex = .82, ncol = if (is.null(legend_ncol)) length(legend_methods) else legend_ncol,
         x.intersp = .8, text.width = NA)
  par(xpd = FALSE)
}

save_plot <- function(name, w, h, expr) {
  png(file.path("plots", paste0(name, ".png")), width = w, height = h, res = 160, bg = "white")
  expr(); dev.off()
  pdf(file.path("plots", paste0(name, ".pdf")), width = w / 160, height = h / 160, bg = "white")
  expr(); dev.off()
  svglite::svglite(file.path("plots", paste0(name, ".svg")), width = w / 160, height = h / 160,
                   bg = "white", system_fonts = list(sans = "Arial"))
  expr(); dev.off()
}

pts_of <- function(rs, methods, stat = "adjR2", band = c("sd", "q", "none")) {
  band <- match.arg(band)
  do.call(rbind, lapply(methods, function(m) {
    r <- rs[rs$method == m, ]
    lo <- hi <- NA
    if (band == "sd") { lo <- r$adjR2 - r$adjR2_sd; hi <- r$adjR2 + r$adjR2_sd }
    if (band == "q")  { lo <- r$adjR2_q05; hi <- r$adjR2_q95 }
    data.frame(m = sub("^(diet|dietGeo|nunavut)_", "", m), y = r[[stat]], lo = lo, hi = hi)
  }))
}
pct <- function(x) paste0(round(100 * x), "%")
f2 <- function(x) sprintf("%.2f", x)
four <- c("naive", "block", "means")
r2m <- c("naive", "means")   # R2 plots: ID permutation has the same R2 as the naive RDA

grid <- expand.grid(ry = c(0, .5, .8), rd = c(0, .5, .8))[, 2:1]
nullE <- function(rd, ry) G[G$dist == "euclid" & G$R2true == 0 & G$k == 4 & G$p == 5 &
                              G$rho_d == rd & G$rho_y == ry, ]
glab <- function(rd, ry) paste0(rd, " / ", ry)

# Fig 1: rejection, null grid
save_plot("sim_fig1_null_reject", 1500, 640, function() {
  groups <- lapply(seq_len(nrow(grid)), function(i) {
    rs <- nullE(grid$rd[i], grid$ry[i])
    list(label = glab(grid$rd[i], grid$ry[i]), pts = pts_of(rs, four, "reject", "none"))
  })
  dotplot(groups, c(0, 1), "Datasets with p < 0.05", refs = .05, reflab = "5%",
          ticks = seq(0, 1, .25), fmt = pct)
  mtext("Within-individual correlation: diet ρd / microbiome ρy",
        side = 1, line = 3, col = MUTED, cex = .85)
})

# Fig 2: adjR2, null grid (true diet effect = 0; vary within-individual correlations)
save_plot("sim_fig2_null_adjR2", 1500, 700, function() {
  groups <- lapply(seq_len(nrow(grid)), function(i) {
    rs <- nullE(grid$rd[i], grid$ry[i])
    list(label = grid$rd[i], sub = grid$ry[i], pts = pts_of(rs, r2m))
  })
  dotplot(groups, c(-.05, .15), expression("Adjusted " * R^2 * " for diet (mean ± SD)"), refs = 0,
          reflab = "truth = 0", ticks = seq(-.05, .15, .05), fmt = f2,
          main = "No diet effect (true diet–microbiome R² = 0)", left_mar = 10)
  # row labels, written once to the left of the group values
  x0 <- par("usr")[1] - 0.08
  mtext("cor(diet | ID)", side = 1, line = .6, at = x0, adj = 1, col = INK, cex = .85)
  mtext("cor(microbiome | ID)", side = 1, line = 1.7, at = x0, adj = 1, col = MUTED, cex = .85)
})

# Fig 3: k and p (both shown under every group)
save_plot("sim_fig3_k_p", 1500, 640, function() {
  sel <- list(c(1, 5), c(2, 5), c(4, 5), c(8, 5), c(4, 2), c(4, 5), c(4, 10))   # (k, p)
  groups <- lapply(sel, function(s) {
    rs <- G[G$dist == "euclid" & G$R2true == 0 & G$rho_d == .8 & G$rho_y == .8 & G$k == s[1] & G$p == s[2], ]
    list(label = s[1], sub = s[2], pts = pts_of(rs, r2m))
  })
  dotplot(groups, c(-.05, .25), expression("Adjusted " * R^2 * " (mean ± SD)"), refs = 0,
          reflab = "truth = 0", ticks = seq(0, .25, .05), fmt = f2, left_mar = 11)
  x0 <- par("usr")[1] - 0.08
  mtext("samples per person (k)", side = 1, line = .6, at = x0, adj = 1, col = INK, cex = .85)
  mtext("diet variables (p)", side = 1, line = 1.7, at = x0, adj = 1, col = MUTED, cex = .85)
})

# Fig 4: UniFrac, untransformed vs square-root (as in Dubois et al.), adjusted R2
save_plot("sim_fig4_unifrac", 1600, 680, function() {
  layout(matrix(1:2, 1))
  panel <- function(dist, title) {
    groups <- lapply(c(0, .5, .8), function(rr) {
      rs <- G[G$dist == dist & G$rho_d == rr & G$rho_y == rr, ]
      list(label = rr, sub = rr, pts = pts_of(rs, r2m))
    })
    dotplot(groups, c(-.1, .25), expression("Adjusted " * R^2 * " for diet (mean ± SD)"), refs = 0,
            reflab = "truth = 0", ticks = seq(-.1, .25, .05), fmt = f2, legend_ncol = 1, main = title,
            left_mar = 10)
    x0 <- par("usr")[1] - 0.08
    mtext("cor(diet | ID)", side = 1, line = .6, at = x0, adj = 1, col = INK, cex = .8)
    mtext("cor(microbiome | ID)", side = 1, line = 1.7, at = x0, adj = 1, col = MUTED, cex = .8)
  }
  panel("unifrac_paper", "Unweighted UniFrac")
  panel("unifrac_sqrt_paper", "Square-root UniFrac (Dubois et al.)")
})

# Fig 5: true effect
save_plot("sim_fig5_effect", 1500, 660, function() {
  eff <- expand.grid(rr = c(0, .5, .8), t = c(.05, .11))
  groups <- lapply(seq_len(nrow(eff)), function(i) {
    rs <- G[G$dist == "euclid" & G$R2true == eff$t[i] & G$rho_d == eff$rr[i], ]
    list(label = paste0("ρd = ρy = ", eff$rr[i]), sub = paste0("true R² = ", eff$t[i]), truth = eff$t[i],
         pts = pts_of(rs, c("naive", "means")))
  })
  dotplot(groups, c(0, .25), expression("Adjusted " * R^2 * " (mean ± SD)"), ticks = seq(0, .25, .05),
          fmt = f2, legend_methods = c("naive", "means"))
})

# Fig 6: Dubois-matched, adjusted R2 (UniFrac, as in the paper), three analyses
save_plot("sim_fig6_dubois_matched", 1900, 760, function() {
  layout(matrix(1:3, 1))
  rds <- sort(unique(C$rho_d))
  panel <- function(pre, methods, ref, title) {
    groups <- lapply(rds, function(rd) {
      rs <- C[C$dist == "unifrac" & C$R2true == 0 & C$rho_d == rd, ]
      list(label = rd, pts = pts_of(rs, paste0(pre, "_", methods), band = "q"))
    })
    dotplot(groups, c(-.12, .25), expression("Adjusted " * R^2 * " (mean, 5–95%)"), refs = c(0, ref),
            reflab = c("", paste0("Dubois\n", ref)), ticks = seq(-.1, .25, .05), fmt = f2,
            legend_methods = methods, legend_ncol = 1, main = title, left_mar = 7.5)
    x0 <- par("usr")[1] - 0.08
    mtext("cor(diet | ID)", side = 1, line = .6, at = x0, adj = 1, col = INK, cex = .8)
  }
  panel("diet", c("naive", "means"), .11, "Diet (all samples)")
  panel("dietGeo", c("naive", "means"), .064, "Diet given geography")
  panel("nunavut", c("naive"), .169, "Diet, Nunavut only")
})

# Fig 7: Dubois-matched false positives (UniFrac), diet | geography and Nunavut only
save_plot("sim_fig7_dubois_reject", 1600, 640, function() {
  layout(matrix(1:2, 1))
  rds <- sort(unique(C$rho_d))
  panel <- function(pre, methods, title) {
    groups <- lapply(rds, function(rd) {
      rs <- C[C$dist == "unifrac" & C$R2true == 0 & C$rho_d == rd, ]
      list(label = rd, pts = pts_of(rs, paste0(pre, "_", methods), "reject", "none"))
    })
    dotplot(groups, c(0, 1), "Datasets with p < 0.05", refs = .05, reflab = "5%", ticks = seq(0, 1, .25),
            fmt = pct, legend_methods = methods, legend_ncol = 1, main = title, left_mar = 7.5)
    x0 <- par("usr")[1] - 0.08
    mtext("cor(diet | ID)", side = 1, line = .6, at = x0, adj = 1, col = INK, cex = .8)
  }
  panel("dietGeo", c("naive", "block", "means"), "Diet given geography (truth = 0)")
  panel("nunavut", c("naive", "block"), "Diet, Nunavut only (truth = 0)")
})
cat("plots written to", normalizePath("plots"), "\n")
