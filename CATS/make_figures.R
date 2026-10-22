# Regenerates the figures in figures/ (the examples in README.md).
# run from the CATS folder:  Rscript make_figures.R
library(CATS)
dir.create("figures", showWarnings = FALSE)

## Which design is better
png("figures/design.png", width = 600, height = 600)
rr <- seq(1.1, 1.4, by = 0.01)
res <- sapply(rr, function(r) {
  x <- cats(freq = 0.3, ncases = 2000, ncontrols = 2000, ncases2 = 2000, ncontrols2 = 2000,
            risk = r, pimarkers = 0.01, alpha = 5e-7)
  c(joint = x$P.joint, replication = x$P.rep.study, one = x$P.one.study)
})
plot(rr, res["one", ], type = "l", lwd = 2, col = 3, ylim = 0:1,
     xlab = "RR", ylab = "Power", main = "Power plot")
lines(rr, res["joint", ], lwd = 2, col = 1)
lines(rr, res["replication", ], lwd = 2, col = 2)
legend("bottomright", c("joint analysis", "replication design", "one stage design"),
       col = c(1, 2, 3), lwd = 2, bty = "n")
dev.off()

## Number of SNPs for replication
png("figures/replication_snps.png", width = 600, height = 600)
RR <- seq(1.1, 1.5, by = 0.025)
maf <- c(5, 10, 20, 50) / 100
power.J96 <- power.J1536 <- c()
for (m in maf) {
  J1 <- J2 <- c()
  for (r in RR) {
    J1 <- c(J1, cats(freq = m, ncases = 1500, ncontrols = 1500, ncases2 = 2000, ncontrols2 = 2000,
                     risk = r, pimarkers = 96 / 300000, alpha = 0.05 / 300000)$P.joint)
    J2 <- c(J2, cats(freq = m, ncases = 1500, ncontrols = 1500, ncases2 = 2000, ncontrols2 = 2000,
                     risk = r, pimarkers = 1536 / 300000, alpha = 0.05 / 300000)$P.joint)
  }
  power.J96 <- cbind(power.J96, J1)
  power.J1536 <- cbind(power.J1536, J2)
}
col <- seq_along(maf)
matplot(RR, power.J1536, type = "l", lty = 1, lwd = 2, col = col, ylim = 0:1, ylab = "Power",
        main = "Multiplicative model, 1500:1500, 2000:2000")
matlines(RR, power.J96, lty = 2, lwd = 2, col = col)
legend("bottomright", c(paste("MAF =", maf, ", rep = 1536"), paste("MAF =", maf, ", rep = 96")),
       col = col, lwd = 2, bty = "n", lty = rep(1:2, each = length(maf)))
dev.off()

## Heatmap of power
png("figures/heatmap.png", width = 600, height = 600)
rr <- seq(1, 2, by = 0.025)
sc <- super.cats(rr, by = length(rr), ncases = 765, ncontrols = 1274, ncases2 = 100, ncontrols2 = 100,
                 alpha = 0.001, prevalence = 0.01)
plot(sc, main = "power", file = NULL)
dev.off()

## Design and MAF
png("figures/design_maf.png", width = 600, height = 600)
rr <- seq(1, 2, by = 0.05)
maf <- c(0.05, 0.1, 0.2, 0.5)
c2 <- curve.cats(rr, maf, ncases = 600, ncontrols = 600, ncases2 = 600, ncontrols2 = 600,
                 alpha = 0.000001, prevalence = 0.01)
plot(c2, type = "One", main = "power2", ylab = "Power", xlab = "RR", file = NULL, col = 1:4)
lines.cats(c2, type = "Replication", lty = 3)
lines.cats(c2, type = "Joint", lty = 2)
lines.cats(c2, type = "First", lty = 4)
legend("left", c("One stage", "Joint", "Replication", "First stage"), lty = 1:4, bty = "n")
dev.off()
