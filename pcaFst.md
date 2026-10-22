# PCA and Fst in R

Two functions in [pcaFst.R](pcaFst.R), moved here from the old popgen.dk software wiki.

```r
source("https://raw.githubusercontent.com/aalbrechtsen/Rfun/master/pcaFst.R")
```

## PCA (Eigensoft/Eigenstrat)

`eigenstrat(geno, maxMis = 0, minMaf = 0.01)` does a principal component analysis of a genotype matrix in the same way as Eigensoft/smartpca.

- `geno` is a **SNP x individual** matrix of genotypes coded 0, 1, 2 (`NA` for missing).
- Sites with a minor allele frequency below `minMaf`, or more than `maxMis` missing genotypes, are removed.

It returns the eigenvalues and eigenvectors, the number of SNPs and individuals used, and the Tracy-Widom statistic (`TW`) for the top eigenvalue. The test assumes there is no LD between the SNPs. `plot()` draws PC1 against PC2, and `print()` shows the test statistic.

```r
ind <- c(100, 100, 100)
snp <- 10000
freq <- cbind(runif(snp), runif(snp), runif(snp))
l <- lapply(1:length(ind), function(x) matrix(rbinom(snp * ind[x], 2, freq[, x]), snp))
geno <- do.call(cbind, l)
e <- eigenstrat(geno)
plot(e, col = rep(1:length(ind), ind), xlab = "PC1", ylab = "PC2")
```

## Fst (Weir and Cockerham 1984)

`WC84(x, pop)` estimates Fst like the faster `fstat` from the geneland package, but it also returns the variance components, both per SNP and as a joint estimate.

- `x` is an **individual x SNP** matrix of genotypes coded 0, 1, 2. SNPs with any missing genotype are removed.
- `pop` is the population label of each individual.

| Returned value | Meaning |
|---|---|
| `theta` | Fst per SNP |
| `theta_w` | Fst over all SNPs (ratio of sums) |
| `F` | F_IT per SNP |
| `f` | F_IS per SNP |
| `a`, `b`, `c` | variance between populations, between individuals within populations, and within individuals |
| `total` | total variance (`a + b + c`) |

```r
nsnp <- 10000
x <- matrix(rbinom(160 * nsnp, 2, 0.02), 160)
pop <- rep(1:4, 40)
res <- WC84(x, pop)
res$theta_w
```

Weir, B. S. and Cockerham, C. C. (1984). Estimating F-statistics for the analysis of population structure. *Evolution* 38(6): 1358–1370.

## LD pruning

The old wiki also linked to an R package for LD pruning (`pruning_0.51.tar.gz`, 2007). It no longer compiles with current R, so it is not included here.
