# Rfun

Random R code: small functions and packages for population genetics and association studies.

| What | Description |
|---|---|
| [CATS](CATS/README.md) | R package for the power of two-stage genome-wide association designs (one-stage, replication and joint analysis). Moved here from the popgen.dk software wiki. |
| [pcaFst.R](pcaFst.R) | PCA of a genotype matrix (Eigenstrat style, with Tracy-Widom test) and Weir & Cockerham Fst with variance components. See [pcaFst.md](pcaFst.md). Moved here from the popgen.dk software wiki. |
| [sim_repeated_rda](sim_repeated_rda/README.md) | Simulations showing that distance-based RDA with repeated samples per person, treated as independent, inflates p-values and adjusted R² for diet–microbiome associations; matched to Dubois et al. (2017). Includes an HTML report. |
| [aaplot.R](aaplot.R) | Plot admixture proportions and evalAdmix residual correlations. See the [aaplot example](aaplot_github/README.md). |
| [online.R](online.R) | Everyday helpers: colour palettes, QQ and Manhattan plots, PLINK file reading and writing, LocusZoom-style plots and more. |
| [NGS.R](NGS.R) | Genotype likelihoods from read counts and EM estimation of allele frequencies. |
| [latentVarRegression.R](latentVarRegression.R) | Regression with a latent variable (EM) for quantitative and binary traits, plus score tests. |
| [plot2dSFS.R](plot2dSFS.R) | Plot a 2D site frequency spectrum. |
| [power.R](power.R) | Allele frequencies in cases and controls for a given relative risk and prevalence. |

Most files are plain R scripts; load them with `source()`, for example:

```r
source("https://raw.githubusercontent.com/aalbrechtsen/Rfun/master/online.R")
```

CATS is a package with C code; see its [README](CATS/README.md) for how to install it.
