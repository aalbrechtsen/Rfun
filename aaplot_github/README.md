# aaplot example

This page shows the plotting workflow from [`../aaplot.R`](../aaplot.R) using the example data bundled with [`../evalAdmix`](../evalAdmix).

The data are the `admixTjeck2` PLINK files, K=3 admixture proportions, and K=3 ancestral allele frequencies from `evalAdmix/data`. The residual correlation matrix was generated with:

```bash
cd evalAdmix
make
./evalAdmix \
  -plink data/admixTjeck2 \
  -fname data/admixTjeck2.3.P \
  -qname data/admixTjeck2.3.Q \
  -P 2 \
  -o ../aaplot_github/figures/admixTjeck2.corres.txt
```

The figures can be regenerated from the repository root with:

```bash
Rscript aaplot_github/make_figures.R
```

## Individual residual correlations above

The admixture barplot is shown below the pairwise residual correlation diamonds.

![Admixture proportions with evalAdmix residual correlations above](figures/aaplot_correlations_above.png)

## Individual residual correlations below

The same residual correlations can be placed below the admixture barplot by changing the `from` and `to` coordinates supplied to `addKey()` and `addCor()`.

![Admixture proportions with evalAdmix residual correlations below](figures/aaplot_correlations_below.png)

## Individual and population mean correlations

This view combines individual residual correlations above the admixture plot with population-mean residual correlations below it.

![Admixture proportions with individual and mean evalAdmix residual correlations](figures/aaplot_individual_and_mean_correlations.png)
