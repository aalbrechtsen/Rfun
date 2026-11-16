# Repeated samples in diet–microbiome RDA

Simulations of what happens when distance-based redundancy analysis (RDA) treats repeated samples from the same people as independent, as in Dubois et al. (2017, *Microbiome* 5:151, [PMC5689144](https://pmc.ncbi.nlm.nih.gov/articles/PMC5689144/)). They related diet to square-root unweighted UniFrac distances using 128 samples from 24 people, with free permutations, and reported adjusted R² of 0.11 (diet), 0.064 (diet given geography) and 0.169 (Nunavut only).

When both diet and the microbiome are stable within people, cor(diet | ID) > 0 and cor(microbiome | ID) > 0, such an RDA finds a diet that has no effect significant in nearly every dataset. Its adjusted R² is also biased upwards. With the paper's design, values like those reported arise with no diet effect at all. Permuting whole individuals gives correct p-values, and an RDA of individual means gives unbiased R².

The full write-up with figures is in [report/rda_repeated_report.html](report/rda_repeated_report.html) (download and open in a browser).

## Methods compared

| Name in figures | What it is |
|---|---|
| RDA (Dubois et al.) | RDA of diet on all samples, rows permuted freely, as in the paper |
| RDA (ID permutation) | Same RDA and R², but permutations move whole individuals (only among people with the same number of samples) |
| Individual means | One row per person: the centroid of each person's samples in distance space against their average diet (n = 24) |

## Files

| File | Makes |
|---|---|
| `sim_repeated_rda_lib.R` | Shared functions: dbRDA (McArdle & Anderson), adjusted R², permutation tests, unweighted UniFrac |
| `sim_repeated_rda.R` | Report Figures 1–2: 24 people, equal samples per person, Euclidean distance |
| `sim_k1.R` | One sample per person (k = 1) for Figure 2 |
| `sim_unifrac_paper.R` | Figure 3: the paper's sample sizes (15 × 4, 5 × 8, 4 × 7 samples) and 16 diet variables, UniFrac and square-root UniFrac |
| `sim_unifrac_sqrt.R` | Square-root UniFrac with 4 samples per person and 5 diet variables (not shown in the report) |
| `sim_dubois_calibrated.R` | Figures 4–5 and Table 1: matched to Dubois et al. (geography, Nunavut only, individual 45% and geography 4% of UniFrac variation) |
| `plot_sim_results.R` | All figures as PNG, PDF and SVG in `plots/` |
| `report/build_report.py` | Builds `report/rda_repeated_report.html` from `report/report_template.html`, the SVG figures and the results tables |
| `*_results.tsv` | Simulation results (rejection rate, mean R², adjusted R² and its spread per scenario and method) |

## Running

The scripts use `parallel::mclapply` with up to 80 cores. Set `OPENBLAS_NUM_THREADS=1` so the workers don't also start threaded BLAS. Run them from this folder in this order, because the later ones add rows to `sim_repeated_rda_results.tsv`:

```bash
OPENBLAS_NUM_THREADS=1 Rscript --vanilla sim_repeated_rda.R
OPENBLAS_NUM_THREADS=1 Rscript --vanilla sim_k1.R
OPENBLAS_NUM_THREADS=1 Rscript --vanilla sim_unifrac_sqrt.R
OPENBLAS_NUM_THREADS=1 Rscript --vanilla sim_unifrac_paper.R
OPENBLAS_NUM_THREADS=1 Rscript --vanilla sim_dubois_calibrated.R
Rscript --vanilla plot_sim_results.R
python3 report/build_report.py
```

Needs R with `ape`, `svglite` and `parallel`; the dbRDA was checked against vegan 2.6-4 (`rda`, `varpart`, `RsquareAdj`). Tested with R 4.3.3.

Plot file names are from an earlier numbering. Report Figure 1 is `sim_fig2_null_adjR2`, Figure 2 `sim_fig3_k_p`, Figure 3 `sim_fig4_unifrac`, Figure 4 `sim_fig6_dubois_matched` and Figure 5 `sim_fig7_dubois_reject`.
