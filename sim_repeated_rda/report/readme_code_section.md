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
| `report/build_readme.py` | Builds this README from the same template, with the PNG figures (the Code section comes from `report/readme_code_section.md`) |
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
python3 report/build_readme.py
```

Needs R with `ape`, `svglite` and `parallel`; the dbRDA was checked against vegan 2.6-4 (`rda`, `varpart`, `RsquareAdj`). Tested with R 4.3.3.

Plot file names are from an earlier numbering. Report Figure 1 is `sim_fig2_null_adjR2`, Figure 2 `sim_fig3_k_p`, Figure 3 `sim_fig4_unifrac`, Figure 4 `sim_fig6_dubois_matched` and Figure 5 `sim_fig7_dubois_reject`.
