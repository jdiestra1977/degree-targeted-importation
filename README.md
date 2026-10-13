# Condensation and optimal targeting of repeated forcing in heterogeneous networks

Code, simulation outputs and manuscript sources for the article submitted to
*Physical Review E* by Jose Herrera-Diestra (Department of Mathematics,
Tarleton State University).

An external process seeds a network repeatedly through a subset of entry
nodes. Each seeding event lands on entry node *i* with probability
proportional to w_i^gamma, where w is a heavy-tailed importance (the degree,
in the epidemic application). For tails p(w) ~ w^-a and responses ~ w^eta,
the targeted response diverges at gamma_1 = a - 1 - eta and the forcing
condenses at gamma_2 = a - 1. Under repeated forcing the burden-maximizing
exponent gamma* is positive, and depletion of nodes already reached bounds it
by gamma_2 in large populations, with gamma* ~ K/x for intense forcing and
gamma_2 - gamma* ~ (a - 1)/ln(1/x) for sparse forcing, where x is the number
of seeding events per entry node. A degree cut-off rounds the transitions
below a crossover size n_c. The dynamics are a stochastic SIR epidemic on
bipartite networks calibrated to published survey data; the pathogen is
hypothetical.

## Repository layout

```
code/                 R code (shared libraries + numbered analysis scripts)
data/                 the one external input file (incidence curve)
results/              CSV outputs of the scripts, behind every number in the article
figures/paper/        figure files used in the article and Supplemental Material
figures/schematic/    TikZ source of Fig. 1
manuscript/pre/       LaTeX sources of the article and Supplemental Material
manuscript/references.bib
```

### Shared code (sourced by the scripts)

| File | Contents |
|---|---|
| `code/network.R` | Bipartite configuration networks calibrated to the NSFG 2006-2008 partner distributions with a power-law tail; entry-node (bridge) selection; optional commercial-sex core |
| `code/theory.R` | Analytical results: linear response, targeting, depletion, Gibbs-measure limit laws, asymptotics of gamma*, major outbreaks |
| `code/sim.R` | Discrete-time stochastic SIR, matched exactly to the theory; reads the forcing profile |

### Analysis scripts and what they produce

Scripts are numbered in the order they were written. Numbers 01, 02 and 10
belonged to preliminary analyses that are not part of this article, so they
are absent.

| Script | Article | Runtime |
|---|---|---|
| `03_ensemble_validation.R`, then `03b_summarise_forced.R` | Linear response and repeated forcing on network ensembles (Secs. III A, III C; Figs. 4a, S1, S7) | 20-40 min |
| `04_gamma_star.R` | Optimal exponent gamma* on survey-calibrated networks (Sec. III C; Fig. 4b,c) | 1-3 min |
| `05_supercritical.R` | Major outbreaks above threshold (Supplemental Material; Fig. S2) | 15-40 min |
| `07_scaling.R` | Targeting transition and finite-size scaling (Sec. III B; Figs. 3, S3, S4) | minutes |
| `08_gamma_star_limit.R`, `08b_log_law.R` | gamma* in the large-population limit, K/x and logarithmic asymptotes (Sec. III D; Figs. 4b, S5) | 10-20 min |
| `09_core_network_check.R` | Calibration of the commercial-sex core (Supplemental Material) | 5-15 min |
| `10b_core_realistic_tau.R` | Illustration on the core network (Sec. III G; Fig. S8) | 45-90 min |
| `11_generalization_crossover.R`, then `11b_smallx_check.R` | General response exponent eta, small-x asymptotics, cut-off crossover (Secs. III E, III F; Figs. 5, 6, S6) | 10-25 min |
| `06_figures.R` | All figures, from the saved results (interactive, section by section) | seconds |

### Figure files

File names follow the order in which the figures were made, not their final
numbers.

| Figure | File in `figures/paper/` |
|---|---|
| Fig. 1 | `fig_schematic.pdf` (from `figures/schematic/fig_schematic.tex`) |
| Fig. 2 | `fig2_targeting_gain.pdf` |
| Fig. 3 | `fig_transition.pdf` |
| Fig. 4 | `fig3_depletion_gamma_star.pdf` |
| Fig. 5 | `figS_gen_theta_eta.pdf` |
| Fig. 6 | `figS_gen_crossover.pdf` |
| Fig. S1 | `fig1_secondary_cases.pdf` |
| Fig. S2 | `fig4_major_outbreaks.pdf` |
| Fig. S3 | `figS_transition.pdf` |
| Fig. S4 | `figS_scaling_Y2_collapse.pdf`, `figS_scaling_mu.pdf` |
| Fig. S5 | `figS_scaling_gamma_star_limit.pdf` |
| Fig. S6 | `figS_gen_gamma_star_eta.pdf` |
| Fig. S7 | `figS1_network_calibration.pdf` |
| Fig. S8 | `fig5_core_network.pdf` |

Scripts 03-05 also write quick diagnostic plots to `figures/*.pdf`, and
`06_figures.R` writes a few diagnostic figures that are not in the article.
These are not tracked.

## Requirements

R 4.5.2 with dplyr 1.2.1, tidyr 1.3.2, ggplot2 4.0.2, patchwork 1.3.2,
scales 1.4.0 and igraph 2.2.3 (earlier recent versions should work). LaTeX
with REVTeX 4.2 and TikZ for the manuscript.

## Reproducing the results

All scripts run from the repository root (in RStudio, open
`degree-targeted-importation.Rproj`). The outputs are already in `results/`,
so the figures can be rebuilt without rerunning the simulations.

```
Rscript code/03_ensemble_validation.R      # ~20-40 min
Rscript code/03b_summarise_forced.R
Rscript code/04_gamma_star.R
Rscript code/05_supercritical.R            # ~15-40 min
Rscript code/07_scaling.R
Rscript code/08_gamma_star_limit.R
Rscript code/08b_log_law.R
Rscript code/09_core_network_check.R
Rscript code/10b_core_realistic_tau.R      # ~45-90 min
Rscript code/11_generalization_crossover.R
Rscript code/11b_smallx_check.R
```

Figures: run `code/06_figures.R` section by section, starting with section 0.

Manuscript, from `manuscript/pre/`:

```
pdflatex main && bibtex main && pdflatex main && pdflatex main
pdflatex supplement && bibtex supplement && pdflatex supplement && pdflatex supplement
```

The Supplemental Material includes its proofs from `pre_math_body.tex` and
`math/*.tex`.

Schematic, from `figures/schematic/`: `pdflatex fig_schematic.tex`, then copy
`fig_schematic.pdf` to `figures/paper/`.

## Data sources

All inputs are published and are either typed into `code/network.R` with
their source or included in `data/`.

- Past-year opposite-sex partners: Chandra et al. (2011), National Health
  Statistics Reports 36 (NSFG 2006-2008).
- Tail of the partner distribution: Liljeros et al. (2001), Nature 411, 907
  (cumulative exponents 2.31 for men and 2.54 for women).
- Bisexual identification among men (2.0%): Copen et al. (2016), National
  Health Statistics Reports 88.
- Commercial-sex core: Brewer et al. (2000), PNAS 97, 12385.
- Forcing profile (`data/msm_incidence_NYC.csv`): 17 weekly case counts of
  the 2022 mpox outbreak in New York City, from the NYC Department of Health
  and Mental Hygiene, https://github.com/nychealth/monkeypox-data. Only the
  shape of the curve is used (normalized to a maximum of 1, each week spread
  over 7 days).

## License

Code: MIT (see `LICENSE`). Cite the article and this repository (see
`CITATION.cff`) if you use them.
