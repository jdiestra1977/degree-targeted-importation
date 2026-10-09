# Condensation and optimal targeting of repeated forcing in heterogeneous networks

Code, simulation outputs and manuscript sources for the article submitted to
*Physical Review E* by Jose Herrera-Diestra (Department of Mathematics,
Tarleton State University).

An external process seeds a network repeatedly through a subset of entry
nodes. Each seeding event lands on entry node *i* with probability
proportional to w_i^gamma, where w is a heavy-tailed importance (the degree,
in the epidemic application). For tails p(w) ~ w^-a and responses ~ w^eta,
the targeted response diverges at gamma_1 = a - 1 - eta and the forcing
condenses at gamma_2 = a - 1. Depletion of nodes already reached makes the
burden-maximizing exponent gamma* positive; in large populations
gamma* ~ K/x for intense forcing and gamma_2 - gamma* ~ (a - 1)/ln(1/x) for
sparse forcing, where x is the number of seeding events per entry node. A
degree cut-off rounds the transitions below a crossover size n_c. The
dynamics are a stochastic SIR epidemic on bipartite networks calibrated to
published survey data; the pathogen is hypothetical.

## Repository layout

| Path | Contents |
|---|---|
| `code/network.R` | Bipartite configuration networks calibrated to the NSFG 2006-2008 partner distributions with a power-law tail; entry-node selection; optional commercial-sex core |
| `code/theory.R` | All analytical results: linear response, targeting gain, depletion, Gibbs measure, large-population limit, asymptotes, major outbreaks |
| `code/sim.R` | Discrete-time stochastic SIR matched exactly to the theory |
| `code/03_ensemble_validation.R`, `03b_summarise_forced.R` | Linear response and repeated forcing over network ensembles |
| `code/04_gamma_star.R` | Optimal exponent gamma* on survey-calibrated networks (theory) |
| `code/05_supercritical.R` | Major outbreaks above threshold (synthetic negative-binomial networks) |
| `code/07_scaling.R` | Targeting transition: participation ratio, growth exponent, finite-size scaling |
| `code/08_gamma_star_limit.R`, `08b_log_law.R` | gamma* in the large-population limit and the K/x asymptote |
| `code/09_core_network_check.R` | Calibration of the commercial-sex core |
| `code/10b_core_realistic_tau.R` | Illustration on the core network |
| `code/11_generalization_crossover.R`, `11b_smallx_check.R` | General response exponent eta, small-x asymptotics, cut-off crossover |
| `code/06_figures.R` | All figures (interactive, section by section) |
| `data/msm_incidence_NYC.csv` | Weekly 2022 mpox incidence in New York City, used only as the shape of the forcing profile |
| `results/` | Simulation and theory outputs (CSV) behind every number in the article |
| `figures/paper/` | Figures used in the article and Supplemental Material |
| `figures/schematic/fig_schematic.tex` | TikZ source of the model schematic (Fig. 1) |
| `manuscript/pre/main.tex` | Article (REVTeX 4.2) |
| `manuscript/pre/supplement.tex` | Supplemental Material |
| `manuscript/pre/pre_math_report.tex` | Extended methods and results: every derivation followed by its numerical evidence |
| `manuscript/pre/math/` | Assumptions, statements and proofs, shared by the Supplemental Material and the extended report |
| `manuscript/references.bib` | Bibliography |

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
- Weekly incidence: NYC Department of Health and Mental Hygiene,
  https://github.com/nychealth/monkeypox-data

## License

Code: MIT (see `LICENSE`). Cite the article and this repository (see
`CITATION.cff`) if you use them.
