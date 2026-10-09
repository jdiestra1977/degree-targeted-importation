###############################################################
# 05_supercritical.R
# Result (4): above threshold.
#   Part A  P(major outbreak | one importation into a man of degree k)
#           = 1 - h(k, u_f), u from the bipartite fixed point with the
#           exact (geometric) infectious-period distribution.
#   Part B  P(at least one major outbreak) under continuous NYC-2022
#           importation, as a function of gamma, vs
#           1 - prod_i [e^{-c pi_i} + (1 - e^{-c pi_i}) h(k_i, u_f)].
#
# The survey-calibrated (NSFG) networks are subcritical for any tau, so
# this uses NB(1.2, 0.7) networks (kappa ~ 3, threshold tau ~ 0.33) as a
# stand-in for a heterosexual network that contains a high-activity
# core. Ensemble of networks; a run counts as "major" once it reaches
# CUT_FRAC x (theoretical major-outbreak size for that network and tau)
# cases, where it is stopped. (A fixed cut-off misclassifies major
# outbreaks near threshold, where they are small.)
#
# Run from the repo root:   Rscript code/05_supercritical.R
# Output: results/05_*.csv, figures/05_*.pdf     Runtime: ~15-40 min.
###############################################################

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source("code/network.R"); source("code/theory.R"); source("code/sim.R")

RHO     <- 1 / 14
R_DAY   <- 1 - exp(-RHO)
P_BI    <- 0.02
CUT_FRAC <- 0.5
MIN_CUT  <- 100L
cutoff_for <- function(pm, pf, net, tau, fp)
  max(MIN_CUT, round(CUT_FRAC * major_outbreak_size(pm, pf, net$n_m, net$n_f, tau, fp$u_m, fp$u_f)))
NETTYPE <- "NB(1.2,0.7)"
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
set.seed(2028)

# ---- Part A: P(major | degree k) ---------------------------------------------

N_NET   <- 20L
PER_NET <- 50L
TAUS_A  <- c(0.4, 0.5, 0.6, 0.8)
K_MAX_A <- 10L

A_rows <- list()
for (s in seq_len(N_NET)) {
  net <- make_network(NETTYPE, seed = 7000L + s)
  pm <- realised_pmf(net, "M"); pf <- realised_pmf(net, "F")
  men <- seq_len(net$n_m)
  for (tau in TAUS_A) {
    p  <- 1 - exp(-beta_for_tau(tau, RHO))
    fp <- major_outbreak_fixed_point(pm, pf, p, R_DAY)
    cut <- cutoff_for(pm, pf, net, tau, fp)
    for (k in 1:K_MAX_A) {
      pool <- men[net$deg[men] == k]
      if (!length(pool)) next
      seeds <- pool[sample.int(length(pool), PER_NET, replace = TRUE)]
      maj <- vapply(seeds, function(z) sim_cluster(net, z, p, R_DAY, stop_at = cut) >= cut,
                    logical(1))
      A_rows[[length(A_rows) + 1]] <- data.frame(
        net = s, tau = tau, k = k, runs = PER_NET, majors = sum(maj), cutoff = cut,
        R = tau * sqrt(pmf_kappa(pm) * pmf_kappa(pf)),
        theory = p_major_from_man(k, fp$u_f, p, R_DAY))
    }
  }
  cat(sprintf("Part A: %d / %d networks\n", s, N_NET))
}
A_raw <- bind_rows(A_rows)
A_sum <- A_raw %>% group_by(tau, k) %>%
  summarise(R = mean(R), cutoff = median(cutoff), runs = sum(runs), sim = sum(majors) / runs,
            se = sqrt(sim * (1 - sim) / runs), theory = mean(theory), .groups = "drop") %>%
  mutate(z = (sim - theory) / sqrt(theory * (1 - theory) / runs))   # binomial SE under the theory (no 0/0 when sim = 0 or 1)
write.csv(A_raw, "results/05_pmajor_by_degree_raw.csv", row.names = FALSE)
write.csv(A_sum, "results/05_pmajor_by_degree.csv", row.names = FALSE)
print(as.data.frame(A_sum), digits = 3)

fA <- ggplot(A_sum, aes(k, sim, color = factor(tau))) +
  geom_pointrange(aes(ymin = sim - 2 * se, ymax = sim + 2 * se), size = 0.3) +
  geom_line(aes(y = theory), linewidth = 0.8) +
  scale_color_viridis_d(name = expression(tau), end = 0.85) +
  labs(x = "Degree k of the imported case", y = "P(major outbreak)",
       title = sprintf("NB(1.2,0.7) networks: simulation (points, %d networks) vs 1 - h(k, u_f) (lines)", N_NET)) +
  theme_bw(base_size = 11)
ggsave("figures/05_pmajor_by_degree.pdf", fA, width = 7.5, height = 4.8)

# ---- Part B: P(any major outbreak) under forcing, vs gamma -------------------

N_NET_F <- 20L
REPS_F  <- 25L
TAUS_B  <- c(0.4, 0.6)
C_LVLS  <- c(1, 5, 20)
GB      <- seq(-2, 6, by = 0.5)
alpha   <- load_alpha("data/msm_incidence_NYC.csv")

B_rows <- list()
for (s in seq_len(N_NET_F)) {
  net <- make_network(NETTYPE, seed = 8000L + s)
  B   <- choose_bridges(net, P_BI, seed = 9000L + s)
  S   <- net_summary(net, B)
  for (tau in TAUS_B) {
    p  <- 1 - exp(-beta_for_tau(tau, RHO))
    fp <- major_outbreak_fixed_point(S$pm, S$pf, p, R_DAY)
    cut <- cutoff_for(S$pm, S$pf, net, tau, fp)
    for (cc in C_LVLS) {
      beta_ext <- cc / (S$nB * sum(alpha))
      th <- p_any_major_forced(S$kB, GB, cc, fp$u_f, p, R_DAY, w = "k")
      for (gi in seq_along(GB)) {
        pi <- target_weights(S$kB, GB[gi], "k")
        x  <- t(replicate(REPS_F, sim_forced(net, B, pi, alpha, beta_ext, p, R_DAY, stop_at = cut)))
        B_rows[[length(B_rows) + 1]] <- data.frame(
          net = s, nB = S$nB, tau = tau, c = cc, gamma = GB[gi],
          runs = REPS_F, majors = sum(x[, "major"]), theory = th[gi], cutoff = cut)
      }
    }
  }
  cat(sprintf("Part B: %d / %d networks\n", s, N_NET_F))
}
B_raw <- bind_rows(B_rows)
B_sum <- B_raw %>% group_by(tau, c, gamma) %>%
  summarise(cutoff = median(cutoff), runs = sum(runs), sim = sum(majors) / runs, se = sqrt(sim * (1 - sim) / runs),
            theory = mean(theory), .groups = "drop")
B_sum <- B_sum %>% mutate(z = (sim - theory) / sqrt(pmax(theory * (1 - theory), 1e-12) / runs))
write.csv(B_raw, "results/05_forced_major_raw.csv", row.names = FALSE)
write.csv(B_sum, "results/05_forced_major.csv", row.names = FALSE)
print(as.data.frame(B_sum), digits = 3)

fB <- ggplot(B_sum, aes(gamma, sim, color = factor(tau))) +
  geom_pointrange(aes(ymin = sim - 2 * se, ymax = sim + 2 * se), size = 0.25) +
  geom_line(aes(y = theory), linewidth = 0.8) +
  facet_wrap(~ paste("c =", c)) +
  scale_color_viridis_d(name = expression(tau), end = 0.8) +
  labs(x = expression(gamma ~ "(w = k)"), y = "P(at least one major outbreak)",
       title = "Continuous importation, NB(1.2,0.7) ensemble: simulation (points) vs theory (lines)") +
  theme_bw(base_size = 11)
ggsave("figures/05_forced_major.pdf", fB, width = 11, height = 4.5)
cat("\nSaved results/05_*.csv and figures/05_*.pdf\n")
