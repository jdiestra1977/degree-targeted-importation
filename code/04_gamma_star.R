###############################################################
# 04_gamma_star.R   (theory only; no epidemic simulation)
# Result (3): the targeting exponent gamma* that maximises expected
# burden C(gamma) under continuous importation with depletion, as a
# function of spillover intensity c (expected importations without
# depletion) and transmissibility.
#
#   C(gamma) = sum_i (1 - exp(-c pi_i(gamma))) (1 + A k_i),  w = k
#
# Analytical facts checked numerically here:
#   - gamma* > 0 always (dC/dgamma at 0 = c e^{-c/|B|} A Cov_B(log k, k) > 0)
#   - gamma* -> large as c -> 0, gamma* -> 0 as c grows
#   - realised importations I(gamma) are maximal at gamma = 0
# gamma* equal to the grid maximum (30) means "unbounded" (C still rising).
#
# Bridge sets: 100 networks per degree distribution (NSFG 25-44,
# NSFG 15-44, and NB(1.2, 0.7) as a more heterogeneous reference),
# p_bi = 0.02, bridges = bisexual men with k >= 1.
# Transmissibility: NSFG networks are subcritical for any tau, so tau
# is set directly (0.3 / 0.6 / 0.9); for NB, tau is set per network so
# that R = tau sqrt(kappa_m kappa_f) = 0.3 / 0.6 / 0.9.
#
# Run from the repo root:   Rscript code/04_gamma_star.R
# Output: results/04_*.csv, figures/04_*.pdf     Runtime: ~1-3 min.
###############################################################

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
source("code/network.R"); source("code/theory.R")

RHO   <- 1 / 14
P_BI  <- 0.02
N_NET <- 100L
TYPES <- c("NSFG 25-44", "NSFG 15-44", "NB(1.2,0.7)")
LEVELS <- c(0.3, 0.6, 0.9)            # tau (NSFG) or target R (NB)
C_GRID <- 10^seq(0, 3.5, by = 0.1)
G_GRID <- seq(-1, 30, by = 0.05)
TAU_MAX <- exp(-RHO)                  # = 1 - r
TOL     <- 1e-9                # relative tolerance for "flat" curves
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

rows <- list(); bridge_deg <- list()
for (ty in TYPES) {
  for (s in seq_len(N_NET)) {
    net <- make_network(ty, seed = 5000L + s)
    B   <- choose_bridges(net, P_BI, seed = 6000L + s)
    S   <- net_summary(net, B)
    if (length(unique(S$kB)) < 2) next
    bridge_deg[[length(bridge_deg) + 1]] <- data.frame(type = ty, net = s, k = S$kB)
    PI  <- pi_matrix(S$kB, G_GRID, "k")
    i0  <- which.min(abs(G_GRID))
    for (lv in LEVELS) {
      tau <- if (ty == "NB(1.2,0.7)") lv / sqrt(S$kappa_m * S$kappa_f) else lv
      if (tau >= TAU_MAX) next
      A   <- A_factor(tau, S$kappa_m, S$kappa_f)
      if (is.na(A)) next
      val <- 1 + A * S$kB
      for (cc in C_GRID) {
        Q  <- 1 - exp(-cc * PI)
        Cg <- colSums(Q * val)
        Ig <- colSums(Q)
        j  <- which.max(Cg)
        # When c >> |B| every bridge is imported for gamma near 0 and C is
        # flat to machine precision; then gamma* is not identifiable
        # (targeting near uniform is irrelevant) and is reported as NA.
        flat <- (Cg[j] - Cg[i0]) <= TOL * Cg[i0]
        # Exact derivative at gamma = 0 (Corollary): c e^{-c/|B|} A Cov_B(log k, k)
        n  <- S$nB
        dC0 <- cc * exp(-cc / n) * A * mean((log(S$kB) - mean(log(S$kB))) * (S$kB - mean(S$kB)))
        rows[[length(rows) + 1]] <- data.frame(
          type = ty, net = s, level = lv, tau = tau,
          R = tau * sqrt(S$kappa_m * S$kappa_f), A = A,
          nB = n, kmax = max(S$kB), c = cc,
          gamma_star = if (flat) NA_real_ else G_GRID[j],
          gain_opt = Cg[j] / Cg[i0],                       # C(gamma*) / C(0)
          gain_hub = Cg[length(G_GRID)] / Cg[i0],          # C(gamma = 30) / C(0)
          dC_at_0_exact = dC0,
          I_max_at_0 = (max(Ig) - Ig[i0]) <= TOL * Ig[i0]) # TRUE if gamma = 0 maximises I (within TOL)
      }
    }
  }
  cat(sprintf("%s: done\n", ty))
}
res <- bind_rows(rows)
write.csv(res, "results/04_gamma_star_by_network.csv", row.names = FALSE)

cat("\nChecks over all networks, levels and c:\n")
cat(sprintf("  exact dC/dgamma at 0: min %.3g (should be >= 0; ~0 only when c >> |B|)\n", min(res$dC_at_0_exact)))
cat(sprintf("  gamma* identifiable in %d of %d cases; min identifiable gamma* = %.2f (should be > 0)\n",
            sum(!is.na(res$gamma_star)), nrow(res), min(res$gamma_star, na.rm = TRUE)))
cat(sprintf("  gamma = 0 maximises realised importations in %d of %d cases (should be all)\n",
            sum(res$I_max_at_0), nrow(res)))

summ <- res %>% group_by(type, level, c) %>%
  summarise(R_med = median(R), A_med = median(A),
            frac_identifiable = mean(!is.na(gamma_star)),
            gstar_med = median(gamma_star, na.rm = TRUE),
            gstar_q25 = quantile(gamma_star, 0.25, na.rm = TRUE),
            gstar_q75 = quantile(gamma_star, 0.75, na.rm = TRUE),
            gain_opt_med = median(gain_opt), gain_hub_med = median(gain_hub), .groups = "drop")
write.csv(summ, "results/04_gamma_star_summary.csv", row.names = FALSE)
print(as.data.frame(summ %>% filter(c %in% C_GRID[c(1, 11, 21, 31, 36)])), digits = 3)

lvl_lab <- function(ty) ifelse(ty == "NB(1.2,0.7)", "R", "tau")
summ <- summ %>% mutate(level_lbl = paste0(lvl_lab(type), " = ", level))

f1 <- ggplot(summ, aes(c, gstar_med, color = factor(level), fill = factor(level))) +
  geom_ribbon(aes(ymin = gstar_q25, ymax = gstar_q75), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  scale_x_log10() +
  facet_wrap(~ type, scales = "free_y") +
  scale_color_viridis_d(name = "tau (NSFG) / R (NB)", end = 0.85) +
  scale_fill_viridis_d(name = "tau (NSFG) / R (NB)", end = 0.85) +
  labs(x = "c, expected importations without depletion (log)",
       y = expression(gamma^"*" ~ "(median, IQR over networks)"),
       title = "Burden-maximising targeting exponent") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave("figures/04_gamma_star.pdf", f1, width = 11, height = 4.5)

f2 <- summ %>% pivot_longer(c(gain_opt_med, gain_hub_med), names_to = "what", values_to = "gain") %>%
  mutate(what = recode(what, gain_opt_med = "at gamma*", gain_hub_med = "extreme hub targeting (gamma = 30)")) %>%
  ggplot(aes(c, gain, color = factor(level), linetype = what)) +
  geom_line(linewidth = 0.8) + geom_hline(yintercept = 1, linetype = "dotted") +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ type, scales = "free_y") +
  scale_color_viridis_d(name = "tau (NSFG) / R (NB)", end = 0.85) +
  labs(x = "c (log)", y = "Burden relative to uniform targeting (log)", linetype = NULL,
       title = "How much targeting can raise (or lower) the burden") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave("figures/04_gain_vs_c.pdf", f2, width = 11, height = 4.5)

# Example C(gamma) curves for one network of each type, middle level
ex <- list()
for (ty in TYPES) {
  net <- make_network(ty, seed = 5001L); B <- choose_bridges(net, P_BI, seed = 6001L)
  S <- net_summary(net, B)
  tau <- if (ty == "NB(1.2,0.7)") 0.6 / sqrt(S$kappa_m * S$kappa_f) else 0.6
  A <- A_factor(tau, S$kappa_m, S$kappa_f)
  G <- seq(-1, 15, by = 0.05); PI <- pi_matrix(S$kB, G, "k")
  for (cc in c(1, 10, 100, 1000)) {
    Cg <- colSums((1 - exp(-cc * PI)) * (1 + A * S$kB))
    ex[[length(ex) + 1]] <- data.frame(type = ty, c = cc, gamma = G, rel = Cg / Cg[which.min(abs(G))])
  }
}
f3 <- ggplot(bind_rows(ex), aes(gamma, rel, color = factor(c))) +
  geom_line(linewidth = 0.8) + geom_hline(yintercept = 1, linetype = "dotted") +
  facet_wrap(~ type, scales = "free_y") +
  scale_color_viridis_d(name = "c", option = "C", end = 0.85) +
  labs(x = expression(gamma), y = "C(gamma) / C(0)",
       title = "Burden vs targeting for one network (tau = 0.6 NSFG; R = 0.6 NB)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave("figures/04_example_curves.pdf", f3, width = 11, height = 4.5)
cat("\nSaved results/04_*.csv and figures/04_*.pdf\n")
