###############################################################
# 03_ensemble_validation.R
# Results (1)-(3) validated over an ENSEMBLE of independently built
# networks (the analytical formulas describe the configuration-model
# ensemble, not one particular network).
#
# Part A  Secondary cases per importation, s(k) = A k.
#         N_NET networks per survey distribution; in each, PER_NET runs
#         per degree k seeded in random men of degree k.
#         Theory: (i) ensemble formula with kappas from the pooled degree
#         distribution of all networks; (ii) node-level formula using
#         each seed's actual partners (same pooled kappas beyond them).
#         SEs: naive (all runs independent) and network-clustered
#         (between-network spread of per-network means).
#
# Part B  Continuous importation with depletion (NYC 2022 forcing),
#         NSFG 25-44, bridges p_bi = 0.02 with k >= 1, weights w = k.
#         For each network the formulas use that network's bridges and
#         kappas; predictions and simulations are averaged over networks.
#         Theory: average-degree C(gamma) and node-level C(gamma).
#
# Run from the repo root:   Rscript code/03_ensemble_validation.R
# Output: results/03_*.csv, figures/03_*.pdf     Runtime: ~20-40 min.
###############################################################

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
source("code/network.R"); source("code/theory.R"); source("code/sim.R")

RHO   <- 1 / 14
R_DAY <- 1 - exp(-RHO)
P_BI  <- 0.02
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
set.seed(2027)

# ---- Part A: cluster sizes over a network ensemble ---------------------------

N_NET   <- 100L
PER_NET <- 30L
TAUS    <- c(0.2, 0.4, 0.6, 0.8, 0.9)
NETS    <- c("NSFG 25-44", "NSFG 15-44")
NET_SEEDS <- 1000L + seq_len(N_NET)

A_rows <- list()
for (nm in NETS) {
  # Pass 1: pooled degree distribution of the ensemble (networks are
  # rebuilt identically in pass 2 from the same seeds)
  dm <- integer(0); df <- integer(0); kap <- matrix(NA, N_NET, 2)
  for (s in seq_len(N_NET)) {
    net <- make_network(nm, seed = NET_SEEDS[s])
    dm <- c(dm, net$deg[net$sex == "M"]); df <- c(df, net$deg[net$sex == "F"])
    kap[s, ] <- c(pmf_kappa(realised_pmf(net, "M")), pmf_kappa(realised_pmf(net, "F")))
  }
  km <- pmf_kappa(tabulate(dm + 1L) / length(dm))
  kf <- pmf_kappa(tabulate(df + 1L) / length(df))
  cat(sprintf("\n%s ensemble (%d networks): pooled kappa_m %.3f, kappa_f %.3f, tau threshold %.3f\n",
              nm, N_NET, km, kf, tau_threshold(km, kf)))
  cat(sprintf("  per-network kappa_m range %.2f-%.2f, kappa_f range %.2f-%.2f\n",
              min(kap[, 1]), max(kap[, 1]), min(kap[, 2]), max(kap[, 2])))
  write.csv(data.frame(network = nm, net_seed = NET_SEEDS, kappa_m = kap[, 1], kappa_f = kap[, 2]),
            sprintf("results/03_kappa_by_network_%s.csv", gsub(" ", "_", nm)), row.names = FALSE)

  # Pass 2: simulations
  for (s in seq_len(N_NET)) {
    net <- make_network(nm, seed = NET_SEEDS[s])
    men <- seq_len(net$n_m)
    for (tau in TAUS) {
      if (tau^2 * km * kf >= 1) next
      p <- 1 - exp(-beta_for_tau(tau, RHO))
      for (k in sort(unique(net$deg[men][net$deg[men] >= 1]))) {
        pool  <- men[net$deg[men] == k]
        seeds <- pool[sample.int(length(pool), PER_NET, replace = TRUE)]
        sec   <- vapply(seeds, function(z) sim_cluster(net, z, p, R_DAY), numeric(1))
        A_rows[[length(A_rows) + 1]] <- data.frame(
          network = nm, net_seed = NET_SEEDS[s], tau = tau, k = k, n = PER_NET,
          sum = sum(sec), sumsq = sum(sec^2), n0 = sum(sec == 0),
          node_theory_sum = sum(mean_secondary_node(net, seeds, tau, km, kf)))
      }
    }
    if (s %% 10 == 0) cat(sprintf("  %s: %d / %d networks done\n", nm, s, N_NET))
  }
  assign(paste0("KAP_", gsub("[^A-Za-z0-9]", "_", nm)), c(km, kf))
}
A_raw <- bind_rows(A_rows)
write.csv(A_raw, "results/03_cluster_ensemble_raw.csv", row.names = FALSE)

kap_of <- function(nm) get(paste0("KAP_", gsub("[^A-Za-z0-9]", "_", nm)))
A_sum <- A_raw %>%
  group_by(network, tau, k) %>%
  summarise(n_nets = n_distinct(net_seed), runs = sum(n),
            sim_mean = sum(sum) / runs,
            sim_se_naive = sqrt((sum(sumsq) / runs - sim_mean^2) / runs),
            sim_se_net = if (n_nets > 1) sd(sum / n) / sqrt(n_nets) else NA_real_,
            sim_p0 = sum(n0) / runs,
            theory_node = sum(node_theory_sum) / runs, .groups = "drop") %>%
  rowwise() %>%
  mutate(theory_ens = mean_secondary(k, tau, kap_of(network)[1], kap_of(network)[2]),
         p = 1 - exp(-beta_for_tau(tau, RHO)),
         theory_p0 = h_fun(k, 0, p, R_DAY)) %>%
  ungroup() %>%
  mutate(z_ens = (sim_mean - theory_ens) / sim_se_net,
         z_node = (sim_mean - theory_node) / sim_se_net)
write.csv(A_sum, "results/03_cluster_ensemble.csv", row.names = FALSE)
print(as.data.frame(A_sum %>% select(network, tau, k, n_nets, runs, sim_mean, sim_se_net,
                                     theory_ens, z_ens, theory_node, z_node, sim_p0, theory_p0)),
      digits = 3)

figA <- ggplot(A_sum, aes(k, sim_mean, color = factor(tau))) +
  geom_pointrange(aes(ymin = sim_mean - 2 * sim_se_net, ymax = sim_mean + 2 * sim_se_net), size = 0.25) +
  geom_line(aes(y = theory_ens), linewidth = 0.7) +
  geom_point(aes(y = theory_node), shape = 4, size = 2) +
  facet_wrap(~ network, scales = "free") +
  scale_color_viridis_d(name = expression(tau), end = 0.9) +
  labs(x = "Degree k of the imported case", y = "Mean secondary cases",
       title = sprintf("Ensemble of %d networks: simulation (+-2 network-clustered SE), s(k) (line), node-level (x)", N_NET)) +
  theme_bw(base_size = 11)
ggsave("figures/03_cluster_ensemble.pdf", figA, width = 11, height = 4.8)

# ---- Part B: forced importation over a network ensemble ----------------------

N_NET_F <- 50L
REPS_F  <- 10L
C_LVLS  <- c(10, 50, 200)
TAUS_B  <- c(0.6, 0.9)
GB      <- seq(-2, 6, by = 0.5)
alpha   <- load_alpha("data/msm_incidence_NYC.csv")

B_rows <- list()
for (s in seq_len(N_NET_F)) {
  net <- make_network("NSFG 25-44", seed = 2000L + s)
  B   <- choose_bridges(net, P_BI, seed = 3000L + s)
  S   <- net_summary(net, B)
  for (tau in TAUS_B) {
    p  <- 1 - exp(-beta_for_tau(tau, RHO))
    sB <- mean_secondary_node(net, B, tau, S$kappa_m, S$kappa_f)
    for (cc in C_LVLS) {
      beta_ext <- cc / (S$nB * sum(alpha))
      th_avg  <- forced_expectations(S$kB, GB, cc, tau, S$kappa_m, S$kappa_f, w = "k")
      th_node <- forced_expectations_node(S$kB, sB, GB, cc, w = "k")
      for (gi in seq_along(GB)) {
        pi <- target_weights(S$kB, GB[gi], "k")
        x  <- t(replicate(REPS_F, sim_forced(net, B, pi, alpha, beta_ext, p, R_DAY)))
        B_rows[[length(B_rows) + 1]] <- data.frame(
          net_seed = 2000L + s, nB = S$nB, kmax_B = max(S$kB), tau = tau, c = cc, gamma = GB[gi],
          sim_imp = mean(x[, "importations"]), sim_cases = mean(x[, "total"]),
          th_imp = th_avg[gi, "importations"], th_cases_avg = th_avg[gi, "cases"],
          th_cases_node = th_node[gi, "cases"])
      }
    }
  }
  if (s %% 5 == 0) cat(sprintf("Part B: %d / %d networks done\n", s, N_NET_F))
}
B_raw <- bind_rows(B_rows)
write.csv(B_raw, "results/03_forced_ensemble_raw.csv", row.names = FALSE)

# (SEs first: summarise() works sequentially, so computing the means
#  first would overwrite sim_imp / sim_cases before sd() sees them)
B_sum <- B_raw %>% group_by(tau, c, gamma) %>%
  summarise(se_imp = sd(sim_imp) / sqrt(n()), se_cases = sd(sim_cases) / sqrt(n()),
            across(c(sim_imp, sim_cases, th_imp, th_cases_avg, th_cases_node), mean),
            .groups = "drop")
write.csv(B_sum, "results/03_forced_ensemble.csv", row.names = FALSE)
print(as.data.frame(B_sum), digits = 3)

# gamma* per network (simulation grid vs node-level theory grid)
gstar <- B_raw %>% group_by(net_seed, tau, c) %>%
  summarise(gstar_sim = gamma[which.max(sim_cases)], gstar_node = gamma[which.max(th_cases_node)],
            gstar_avg = gamma[which.max(th_cases_avg)], .groups = "drop")
write.csv(gstar, "results/03_gamma_star_by_network.csv", row.names = FALSE)

figB <- B_sum %>%
  pivot_longer(c(th_cases_avg, th_cases_node), names_to = "theory", values_to = "th") %>%
  mutate(theory = recode(theory, th_cases_avg = "average-degree", th_cases_node = "node-level")) %>%
  ggplot(aes(gamma)) +
  geom_line(aes(y = th, color = factor(tau), linetype = theory), linewidth = 0.7) +
  geom_pointrange(data = B_sum, aes(y = sim_cases, ymin = sim_cases - 2 * se_cases,
                                    ymax = sim_cases + 2 * se_cases, color = factor(tau)), size = 0.25) +
  facet_wrap(~ paste("c =", c), scales = "free_y") +
  scale_color_viridis_d(name = expression(tau), end = 0.8) +
  labs(x = expression(gamma ~ "(w = k)"), y = "Expected total cases", linetype = "Theory",
       title = sprintf("Forced importation, ensemble of %d networks: simulation (points) vs theory (lines)", N_NET_F)) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave("figures/03_forced_ensemble.pdf", figB, width = 11, height = 4.8)
cat("\nSaved results/03_*.csv and figures/03_*.pdf\n")
