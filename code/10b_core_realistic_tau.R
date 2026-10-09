###############################################################
# 10b_core_realistic_tau.R   (completes 10_multitype.R; final analysis)
# 10_multitype.R chose tau relative to the tree-based threshold, which is
# ~0.02 on the FSW-core network and turned out not to mark a real
# transition (the core is ~25 shared hubs, far from tree-like). This script
# repeats its Parts A and C at FIXED, realistic per-partnership
# transmissibilities tau = 0.1, 0.3, 0.6, 0.9 (mpox: roughly 0.3-0.8).
#
#   Part A  Single importations into random clients and random non-client
#           bridges: mean secondary cases (runs stopped at STOP_AT) and the
#           share reaching STOP_AT (1% of the population).
#   Part C  Continuous NYC-2022 importation into bridges, w = k,
#           gamma in {-2, 0, 1, 2, 3, 4, 6}, c = 10 and 100: expected cases
#           (capped at STOP_AT) and P(reaching STOP_AT).
# Same 10 baseline networks as 10_multitype.R (FSW 0.05%, q = 8, 50k/sex).
#
# Run from the repo root:   Rscript code/10b_core_realistic_tau.R
# Output: results/10b_*.csv     Runtime: ~45-90 min (3x the replicates of the
# first run: 300 single importations and 30 forced runs per point per network).
###############################################################

suppressPackageStartupMessages({ library(dplyr) })
source("code/network.R"); source("code/theory.R"); source("code/sim.R")

RHO     <- 1 / 14
R_DAY   <- 1 - exp(-RHO)
P_BI    <- 0.02
N_NET   <- 10L
STOP_AT <- 1000L
TAUS    <- c(0.1, 0.3, 0.6, 0.9)
N_SEEDS <- 300L                            # single importations per (network, tau, seed type)
REPS    <- 30L                             # forced runs per (network, tau, c, gamma)
dir.create("results", showWarnings = FALSE)
set.seed(2032)
alpha <- load_alpha("data/msm_incidence_NYC.csv")

A_rows <- list(); C_rows <- list()
for (s in seq_len(N_NET)) {
  net <- make_network_core(seed = 9100L + s)          # same networks as 10_multitype.R
  MT  <- multitype_M(net)
  B   <- choose_bridges(net, P_BI, seed = 9200L + s)
  kB  <- net$deg[B]
  clients    <- which(MT$m > 0)
  nonclientB <- B[MT$m[B] == 0]
  for (tau in TAUS) {
    p <- 1 - exp(-beta_for_tau(tau, RHO))
    # Part A
    for (who in c("client", "non-client bridge")) {
      pool  <- if (who == "client") clients else nonclientB
      seeds <- pool[sample.int(length(pool), N_SEEDS, replace = TRUE)]
      sec <- vapply(seeds, function(z) sim_cluster(net, z, p, R_DAY, stop_at = STOP_AT), numeric(1))
      A_rows[[length(A_rows) + 1]] <- data.frame(net = s, tau = tau, seed_type = who,
        mean_capped = mean(sec), median = median(sec), frac_reach = mean(sec >= STOP_AT))
    }
    # Part C
    for (cc in c(10, 100)) {
      beta_ext <- cc / (length(B) * sum(alpha))
      for (g in c(-2, 0, 1, 2, 3, 4, 6)) {
        pi <- target_weights(kB, g, "k")
        x  <- t(replicate(REPS, sim_forced(net, B, pi, alpha, beta_ext, p, R_DAY, stop_at = STOP_AT)))
        C_rows[[length(C_rows) + 1]] <- data.frame(net = s, tau = tau, c = cc, gamma = g,
          cases_capped = mean(x[, "total"]), p_reach = mean(x[, "major"]),
          importations = mean(x[, "importations"]))
      }
    }
  }
  cat(sprintf("network %d / %d done\n", s, N_NET))
}
A <- bind_rows(A_rows); C <- bind_rows(C_rows)
write.csv(A, "results/10b_single_importation_raw.csv", row.names = FALSE)
write.csv(C, "results/10b_targeting_raw.csv", row.names = FALSE)

As <- A %>% group_by(tau, seed_type) %>%
  summarise(mean_capped = mean(mean_capped), median = median(median),
            frac_reach = mean(frac_reach), .groups = "drop")
Cs <- C %>% group_by(tau, c, gamma) %>%
  summarise(se_cases = sd(cases_capped) / sqrt(n()),
            across(c(cases_capped, p_reach, importations), mean), .groups = "drop")
write.csv(As, "results/10b_single_importation.csv", row.names = FALSE)
write.csv(Cs, "results/10b_targeting.csv", row.names = FALSE)
options(width = 160)
print(as.data.frame(As), digits = 3)
print(as.data.frame(Cs), digits = 3)
cat("\nSaved results/10b_*.csv\n")
