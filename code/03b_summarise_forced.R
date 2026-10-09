###############################################################
# 03b_summarise_forced.R
# Rebuilds the summaries and figures of 03_ensemble_validation.R from
# its saved results, without re-simulating:
#   Part B from results/03_forced_ensemble_raw.csv (the first version of
#     03 computed the SEs after the means had overwritten the columns,
#     giving NA);
#   Part A: flags degrees present in fewer than MIN_NETS networks and
#     redraws the figure without them.
#
# SEs are across networks: each raw row is one network's mean over
# REPS_F runs, so sd / sqrt(number of networks) includes both
# within- and between-network variation.
#
# Run from the repo root:   Rscript code/03b_summarise_forced.R
###############################################################

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

B_raw <- read.csv("results/03_forced_ensemble_raw.csv")
B_sum <- B_raw %>% group_by(tau, c, gamma) %>%
  summarise(n_nets = n(),
            se_imp = sd(sim_imp) / sqrt(n()), se_cases = sd(sim_cases) / sqrt(n()),
            across(c(sim_imp, sim_cases, th_imp, th_cases_avg, th_cases_node), mean),
            .groups = "drop") %>%
  mutate(z_imp = (sim_imp - th_imp) / se_imp,
         z_avg = (sim_cases - th_cases_avg) / se_cases,
         z_node = (sim_cases - th_cases_node) / se_cases)
write.csv(B_sum, "results/03_forced_ensemble.csv", row.names = FALSE)
print(as.data.frame(B_sum), digits = 3)

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
       title = sprintf("Forced importation, ensemble of %d networks: simulation (points) vs theory (lines)",
                       max(B_sum$n_nets))) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave("figures/03_forced_ensemble.pdf", figB, width = 11, height = 4.8)
cat("\nSaved results/03_forced_ensemble.csv and figures/03_forced_ensemble.pdf\n")

# ---- Part A: flag degrees seen in too few networks ---------------------------
# Very high degrees occur in only a few networks; there the "ensemble"
# is one or two specific men and the network-clustered SE is undefined
# (n_nets = 1) or unreliable. Rows are kept but flagged, and excluded
# from the z-score summary and the figure.
MIN_NETS <- 5L
A_sum <- read.csv("results/03_cluster_ensemble.csv") %>%
  mutate(enough_networks = n_nets >= MIN_NETS)
write.csv(A_sum, "results/03_cluster_ensemble.csv", row.names = FALSE)

cat(sprintf("\nPart A: %d of %d (network, tau, k) rows have k seen in >= %d networks\n",
            sum(A_sum$enough_networks), nrow(A_sum), MIN_NETS))
print(as.data.frame(A_sum %>% filter(enough_networks) %>%
  group_by(network, tau) %>%
  summarise(k_max_used = max(k), median_abs_z_ens = median(abs(z_ens)),
            median_abs_z_node = median(abs(z_node)),
            max_abs_dP0 = max(abs(sim_p0 - theory_p0)), .groups = "drop")), digits = 3)
print(as.data.frame(A_sum %>% filter(enough_networks) %>%
  select(network, tau, k, n_nets, runs, sim_mean, sim_se_net, theory_ens, z_ens, theory_node, z_node)),
  digits = 3)

figA <- ggplot(filter(A_sum, enough_networks), aes(k, sim_mean, color = factor(tau))) +
  geom_pointrange(aes(ymin = sim_mean - 2 * sim_se_net, ymax = sim_mean + 2 * sim_se_net), size = 0.25) +
  geom_line(aes(y = theory_ens), linewidth = 0.7) +
  geom_point(aes(y = theory_node), shape = 4, size = 2) +
  facet_wrap(~ network, scales = "free") +
  scale_color_viridis_d(name = expression(tau), end = 0.9) +
  labs(x = "Degree k of the imported case", y = "Mean secondary cases",
       title = sprintf("Network ensemble (degrees present in >= %d networks): simulation (+-2 SE), s(k) (line), node-level (x)",
                       MIN_NETS)) +
  theme_bw(base_size = 11)
ggsave("figures/03_cluster_ensemble.pdf", figA, width = 11, height = 4.8)
cat("Saved results/03_cluster_ensemble.csv (with enough_networks flag) and figures/03_cluster_ensemble.pdf\n")
