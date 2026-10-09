###############################################################
# 09_core_network_check.R   (network structure only; no epidemics)
# NSFG 25-44 network with a female-sex-worker (FSW) core that absorbs the
# men's reporting surplus (Brewer et al. 2000), vs the trimmed network.
#
# For each FSW prevalence (per woman), client concentration q and network size, over
# N_SEEDS networks, reports:
#   - n_fsw and the IMPLIED FSW partner numbers (surplus / n_fsw), to compare
#     with Brewer et al.: median ~206, mean ~694 partners per year
#     (103 / 347 per 6 months, doubled);
#   - share of men who are clients (>= 1 FSW partner), to compare with the
#     GSS: 0.6-1% of men per year (2% in the largest cities; surveys likely
#     under-count clients); FSW partners per client and its size-biased
#     excess (what links FSW hubs to each other);
#   - kappa_m, kappa_f (all women incl. FSW), kappa_f of household women,
#     combined kappa, the transmissibility threshold tau_c = 1/kappa and R
#     at tau = 0.3 and 0.6 (average-degree theory; see caveat below);
#   - giant component share, number of FSW in it;
#   - bridges (p_bi = 0.02, k >= 1): how many, how many are FSW clients,
#     how many are in the giant component.
# Caveat: with a handful of hubs of degree ~100-1000, kappa_f is dominated
# by a few nodes and the average-degree theory is not expected to hold; the
# epidemic validation on this network is the next step.
#
# Run from the repo root:   Rscript code/09_core_network_check.R
# Output: results/09_core_network_check.csv (per network) and _summary.csv
# Runtime: ~5-15 min.
###############################################################

suppressPackageStartupMessages({ library(dplyr); library(igraph) })
source("code/network.R"); source("code/theory.R")
dir.create("results", showWarnings = FALSE)

FRACS   <- c(0.0002, 0.0005, 0.0017)       # FSW per woman (low / baseline / high)
QS      <- c(1, 2, 4, 8, 16)               # client concentration (weight (k - 1)^q)
SIZES   <- c(20000L, 50000L)               # people per sex
N_SEEDS <- 10L
P_BI    <- 0.02

describe <- function(net, label, frac, q = NA) {
  g  <- graph_from_edgelist(net$edges, directed = FALSE)
  g  <- add_vertices(g, max(0, net$N - vcount(g)))
  cm <- components(g); in_gc <- cm$membership == which.max(cm$csize)
  men <- seq_len(net$n_m)
  fsw <- if (is.null(net$fsw)) rep(FALSE, net$N) else net$fsw
  m_fsw <- vapply(men, function(i) sum(fsw[net$adj[[i]]]), numeric(1))   # FSW partners per man
  is_client <- m_fsw > 0
  pm <- realised_pmf(net, "M"); pf <- realised_pmf(net, "F")
  km <- pmf_kappa(pm); kf <- pmf_kappa(pf)
  kf_house <- { k <- net$deg[net$sex == "F" & !fsw]; sum(k * (k - 1)) / sum(k) }
  B  <- choose_bridges(net, P_BI)
  data.frame(
    network = label, fsw_frac = frac, q = q, n_per_sex = net$n_m, n_fsw = sum(fsw),
    fsw_deg_mean = if (any(fsw)) mean(net$deg[fsw]) else NA,
    fsw_deg_median = if (any(fsw)) median(net$deg[fsw]) else NA,
    fsw_deg_max = if (any(fsw)) max(net$deg[fsw]) else NA,
    pct_men_clients = 100 * mean(is_client),
    fsw_per_client = if (any(is_client)) mean(m_fsw[is_client]) else NA,
    # excess FSW partners of a client reached along an FSW partnership (size-biased);
    # this is what couples FSW hubs to each other
    fsw_per_client_excess = if (any(is_client)) sum(m_fsw * (m_fsw - 1)) / sum(m_fsw) else NA,
    max_fsw_per_client = max(m_fsw),
    kappa_m = km, kappa_f = kf, kappa_f_household = kf_house, kappa = sqrt(km * kf),
    tau_c = tau_threshold(km, kf),
    R_tau03 = R_two_step(0.3, km, kf), R_tau06 = R_two_step(0.6, km, kf),
    pct_in_gc = 100 * mean(in_gc), fsw_in_gc = sum(fsw & in_gc),
    n_bridges = length(B), bridges_clients = sum(is_client[B]), bridges_in_gc = sum(in_gc[B]))
}

rows <- list()
for (n in SIZES) {
  for (s in seq_len(N_SEEDS)) {
    rows[[length(rows) + 1]] <- describe(make_network("NSFG 25-44", n_per_sex = n, seed = s),
                                         "NSFG 25-44 (trimmed)", 0)
    for (fr in FRACS) for (q in QS) {
      rows[[length(rows) + 1]] <- describe(make_network_core("25-44", n_per_sex = n, fsw_frac = fr, q = q, seed = s),
                                           "NSFG 25-44 + FSW core", fr, q)
    }
  }
  cat(sprintf("size %d done\n", n))
}
res <- bind_rows(rows)
write.csv(res, "results/09_core_network_check.csv", row.names = FALSE)

summ <- res %>% group_by(network, fsw_frac, q, n_per_sex) %>%
  summarise(across(c(n_fsw, fsw_deg_mean, fsw_deg_median, fsw_deg_max, pct_men_clients,
                     fsw_per_client, fsw_per_client_excess, max_fsw_per_client,
                     kappa_m, kappa_f, kappa_f_household, kappa, tau_c, R_tau03, R_tau06,
                     pct_in_gc, fsw_in_gc, n_bridges, bridges_clients, bridges_in_gc), median),
            .groups = "drop")
write.csv(summ, "results/09_core_network_check_summary.csv", row.names = FALSE)
options(width = 200)
cat("\nMedians over", N_SEEDS, "networks (compare FSW degrees with Brewer: median ~206, mean ~694 per year;",
    "clients with GSS: 0.6-2% of men)\n")
print(as.data.frame(summ), digits = 3)
