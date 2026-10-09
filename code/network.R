###############################################################
# network.R
# Bipartite (men-women) configuration-model contact networks with
# survey-calibrated degree distributions, and bridge selection.
#
# Degree = number of opposite-sex partners in the past 12 months.
#
# Data:
#   NSFG 2006-2008, Chandra et al. (2011) National Health Statistics
#   Reports no. 36, Tables 1 (women) and 2 (men). Percent with
#   0 / 1 / 2 / 3 / 4+ opposite-sex partners in the past 12 months
#   ("0" = never had opposite-sex contact + had it but not in the last
#   12 months; "did not report" dropped and the rest renormalised).
#   The 4+ class is spread over k = 4..KMAX with a power-law tail,
#   P(k) ~ k^-(a + 1), where a is the CUMULATIVE exponent for past-year
#   partners in Liljeros et al. (2001) Nature 411:907, Fig. 1a
#   (men a = 2.31 +- 0.2, k > 5; women a = 2.54 +- 0.2, k > 4).
#
# Balancing: men report more partners than women, so stub totals
# differ. The side with the excess loses only "extra" stubs (a node of
# degree k can lose at most k - 1, chosen in proportion to k - 1), so
# nodes with 0 or 1 partner keep their reported degree.
###############################################################

KMAX <- 25L

NSFG_PCT <- list(
  `25-44` = list(M = c(2.3 + 6.3, 75.4, 7.0, 2.9, 5.2),
                 F = c(1.6 + 6.6, 82.3, 5.4, 1.7, 1.9)),
  `15-44` = list(M = c(11.2 + 6.6, 62.5, 8.6, 3.9, 6.0),
                 F = c(11.0 + 6.1, 69.0, 7.6, 2.5, 2.9))
)
LILJEROS_CUM_ALPHA <- c(M = 2.31, F = 2.54)

# pmf over k = 0..kmax (element i is P(k = i - 1))
nsfg_pmf <- function(pct, alpha_cum, kmax = KMAX) {
  p <- pct / sum(pct)
  w <- (4:kmax)^(-(alpha_cum + 1))
  c(p[1:4], p[5] * w / sum(w))
}
nb_pmf <- function(mu, size, kmax = KMAX) {          # truncated by capping at kmax
  p <- dnbinom(0:kmax, size = size, mu = mu)
  p[kmax + 1] <- p[kmax + 1] + pnbinom(kmax, size = size, mu = mu, lower.tail = FALSE)
  p
}
sample_pmf <- function(n, pmf) as.integer(sample.int(length(pmf), n, replace = TRUE, prob = pmf) - 1L)

# Remove `excess` stubs from deg, at most k - 1 per node, prob. proportional to k - 1
trim_extra <- function(deg, excess) {
  extra <- rep.int(seq_along(deg), times = pmax(deg - 1L, 0L))
  if (length(extra) < excess) stop("Not enough multi-partner stubs to balance.")
  deg - tabulate(extra[sample.int(length(extra), excess)], nbins = length(deg))
}

# Build the network. Nodes 1..n_m are men, n_m+1..n_m+n_f women.
# Multi-edges (the same couple paired twice) are collapsed to one
# partnership, so realised degrees can be marginally below the drawn ones.
build_network <- function(n_m, n_f, pmf_m, pmf_f, seed = 1L) {
  set.seed(seed)
  dm <- sample_pmf(n_m, pmf_m)
  df <- sample_pmf(n_f, pmf_f)
  if (sum(dm) > sum(df)) dm <- trim_extra(dm, sum(dm) - sum(df))
  if (sum(df) > sum(dm)) df <- trim_extra(df, sum(df) - sum(dm))

  u <- rep.int(seq_len(n_m), dm)
  v <- n_m + sample(rep.int(seq_len(n_f), df))
  e <- unique(cbind(u, v))                                   # collapse multi-edges
  N <- n_m + n_f
  adj <- split(c(e[, 2], e[, 1]), factor(c(e[, 1], e[, 2]), levels = seq_len(N)))
  adj <- lapply(adj, as.integer)
  names(adj) <- NULL
  deg <- lengths(adj)
  list(N = N, n_m = n_m, n_f = n_f,
       sex = c(rep("M", n_m), rep("F", n_f)),
       deg = deg, adj = adj, edges = e)
}

# ── Core group: female sex workers (FSW) ────────────────────────────────────
# Brewer et al. (2000, PNAS 97:12385) showed that the men-women gap in
# reported numbers of opposite-sex partners is accounted for by female sex
# workers, who are under-sampled by household surveys and have very many
# partners (Colorado Springs: median 103, mean 347 partners per 6 months;
# 23 full-time-equivalent FSW per 100,000 population).
# Construction: men and household women keep their NSFG degrees exactly.
# The men's surplus S = sum(k_men) - sum(k_women) is NOT trimmed; those
# stubs are taken from men's "extra" stubs (a man of degree k contributes
# at most k - 1, chosen in proportion to k - 1, as in trim_extra) and paired
# with n_fsw = round(fsw_frac * n_f) FSW nodes. FSW degrees split S with
# log-normal weights (sdlog from Brewer's median/mean: sqrt(2 log(347/103))).
# FSW degrees are NOT capped at KMAX. Partnerships between the same client
# and FSW are collapsed (repeat client).
# Client concentration q: a man with k partners is chosen for FSW
# partnerships with weight proportional to (k - 1)^q (each man still gives
# at most k - 1 stubs, so men's NSFG degrees are unchanged). q = 1 spreads
# the surplus thinly (~2 FSW contacts per client); larger q puts it on fewer,
# more active men. Calibrated with code/09_core_network_check.R.
FSW_FRAC_DEFAULT <- 0.0005      # Brewer's 23 FTE FSW per 100,000 population, applied to women;
                                # reproduces Brewer's FSW partner numbers (see README)
FSW_SDLOG        <- sqrt(2 * log(347 / 103))
CLIENT_Q_DEFAULT <- 8               # ~4.6% of men 25-44 are clients (floor ~4.3%; code/09)

build_network_core <- function(n_m, n_f, pmf_m, pmf_f, fsw_frac = FSW_FRAC_DEFAULT,
                               sdlog = FSW_SDLOG, q = CLIENT_Q_DEFAULT, seed = 1L) {
  set.seed(seed)
  dm <- sample_pmf(n_m, pmf_m)
  df <- sample_pmf(n_f, pmf_f)
  S  <- sum(dm) - sum(df)
  if (S <= 0) stop("Women report at least as many partners as men; no surplus for FSW.")
  n_fsw <- max(1L, round(fsw_frac * n_f))
  extra <- rep.int(seq_len(n_m), pmax(dm - 1L, 0L))
  if (length(extra) < S) stop("Not enough multi-partner stubs among men for the FSW surplus.")
  # per-stub weight (k - 1)^(q - 1) gives a per-man weight proportional to (k - 1)^q
  # Weighted sampling without replacement via exponential keys
  # (Efraimidis-Spirakis; same law as successive weighted draws, but fast).
  w_stub <- (dm[extra] - 1)^(q - 1)
  keys   <- log(runif(length(extra))) / w_stub
  client_stubs <- extra[order(keys, decreasing = TRUE)[seq_len(S)]]   # one entry per FSW partnership
  dm_house <- dm - tabulate(client_stubs, nbins = n_m)
  w     <- rlnorm(n_fsw, 0, sdlog)
  d_fsw <- as.vector(rmultinom(1, S, w / sum(w)))

  u1 <- rep.int(seq_len(n_m), dm_house)                        # household partnerships
  v1 <- n_m + sample(rep.int(seq_len(n_f), df))
  u2 <- sample(client_stubs)                                   # FSW partnerships
  v2 <- n_m + n_f + rep.int(seq_len(n_fsw), d_fsw)
  e  <- unique(rbind(cbind(u1, v1), cbind(u2, v2)))
  N  <- n_m + n_f + n_fsw
  adj <- split(c(e[, 2], e[, 1]), factor(c(e[, 1], e[, 2]), levels = seq_len(N)))
  adj <- lapply(adj, as.integer); names(adj) <- NULL
  list(N = N, n_m = n_m, n_f = n_f + n_fsw, n_household_f = n_f, n_fsw = n_fsw,
       sex = c(rep("M", n_m), rep("F", n_f + n_fsw)),
       fsw = c(rep(FALSE, n_m + n_f), rep(TRUE, n_fsw)),
       deg = lengths(adj), adj = adj, edges = e, surplus = S)
}

make_network_core <- function(age = "25-44", n_per_sex = 50000L, fsw_frac = FSW_FRAC_DEFAULT,
                              q = CLIENT_Q_DEFAULT, seed = 1L) {
  pm <- nsfg_pmf(NSFG_PCT[[age]]$M, LILJEROS_CUM_ALPHA[["M"]])
  pf <- nsfg_pmf(NSFG_PCT[[age]]$F, LILJEROS_CUM_ALPHA[["F"]])
  net <- build_network_core(n_per_sex, n_per_sex, pm, pf, fsw_frac = fsw_frac, q = q, seed = seed)
  net$name <- sprintf("NSFG %s + FSW core (%.2f%%, q = %g)", age, 100 * fsw_frac, q)
  net
}

# Bridges: men independently bisexual with prob p_bi, kept only if they
# have at least one female partner in the network (k >= 1).
choose_bridges <- function(net, p_bi, seed = 2L) {
  set.seed(seed)
  bi <- runif(net$n_m) < p_bi
  which(bi & net$deg[seq_len(net$n_m)] >= 1L)
}

# Realised degree pmfs (k = 0..max) by sex
realised_pmf <- function(net, sex) {
  k <- net$deg[net$sex == sex]
  tabulate(k + 1L, nbins = max(k) + 1L) / length(k)
}

# Convenience: the networks used in the paper
make_network <- function(name = c("NSFG 25-44", "NSFG 15-44", "NB(1.2,0.7)"),
                         n_per_sex = 5000L, seed = 1L) {
  name <- match.arg(name)
  if (name == "NB(1.2,0.7)") {
    pm <- pf <- nb_pmf(1.2, 0.7)
  } else {
    age <- sub("NSFG ", "", name)
    pm <- nsfg_pmf(NSFG_PCT[[age]]$M, LILJEROS_CUM_ALPHA[["M"]])
    pf <- nsfg_pmf(NSFG_PCT[[age]]$F, LILJEROS_CUM_ALPHA[["F"]])
  }
  net <- build_network(n_per_sex, n_per_sex, pm, pf, seed = seed)
  net$name <- name
  net
}
