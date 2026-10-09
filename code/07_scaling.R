###############################################################
# 07_scaling.R   (theory only; no epidemic simulation)
# Is there a "targeting transition"? Finite-size scaling of the
# degree-targeted importation measure  pi_i ~ k_i^gamma  over n bridges
# whose degrees have a power-law tail  p(k) ~ k^-a  (no degree cap).
#
# pi is a Gibbs measure: gamma = inverse temperature, log k = -energy.
# Predictions for n -> infinity (a > 2):
#   gamma_1 = a - 2 : E[k^(gamma+1)] diverges -> mean degree of the
#                     targeted bridge, mu(gamma), grows with n.
#   gamma_2 = a - 1 : E[k^gamma] diverges -> importation CONDENSES on a
#                     few bridges (sum dominated by its largest terms).
#   Growth exponent of mu, theta(gamma) = d log mu / d log n:
#       0                         gamma < a - 2
#       (gamma - a + 2)/(a - 1)   a - 2 < gamma < a - 1
#       1/(a - 1)                 gamma > a - 1   (mu ~ k_max ~ n^{1/(a-1)})
#   Participation ratio Y2 = sum_i pi_i^2 (order parameter):
#       E[Y2] -> 0 for gamma < a - 1,  -> 1 - (a - 1)/gamma for gamma > a - 1
#       (k^gamma has tail index beta = (a - 1)/gamma; for beta < 1 the
#        normalised weights are Poisson-Dirichlet(beta), E[sum w^2] = 1 - beta).
#   Finite n: crossover of width ~ 1/log n around gamma_2; collapse is
#   tested with x = (gamma - gamma_2) log n.
#
# Depletion: does gamma*(c) collapse when plotted against c/n?
#   C(gamma) = sum_i (1 - exp(-c pi_i)) (1 + A k_i)
#
# Degree distributions of bridges (all k >= 1):
#   "power a=2.8", "power a=3.31", "power a=4": P(K >= k) = k^-(a-1)
#   "NSFG men 25-44": NSFG 2006-08 classes k = 1, 2, 3 and a 4+ tail with
#                     the Liljeros men's exponent (a = 3.31), uncapped.
#   "NSFG men, capped 25": as used in the network model (cap removes the
#                     transition; included for contrast).
#
# All sums use the counts of DISTINCT degree values, so n = 10^6 is fast.
#
# Run from the repo root:   Rscript code/07_scaling.R
# Output: results/07_*.csv   (figures: section "Scaling" in 06_figures.R)
# Runtime: a few minutes.
###############################################################

suppressPackageStartupMessages({ library(dplyr) })
source("code/network.R")                    # NSFG_PCT, LILJEROS_CUM_ALPHA
dir.create("results", showWarnings = FALSE)
set.seed(2029)

# ---- Degree samplers (k >= 1) -------------------------------------------------

r_power <- function(n, a) floor(runif(n)^(-1 / (a - 1)))          # P(K >= k) = k^-(a-1)
NSFG_M  <- NSFG_PCT$`25-44`$M[2:5] / sum(NSFG_PCT$`25-44`$M[2:5]) # k = 1, 2, 3, 4+
A_MEN   <- LILJEROS_CUM_ALPHA[["M"]] + 1                           # pmf exponent 3.31
r_nsfg  <- function(n, cap = Inf) {
  cls <- sample.int(4, n, replace = TRUE, prob = NSFG_M)
  k <- cls
  t <- cls == 4
  k[t] <- floor(4 * runif(sum(t))^(-1 / (A_MEN - 1)))             # tail ~ k^-a from 4
  pmin(k, cap)
}
DISTS <- list(
  `power a=2.8`         = list(a = 2.8,   r = function(n) r_power(n, 2.8)),
  `power a=3.31`        = list(a = 3.31,  r = function(n) r_power(n, 3.31)),
  `power a=4`           = list(a = 4,     r = function(n) r_power(n, 4)),
  `NSFG men 25-44`      = list(a = A_MEN, r = function(n) r_nsfg(n)),
  `NSFG men, capped 25` = list(a = A_MEN, r = function(n) r_nsfg(n, cap = 25))
)

# ---- Observables on distinct values ------------------------------------------
# For a sample of degrees: kv = distinct values, nv = counts.
observables <- function(kv, nv, gamma) {
  lk <- log(kv)
  t(sapply(gamma, function(g) {
    x  <- g * lk; m <- max(x); w <- exp(x - m)          # k^g / max(k^g)
    Z  <- sum(nv * w)
    c(gamma = g,
      mu = sum(nv * w * kv) / Z,                        # mean degree of targeted bridge
      Y2 = sum(nv * w^2) / Z^2,                         # participation ratio
      pmax = max(w) / Z)                                # share of the single top bridge
  }))
}

# ---- Part 1: mu, Y2, pi_max vs gamma for n = 10^2 ... 10^6 --------------------

G1   <- seq(0, 8, by = 0.05)
NS   <- 10^(2:6)
REPS <- c(`100` = 2000, `1000` = 1000, `10000` = 1000, `1e+05` = 500, `1e+06` = 200)

rows <- list()
for (dn in names(DISTS)) {
  D <- DISTS[[dn]]
  for (n in NS) {
    R <- REPS[[format(n, scientific = n >= 1e5)]]
    obs <- replicate(R, {
      k  <- D$r(n); tb <- tabulate(k)
      kv <- which(tb > 0); nv <- tb[kv]
      o  <- observables(kv, nv, G1)
      o[, "mu"] <- o[, "mu"] / (sum(nv * kv) / n)       # relative to mu(0) = mean degree
      o
    }, simplify = "array")                               # gamma x obs x rep
    rows[[length(rows) + 1]] <- data.frame(
      dist = dn, a = D$a, n = n, reps = R, gamma = G1,
      mu_ratio_med = apply(obs[, "mu", ], 1, median),
      log_mu_ratio_mean = apply(log(obs[, "mu", ]), 1, mean),
      Y2_mean = apply(obs[, "Y2", ], 1, mean),
      Y2_q25 = apply(obs[, "Y2", ], 1, quantile, 0.25),
      Y2_q75 = apply(obs[, "Y2", ], 1, quantile, 0.75),
      pmax_mean = apply(obs[, "pmax", ], 1, mean))
    cat(sprintf("%-20s n = %7d done (%d reps)\n", dn, n, R))
  }
}
curves <- bind_rows(rows) %>%
  mutate(gamma1 = a - 2, gamma2 = a - 1,
         Y2_limit = ifelse(gamma > a - 1, 1 - (a - 1) / gamma, 0),
         theta_theory = ifelse(gamma < a - 2, 0, ifelse(gamma < a - 1, (gamma - a + 2) / (a - 1), 1 / (a - 1))))
write.csv(curves, "results/07_scaling_curves.csv", row.names = FALSE)

# Growth exponent theta(gamma): slope of mean log mu(gamma)/mu(0) vs log n (n >= 10^3)
theta <- curves %>% filter(n >= 1e3) %>% group_by(dist, a, gamma) %>%
  summarise(theta_fit = coef(lm(log_mu_ratio_mean ~ log(n)))[2],
            theta_theory = first(theta_theory), .groups = "drop")
write.csv(theta, "results/07_scaling_theta.csv", row.names = FALSE)

cat("\nGrowth exponent theta(gamma) at selected gamma (fit vs theory):\n")
print(as.data.frame(theta %>% filter(gamma %in% c(0.5, 1, 1.5, 2, 2.5, 3, 4, 6)) %>%
        select(dist, gamma, theta_fit, theta_theory)), digits = 3)
cat("\nY2 at n = 10^6 vs the limit 1 - (a-1)/gamma:\n")
print(as.data.frame(curves %>% filter(n == 1e6, gamma %in% c(1, 2, 3, 4, 6, 8)) %>%
        select(dist, gamma, Y2_mean, Y2_limit)), digits = 3)

# ---- Part 2: gamma*(c) vs c/n (depletion) -------------------------------------

G2    <- seq(0, 10, by = 0.05)
NS2   <- 10^(2:5)
XS    <- 10^seq(-3, 1.5, by = 0.5)         # x = c / n, importations per bridge
A_VAL <- c(0.5, 2)                         # secondary cases per partner (A = s(k)/k)
REPS2 <- 50

rows2 <- list()
for (dn in names(DISTS)) {
  D <- DISTS[[dn]]
  for (n in NS2) {
    res_n <- matrix(NA_real_, 0, 6)
    for (rr in seq_len(REPS2)) {
      k  <- D$r(n); tb <- tabulate(k)
      kv <- which(tb > 0); nv <- tb[kv]
      X  <- outer(log(kv), G2); X <- sweep(X, 2, apply(X, 2, max)); W <- exp(X)
      P  <- sweep(W, 2, colSums(nv * W), "/")           # pi of ONE bridge of each distinct degree
      for (Av in A_VAL) for (x in XS) {
        Cg <- colSums(nv * (1 - exp(-(x * n) * P)) * (1 + Av * kv))
        j  <- which.max(Cg)
        flat <- (Cg[j] - Cg[1]) <= 1e-9 * Cg[1]
        res_n <- rbind(res_n, c(n, rr, Av, x, if (flat) NA else G2[j], Cg[j] / Cg[1]))
      }
    }
    colnames(res_n) <- c("n", "rep", "A", "x", "gstar", "gain")
    rows2[[length(rows2) + 1]] <- data.frame(dist = dn, a = D$a, res_n)
    cat(sprintf("gamma*: %-20s n = %6d done\n", dn, n))
  }
}
gstar <- bind_rows(rows2) %>%
  group_by(dist, a, n, A, x) %>%
  summarise(frac_ident = mean(!is.na(gstar)), gstar_med = median(gstar, na.rm = TRUE),
            gstar_q25 = quantile(gstar, 0.25, na.rm = TRUE), gstar_q75 = quantile(gstar, 0.75, na.rm = TRUE),
            gain_med = median(gain), .groups = "drop")
write.csv(gstar, "results/07_gamma_star_scaling.csv", row.names = FALSE)
cat("\nSaved results/07_scaling_curves.csv, 07_scaling_theta.csv, 07_gamma_star_scaling.csv\n")
