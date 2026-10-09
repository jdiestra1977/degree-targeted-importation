###############################################################
# 11_generalization_crossover.R   (theory only; no epidemic simulation)
# Paper 1 (PRE) additions. Bridge "importance" w with tail p(w) ~ w^-a,
# targeting pi_i ~ w_i^gamma, response per importation A w^eta.
#
#   Part A  Generalized growth exponent of the targeted response
#             mu_eta(gamma) = sum_i pi_i w_i^eta,  theta_eta = d log mu_eta / d log n
#           Prediction: 0 for gamma < gamma_1 = a - 1 - eta,
#             (gamma - gamma_1)/(a - 1) between, eta/(a - 1) above gamma_2 = a - 1.
#   Part B  gamma*(x) for eta in {0.5, 1, 2}: exact n -> infinity limit vs
#           (i) the small-x asymptotic  argmax G(gamma) (derived, no fit),
#           (ii) its leading order gamma_2 - (a - 1)/ln(1/x),
#           (iii) the large-x asymptote K/x.
#   Part C  Cut-off crossover: degree laws capped at K_max in {25, 50, 100,
#           400}. Prediction: the transition is rounded once the expected
#           largest degree reaches K_max, at n_c = (K_max/k0)^(a-1) / P_tail;
#           capped / uncapped ratios of mu and Y2 collapse against n / n_c.
#
# Run from the repo root:   Rscript code/11_generalization_crossover.R
# Output: results/11_*.csv     Runtime: ~10-25 min (Part B dominates).
###############################################################

suppressPackageStartupMessages({ library(dplyr) })
source("code/network.R"); source("code/theory.R")
dir.create("results", showWarnings = FALSE)
set.seed(2033)

NSFG_M <- NSFG_PCT$`25-44`$M[2:5] / sum(NSFG_PCT$`25-44`$M[2:5])   # bridges: k = 1, 2, 3, 4+
A_MEN  <- LILJEROS_CUM_ALPHA[["M"]] + 1                            # 3.31
r_power <- function(n, a, cap = Inf) pmin(floor(runif(n)^(-1 / (a - 1))), cap)
r_nsfg  <- function(n, cap = Inf) {
  cls <- sample.int(4, n, replace = TRUE, prob = NSFG_M); k <- cls; t <- cls == 4
  k[t] <- floor(4 * runif(sum(t))^(-1 / (A_MEN - 1))); pmin(k, cap)
}

# Observables over distinct values (kv with counts nv)
obs_eta <- function(kv, nv, gamma, eta) {
  lk <- log(kv)
  t(sapply(gamma, function(g) {
    x <- g * lk; w <- exp(x - max(x)); Z <- sum(nv * w)
    c(mu = sum(nv * w * kv^eta) / Z, Y2 = sum(nv * w^2) / Z^2)
  }))
}

# ---- Part A: theta_eta -----------------------------------------------------------
GA   <- seq(0, 6, by = 0.05)
NS   <- 10^(2:6)
REPS <- c(2000, 1000, 500, 200, 100)
ETAS <- c(0.5, 1, 1.5, 2)
A_laws <- list(`power a=3.31` = function(n) r_power(n, 3.31), `NSFG men 25-44` = function(n) r_nsfg(n))
rowsA <- list()
for (ln in names(A_laws)) for (ni in seq_along(NS)) {
  n <- NS[ni]
  arr <- replicate(REPS[ni], {
    k <- A_laws[[ln]](n); tb <- tabulate(k); kv <- which(tb > 0); nv <- tb[kv]
    sapply(ETAS, function(e) { o <- obs_eta(kv, nv, GA, e); log(o[, "mu"] / (sum(nv * kv^e) / n)) })
  }, simplify = "array")                                   # gamma x eta x rep
  for (ei in seq_along(ETAS))
    rowsA[[length(rowsA) + 1]] <- data.frame(law = ln, a = 3.31, eta = ETAS[ei], n = n,
                                             gamma = GA, log_mu_ratio = rowMeans(arr[, ei, ]))
  cat(sprintf("Part A: %s n = %d\n", ln, n))
}
A_df <- bind_rows(rowsA)
thetaA <- A_df %>% filter(n >= 1e3) %>% group_by(law, a, eta, gamma) %>%
  summarise(theta_fit = coef(lm(log_mu_ratio ~ log(n)))[2], .groups = "drop") %>%
  mutate(g1 = a - 1 - eta, g2 = a - 1,
         theta_theory = ifelse(gamma < g1, 0, ifelse(gamma < g2, (gamma - g1) / (a - 1), eta / (a - 1))))
write.csv(A_df, "results/11_partA_curves.csv", row.names = FALSE)
write.csv(thetaA, "results/11_partA_theta.csv", row.names = FALSE)
print(as.data.frame(thetaA %>% filter(gamma %in% c(0.5, 1, 1.5, 2, 3, 5)) %>%
        select(law, eta, gamma, theta_fit, theta_theory)), digits = 3)

# ---- Part B: gamma*(x) with response exponent eta --------------------------------
A_VAL <- 2
XS    <- 10^seq(-8, 1, by = 0.5)
LAWS <- list(
  `power a=2.8`    = bridge_law(numeric(0), 1, 1, 2.8),
  `power a=3.31`   = bridge_law(numeric(0), 1, 1, 3.31),
  `power a=4`      = bridge_law(numeric(0), 1, 1, 4),
  `NSFG men 25-44` = bridge_law(NSFG_M[1:3], NSFG_M[4], 4, A_MEN)
)
rowsB <- list()
for (ln in names(LAWS)) for (eta in c(0.5, 1, 2)) {
  L <- LAWS[[ln]]
  if (eta >= L$a - 1) next                                    # mean response must be finite
  K <- asymptote_K(L, A_VAL, eta)
  for (x in XS) {
    rowsB[[length(rowsB) + 1]] <- data.frame(
      law = ln, a = L$a, eta = eta, x = x, gamma2 = L$a - 1, gamma1 = L$a - 1 - eta, K = K,
      gstar_inf   = gamma_star_limit(L, x, A_VAL, n_grid = 100L, eta = eta),
      gstar_smallx = gamma_star_smallx(L, x, eta),
      gstar_leading = L$a - 1 - (L$a - 1) / log(1 / x),
      gstar_K = K / x)
  }
  cat(sprintf("Part B: %s eta = %g\n", ln, eta))
}
B_df <- bind_rows(rowsB)
write.csv(B_df, "results/11_partB_gamma_star.csv", row.names = FALSE)
print(as.data.frame(B_df %>% filter(x %in% XS[c(1, 5, 9, 13, 17, 19)]) %>%
        select(law, eta, x, gstar_inf, gstar_smallx, gstar_leading, gstar_K)), digits = 3)

# ---- Part C: cut-off crossover ------------------------------------------------------
CAPS <- c(25, 50, 100, 400)
GC   <- c(1.5, 3, 4)                                          # below gamma_2 = 2.31 and above
REPSC <- c(2000, 1000, 500, 200, 100)
C_laws <- list(`power a=3.31` = list(r = function(n, cap) r_power(n, 3.31, cap), k0 = 1, Pt = 1),
               `NSFG men 25-44` = list(r = function(n, cap) r_nsfg(n, cap), k0 = 4, Pt = NSFG_M[4]))
rowsC <- list()
for (ln in names(C_laws)) for (cap in c(CAPS, Inf)) for (ni in seq_along(NS)) {
  n <- NS[ni]; Lr <- C_laws[[ln]]
  arr <- replicate(REPSC[ni], {
    k <- Lr$r(n, cap); tb <- tabulate(k); kv <- which(tb > 0); nv <- tb[kv]
    o <- obs_eta(kv, nv, GC, 1); cbind(mu = o[, "mu"] / (sum(nv * kv) / n), Y2 = o[, "Y2"])
  }, simplify = "array")                                      # gamma x obs x rep
  rowsC[[length(rowsC) + 1]] <- data.frame(law = ln, cap = cap, n = n, gamma = GC,
    n_c = if (is.finite(cap)) (cap / Lr$k0)^(A_MEN - 1) / Lr$Pt else Inf,
    mu_ratio = apply(arr[, "mu", ], 1, mean), Y2 = apply(arr[, "Y2", ], 1, mean))
  cat(sprintf("Part C: %s cap = %s n = %d\n", ln, cap, n))
}
C_df <- bind_rows(rowsC)
unc  <- C_df %>% filter(!is.finite(cap)) %>% select(law, n, gamma, mu_unc = mu_ratio, Y2_unc = Y2)
C_df <- C_df %>% filter(is.finite(cap)) %>% left_join(unc, by = c("law", "n", "gamma")) %>%
  mutate(n_over_nc = n / n_c, mu_rel = mu_ratio / mu_unc, Y2_rel = Y2 / Y2_unc)
write.csv(C_df, "results/11_partC_crossover.csv", row.names = FALSE)
print(as.data.frame(C_df %>% filter(gamma == 3) %>% select(law, cap, n, n_c, n_over_nc, mu_rel, Y2_rel)), digits = 3)
cat("\nSaved results/11_*.csv\n")
