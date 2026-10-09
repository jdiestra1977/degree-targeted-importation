###############################################################
# 08_gamma_star_limit.R   (theory only; no epidemic simulation)
# Firms up the gamma*(c) results of 07_scaling.R.
#
# x = c / n = expected importations per bridge, A = 2, weights w = k.
#   Part A  Finite n: gamma* for n = 10^2 ... 10^6, x = 10^-4 ... 10^1.5,
#           with a fine gamma grid (0.005 below 1, 0.01 up to 10).
#   Part B  n -> infinity at fixed x: gamma*_inf(x) maximises
#             E[(1 - exp(-x k^g / E[k^g])) (1 + A k)]   over g < gamma_2,
#           computed from the exact degree law (no sampling).
#           Question: is gamma* pinned at gamma_2 = a - 1, or does it only
#           approach gamma_2 as x -> 0?
#   Part C  Large-x asymptote gamma* ~ K / x, K the root of
#             E[(log k - <log k>) (1 + A k) k^-K] = 0.
#
# Run from the repo root:   Rscript code/08_gamma_star_limit.R
# Output: results/08_*.csv  (figure: section "Scaling F" in 06_figures.R)
# Runtime: ~10-20 min (Part A dominates).
###############################################################

suppressPackageStartupMessages({ library(dplyr) })
source("code/network.R"); source("code/theory.R")
dir.create("results", showWarnings = FALSE)
set.seed(2030)

A_VAL <- 2
XS    <- 10^seq(-4, 1.5, by = 0.5)

# Degree laws (bridges, k >= 1) and matching samplers
NSFG_M <- NSFG_PCT$`25-44`$M[2:5] / sum(NSFG_PCT$`25-44`$M[2:5])
A_MEN  <- LILJEROS_CUM_ALPHA[["M"]] + 1
LAWS <- list(
  `power a=2.8`    = bridge_law(numeric(0), 1, 1, 2.8),
  `power a=3.31`   = bridge_law(numeric(0), 1, 1, 3.31),
  `power a=4`      = bridge_law(numeric(0), 1, 1, 4),
  `NSFG men 25-44` = bridge_law(NSFG_M[1:3], NSFG_M[4], 4, A_MEN)
)
r_law <- function(n, L) {
  if (L$k0 == 1) return(floor(runif(n)^(-1 / (L$a - 1))))
  cls <- sample.int(4, n, replace = TRUE, prob = c(L$p[1:3], L$w_tail))
  k <- cls; t <- cls == 4
  k[t] <- floor(L$k0 * runif(sum(t))^(-1 / (L$a - 1)))
  k
}

# ---- Part C: large-x asymptote ------------------------------------------------
asym <- data.frame(dist = names(LAWS), a = sapply(LAWS, `[[`, "a"),
                   K = sapply(LAWS, asymptote_K, A = A_VAL))
print(asym, digits = 4)
write.csv(asym, "results/08_asymptote_K.csv", row.names = FALSE)

# ---- Part B: n -> infinity limit curve ----------------------------------------
lim <- bind_rows(lapply(names(LAWS), function(dn) {
  L <- LAWS[[dn]]
  data.frame(dist = dn, a = L$a, x = XS,
             gstar_inf = sapply(XS, function(x) gamma_star_limit(L, x, A_VAL)))
}))
lim$gamma2 <- lim$a - 1
lim$gap_to_gamma2 <- lim$gamma2 - lim$gstar_inf
print(lim, digits = 4)
write.csv(lim, "results/08_gamma_star_limit.csv", row.names = FALSE)

# ---- Part A: finite n, fine grid ------------------------------------------------
G    <- c(seq(0, 1, by = 0.005), seq(1.01, 10, by = 0.01))
NS   <- 10^(2:6)
REPS <- c(60, 60, 60, 40, 20)
rows <- list()
for (dn in names(LAWS)) {
  L <- LAWS[[dn]]
  for (ni in seq_along(NS)) {
    n <- NS[ni]
    for (rr in seq_len(REPS[ni])) {
      k  <- r_law(n, L); tb <- tabulate(k)
      kv <- which(tb > 0); nv <- tb[kv]
      X  <- outer(log(kv), G); X <- sweep(X, 2, apply(X, 2, max)); W <- exp(X)
      P  <- sweep(W, 2, colSums(nv * W), "/")
      v  <- nv * (1 + A_VAL * kv)
      for (x in XS) {
        Cg <- colSums((1 - exp(-(x * n) * P)) * v)
        j  <- which.max(Cg)
        flat <- (Cg[j] - Cg[1]) <= 1e-9 * Cg[1]
        rows[[length(rows) + 1]] <- data.frame(dist = dn, a = L$a, n = n, rep = rr, x = x,
                                              gstar = if (flat) NA else G[j])
      }
    }
    cat(sprintf("Part A: %-16s n = %7d done\n", dn, n))
  }
}
fin <- bind_rows(rows) %>% group_by(dist, a, n, x) %>%
  summarise(frac_ident = mean(!is.na(gstar)), gstar_med = median(gstar, na.rm = TRUE),
            gstar_q25 = quantile(gstar, 0.25, na.rm = TRUE),
            gstar_q75 = quantile(gstar, 0.75, na.rm = TRUE), .groups = "drop") %>%
  left_join(select(lim, dist, x, gstar_inf), by = c("dist", "x")) %>%
  left_join(select(asym, dist, K), by = "dist") %>%
  mutate(gstar_asym = K / x, gamma2 = a - 1)
write.csv(fin, "results/08_gamma_star_finite.csv", row.names = FALSE)

cat("\ngamma* (median) by n vs the n -> infinity limit and the K/x asymptote:\n")
print(as.data.frame(fin %>% select(dist, x, n, gstar_med, gstar_inf, gstar_asym, gamma2) %>%
        filter(x %in% XS[c(1, 3, 5, 7, 9, 11)])), digits = 3)
cat("\nSaved results/08_asymptote_K.csv, 08_gamma_star_limit.csv, 08_gamma_star_finite.csv\n")
