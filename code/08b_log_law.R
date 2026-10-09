###############################################################
# 08b_log_law.R   (theory only; a few minutes)
# Firms up the sparse-spillover limit of gamma*: in the n -> infinity
# limit (exact degree law, no sampling), how does the gap
#   delta(x) = gamma_2 - gamma*_inf(x),   gamma_2 = a - 1,
# close as x = c/n -> 0? 08 suggested delta ~ c_a / ln(1/x) over
# x = 1e-4 ... 1e-2. Here x goes down to 1e-8 and two forms are fitted:
#   (i)  delta = c_a / ln(1/x)
#   (ii) delta = c_a / (ln(1/x) + b)          (allows a constant shift)
# If (i) holds, delta * ln(1/x) is flat in x.
#
# Run from the repo root:   Rscript code/08b_log_law.R
# Output: results/08b_log_law.csv, results/08b_log_law_fits.csv
###############################################################

suppressPackageStartupMessages({ library(dplyr) })
source("code/network.R"); source("code/theory.R")
dir.create("results", showWarnings = FALSE)

A_VAL  <- 2
XS     <- 10^seq(-8, -1, by = 0.5)
NSFG_M <- NSFG_PCT$`25-44`$M[2:5] / sum(NSFG_PCT$`25-44`$M[2:5])
LAWS <- list(
  `power a=2.8`    = bridge_law(numeric(0), 1, 1, 2.8),
  `power a=3.31`   = bridge_law(numeric(0), 1, 1, 3.31),
  `power a=4`      = bridge_law(numeric(0), 1, 1, 4),
  `NSFG men 25-44` = bridge_law(NSFG_M[1:3], NSFG_M[4], 4, LILJEROS_CUM_ALPHA[["M"]] + 1)
)

res <- bind_rows(lapply(names(LAWS), function(dn) {
  L <- LAWS[[dn]]
  g <- sapply(XS, function(x) gamma_star_limit(L, x, A_VAL, n_grid = 120L))
  cat(sprintf("%s done\n", dn))
  data.frame(dist = dn, a = L$a, gamma2 = L$a - 1, x = XS, gstar_inf = g,
             delta = L$a - 1 - g, L1x = log(1 / XS))
})) %>% mutate(delta_times_ln = delta * L1x)
write.csv(res, "results/08b_log_law.csv", row.names = FALSE)

fits <- bind_rows(lapply(split(res, res$dist), function(d) {
  d <- filter(d, x <= 1e-2)
  f1 <- nls(delta ~ c1 / L1x, data = d, start = list(c1 = 2))
  f2 <- tryCatch(nls(delta ~ c2 / (L1x + b), data = d, start = list(c2 = 2, b = 0)),
                 error = function(e) NULL)
  data.frame(dist = d$dist[1], a = d$a[1], gamma2 = d$gamma2[1],
             c_a_simple = coef(f1)[["c1"]],
             rmse_simple = sqrt(mean(residuals(f1)^2)),
             c_a_shifted = if (is.null(f2)) NA else coef(f2)[["c2"]],
             b_shifted   = if (is.null(f2)) NA else coef(f2)[["b"]],
             rmse_shifted = if (is.null(f2)) NA else sqrt(mean(residuals(f2)^2)))
}))
write.csv(fits, "results/08b_log_law_fits.csv", row.names = FALSE)

options(width = 160)
cat("\ndelta * ln(1/x) (flat if delta ~ c_a / ln(1/x)):\n")
print(as.data.frame(res %>% select(dist, x, gstar_inf, delta, delta_times_ln)), digits = 4)
cat("\nFits over x <= 1e-2:\n")
print(fits, digits = 4)
