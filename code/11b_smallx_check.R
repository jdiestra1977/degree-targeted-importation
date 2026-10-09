###############################################################
# 11b_smallx_check.R   (theory only; ~5-15 min)
# Recomputes the small-x prediction of gamma*_inf in results of
# 11_generalization_crossover.R (Part B) after fixing gamma_star_smallx():
# the asymptotic optimum is the interior local maximum of G closest to
# gamma_2 (G diverges spuriously at gamma_1, where the tail approximation
# fails). Only x <= 1e-2 is recomputed (the formula is a small-x result).
#
# Run from the repo root:   Rscript code/11b_smallx_check.R
# Output: overwrites the gstar_smallx column of results/11_partB_gamma_star.csv
###############################################################
suppressPackageStartupMessages({ library(dplyr) })
source("code/network.R"); source("code/theory.R")

NSFG_M <- NSFG_PCT$`25-44`$M[2:5] / sum(NSFG_PCT$`25-44`$M[2:5])
LAWS <- list(
  `power a=2.8`    = bridge_law(numeric(0), 1, 1, 2.8),
  `power a=3.31`   = bridge_law(numeric(0), 1, 1, 3.31),
  `power a=4`      = bridge_law(numeric(0), 1, 1, 4),
  `NSFG men 25-44` = bridge_law(NSFG_M[1:3], NSFG_M[4], 4, LILJEROS_CUM_ALPHA[["M"]] + 1)
)
B <- read.csv("results/11_partB_gamma_star.csv")
B$gstar_smallx <- NA_real_
for (i in which(B$x <= 1e-2)) {
  B$gstar_smallx[i] <- gamma_star_smallx(LAWS[[B$law[i]]], B$x[i], B$eta[i])
}
write.csv(B, "results/11_partB_gamma_star.csv", row.names = FALSE)
options(width = 150)
print(B %>% filter(x <= 1e-2, x %in% 10^c(-8, -6, -4, -3, -2)) %>%
        transmute(law, eta, x, exact = gstar_inf, smallx = gstar_smallx, leading = gstar_leading,
                  err_smallx = smallx - exact, err_leading = leading - exact), digits = 3)
