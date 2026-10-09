###############################################################
# theory.R
# Analytical results for degree-targeted importation into a bipartite
# (men-women) configuration-model network, matched EXACTLY to the
# discrete-time simulator in sim.R (locally tree-like approximation).
#
# Notation
#   p      daily per-partnership transmission probability, p = 1 - exp(-beta)
#   r      daily recovery probability,                     r = 1 - exp(-rho)
#   D      number of days on which an infectious person can transmit.
#          In the simulator a person who has just become infectious faces
#          recovery before their first transmission day, so
#          P(D = d) = r (1 - r)^d,  d = 0, 1, 2, ...
#   tau    per-partnership transmissibility, tau = E[1 - (1 - p)^D]
#            = (1 - r) p / (p + r - p r)          (tau <= 1 - r)
#   kappa  excess degree <k(k-1)>/<k>, by sex (m = men, f = women)
#   R      two-step (man -> woman -> man) reproduction number^(1/2):
#            R^2 = tau^2 kappa_m kappa_f
#
# Results implemented
#   (1) Mean secondary cases from an imported case in a man of degree k
#         s(k) = tau k (1 + tau kappa_f) / (1 - tau^2 kappa_m kappa_f)
#       Linear in k; valid below threshold (R < 1).
#   (2) Targeting. Importation into bridge i is proportional to
#       w(k_i)^gamma. Mean degree of the bridge that receives an
#       importation: mu(gamma) = sum_i pi_i k_i, pi_i ~ w(k_i)^gamma.
#       Secondary cases per importation = A mu(gamma), A = s(k)/k, so
#       the targeting gain in SECONDARY cases is mu(gamma)/mu(0),
#       independent of transmissibility; in total cases (imported +
#       secondary) it is (1 + A mu(gamma)) / (1 + A mu(0)), which rises
#       with tau towards mu(gamma)/mu(0) as R -> 1.
#   (3) Depletion. With expected (undepleted) importations
#       c = beta_ext |B| sum_t alpha_t, bridge i is ever imported with
#       prob 1 - exp(-c pi_i); expected realised importations and cases
#         I(gamma) = sum_i (1 - exp(-c pi_i))
#         C(gamma) = sum_i (1 - exp(-c pi_i)) (1 + A k_i)
#       Concentrating importations (large |gamma|) wastes them on
#       already-infected bridges, so C(gamma) can peak at an
#       intermediate gamma when c is large.
#       Properties (w increasing in k, bridge degrees not all equal):
#         - I(gamma) is maximised at gamma = 0 (Jensen: 1 - e^{-c pi}
#           is concave in pi);
#         - dC/dgamma at 0 = c e^{-c/|B|} A Cov_B(log w(k), k) > 0, so the
#           burden-maximising gamma* is always > 0;
#         - as c -> 0, C ~ c (1 + A mu(gamma)), increasing: gamma* -> inf;
#         - in general dC/dgamma = c sum_i pi_i e^{-c pi_i} (l_i - lbar)
#           (1 + A k_i), l_i = log w(k_i), lbar = sum_i pi_i l_i.
#   (4) Above threshold: probability that an importation into a man of
#       degree k starts a major outbreak, 1 - h(k, u_f), with
#         h(j, x) = E_D[(1 - tau_D (1 - x))^j],  tau_D = 1 - (1 - p)^D
#         u_m = sum_j pm1_j h(j, u_f),  u_f = sum_j pf1_j h(j, u_m)
#       (pm1, pf1 = excess-degree pmfs). Uses the exact geometric D, so
#       it includes the correlation between a person's partnerships.
#       Probability that forcing produces at least one major outbreak:
#         1 - prod_i [exp(-c pi_i) + (1 - exp(-c pi_i)) h(k_i, u_f)]
###############################################################

# ── Transmissibility ────────────────────────────────────────────────────────

tau_discrete <- function(beta, rho) {
  p <- 1 - exp(-beta); r <- 1 - exp(-rho)
  (1 - r) * p / (p + r - p * r)
}
# Inverse: daily beta giving transmissibility tau (requires tau < 1 - r)
beta_for_tau <- function(tau, rho) {
  r <- 1 - exp(-rho)
  if (any(tau >= 1 - r)) stop("tau must be < 1 - r = ", round(1 - r, 4))
  p <- tau * r / ((1 - r) * (1 - tau))
  -log(1 - p)
}

# ── Degree distributions ────────────────────────────────────────────────────

pmf_mean   <- function(pmf) sum((seq_along(pmf) - 1) * pmf)
pmf_kappa  <- function(pmf) { k <- seq_along(pmf) - 1; sum(k * (k - 1) * pmf) / sum(k * pmf) }
excess_pmf <- function(pmf) {           # pm1_j, j = 0..kmax-1
  k <- seq_along(pmf) - 1
  (k * pmf)[-1] / sum(k * pmf)
}

# Everything the formulas need, from a built network and bridge set
net_summary <- function(net, B) {
  pm <- realised_pmf(net, "M"); pf <- realised_pmf(net, "F")
  list(pm = pm, pf = pf, kappa_m = pmf_kappa(pm), kappa_f = pmf_kappa(pf),
       kB = net$deg[B], nB = length(B))
}

# ── (1) Mean secondary cases ────────────────────────────────────────────────

R_two_step <- function(tau, kappa_m, kappa_f) tau * sqrt(kappa_m * kappa_f)
tau_threshold <- function(kappa_m, kappa_f) 1 / sqrt(kappa_m * kappa_f)

A_factor <- function(tau, kappa_m, kappa_f) {
  R2 <- tau^2 * kappa_m * kappa_f
  ifelse(R2 < 1, tau * (1 + tau * kappa_f) / (1 - R2), NA_real_)
}
mean_secondary <- function(k, tau, kappa_m, kappa_f) A_factor(tau, kappa_m, kappa_f) * k

# Node-level ("quenched") version for an imported case in man i, using
# his actual partners j and their degrees k_j; beyond them, the tree
# average is used:
#   s_i = tau * sum_j [1 + tau (k_j - 1) M],  M = (1 + tau kappa_m) / (1 - R^2)
# M is the expected cluster size (him included) started by a man infected
# along a partnership. Averaging k_j - 1 over partners (mean kappa_f)
# recovers s(k) = A k.
mean_secondary_node <- function(net, i, tau, kappa_m, kappa_f) {
  R2 <- tau^2 * kappa_m * kappa_f
  if (R2 >= 1) return(rep(NA_real_, length(i)))
  M <- (1 + tau * kappa_m) / (1 - R2)
  vapply(i, function(ii) { kj <- net$deg[net$adj[[ii]]]; tau * sum(1 + tau * (kj - 1) * M) },
         numeric(1))
}

# ── (2) Targeting ───────────────────────────────────────────────────────────

target_weights <- function(kB, gamma, w = c("k", "log1p")) {
  w <- match.arg(w)
  wk <- if (w == "k") kB else log1p(kB)
  x  <- gamma * log(wk)                      # stable for large |gamma|
  e  <- exp(x - max(x))
  e / sum(e)
}
mu_gamma <- function(kB, gamma, w = "k") sapply(gamma, function(g) sum(target_weights(kB, g, w) * kB))

gain_secondary <- function(kB, gamma, w = "k") mu_gamma(kB, gamma, w) / mu_gamma(kB, 0, w)
gain_total <- function(kB, gamma, tau, kappa_m, kappa_f, w = "k") {
  A <- A_factor(tau, kappa_m, kappa_f)
  (1 + A * mu_gamma(kB, gamma, w)) / (1 + A * mu_gamma(kB, 0, w))
}

# ── (3) Depletion of the bridge pool under continuous forcing ──────────────

c_expected_importations <- function(beta_ext, nB, alpha) beta_ext * nB * sum(alpha)

forced_expectations <- function(kB, gamma, c, tau, kappa_m, kappa_f, w = "k") {
  A <- A_factor(tau, kappa_m, kappa_f)
  t(sapply(gamma, function(g) {
    q <- 1 - exp(-c * target_weights(kB, g, w))          # P(bridge ever imported)
    c(gamma = g, importations = sum(q), cases = sum(q * (1 + A * kB)),
      secondary = sum(q * A * kB))
  }))
}

# Node-level version: sB = node-level secondary cases of each bridge
# (mean_secondary_node), in the same order as kB.
forced_expectations_node <- function(kB, sB, gamma, c, w = "k") {
  t(sapply(gamma, function(g) {
    q <- 1 - exp(-c * target_weights(kB, g, w))
    c(gamma = g, importations = sum(q), cases = sum(q * (1 + sB)))
  }))
}

# Fast C(gamma) on a whole gamma grid for one bridge set: returns a
# matrix of targeting probabilities (bridges x gamma); reuse it for many
# c and tau.
pi_matrix <- function(kB, gamma, w = "k") {
  lw <- log(if (w == "k") kB else log1p(kB))
  X  <- outer(lw, gamma)
  X  <- sweep(X, 2, apply(X, 2, max))
  E  <- exp(X)
  sweep(E, 2, colSums(E), "/")
}
cases_curve <- function(PI, c, value) {        # value_i = 1 + A k_i  or  1 + s_i
  colSums((1 - exp(-c * PI)) * value)
}

# Large-|B| ensemble version: nB bridges with i.i.d. degrees, pmf b over
# k = 1..K (b[j] = P(k = j)). Approximation: treats every degree class
# as present in proportion b_k (breaks down when nB * b_k < 1).
cases_ensemble <- function(gamma, c, nB, b, A, w = "k") {
  k  <- seq_along(b)
  wk <- if (w == "k") k else log1p(k)
  sapply(gamma, function(g) {
    Mg <- sum(b * wk^g)
    nB * sum(b * (1 - exp(-(c / nB) * wk^g / Mg)) * (1 + A * k))
  })
}

# Optimum of a curve on a grid (grid maximum; refine the grid if needed)
argmax_grid <- function(x, y) x[which.max(y)]

# ── Multi-type theory for the FSW-core network ──────────────────────────────
# Partnerships are household (man - non-FSW woman) or commercial (client -
# FSW). A man has h household and m commercial partners (correlated), a
# household woman d partners, an FSW d_F clients. Arrival types (along a
# partnership):  MW = at a household woman, WM = at a man via a household
# partnership, MF = at an FSW, FM = at a client via a commercial partnership.
# M[a, b] = expected number of onward partnerships of type b after arriving
# by type a (size-biased by the arrival partnership):
#   M[MW, WM] = <d(d-1)>/<d>               (household women)
#   M[WM, MW] = <h(h-1)>/<h>,  M[WM, MF] = <h m>/<h>       (men)
#   M[MF, FM] = <d_F(d_F-1)>/<d_F>          (FSW)
#   M[FM, MW] = <m h>/<m>,     M[FM, MF] = <m(m-1)>/<m>    (men)
# Each onward partnership transmits with prob tau, so K = tau M and
#   R = tau rho(M),  tau_c = 1 / rho(M)   (rho = spectral radius)
#   z = (I - tau M)^{-1} 1   expected cluster size started by one arrival
#   s(h, m) = tau (h z_MW + m z_MF)   secondary cases from an imported man
# With no FSW this reduces to s(k) = A k (one household type).
multitype_M <- function(net) {
  fsw <- net$fsw
  men <- seq_len(net$n_m)
  h <- vapply(men, function(i) sum(!fsw[net$adj[[i]]]), numeric(1))
  m <- net$deg[men] - h
  d  <- net$deg[net$sex == "F" & !fsw]
  dF <- net$deg[fsw]
  sb <- function(x) sum(x * (x - 1)) / sum(x)
  M <- matrix(0, 4, 4, dimnames = list(c("MW", "WM", "MF", "FM"), c("MW", "WM", "MF", "FM")))
  M["MW", "WM"] <- sb(d)
  M["WM", "MW"] <- sb(h);              M["WM", "MF"] <- sum(h * m) / sum(h)
  M["MF", "FM"] <- if (length(dF)) sb(dF) else 0
  if (sum(m) > 0) { M["FM", "MW"] <- sum(m * h) / sum(m); M["FM", "MF"] <- sb(m) }
  rho <- max(Mod(eigen(M, only.values = TRUE)$values))
  list(M = M, rho = rho, tau_c = 1 / rho, h = h, m = m)
}
multitype_z <- function(MT, tau) {
  if (tau * MT$rho >= 1) return(rep(NA_real_, 4))
  setNames(as.vector(solve(diag(4) - tau * MT$M, rep(1, 4))), rownames(MT$M))
}
multitype_secondary <- function(MT, i, tau) {       # i = man indices
  z <- multitype_z(MT, tau)
  tau * (MT$h[i] * z[["MW"]] + MT$m[i] * z[["MF"]])
}

# ── Targeting transition: large-|B| limit with a power-law degree tail ──────
# Bridge degree law: point masses pm_low on k = 1..k0-1 and, from k0, the
# discrete tail P(K >= k) = w_tail (k / k0)^-(a-1) (the law of
# floor(k0 * U^{-1/(a-1)})). Summed exactly up to kcut; beyond kcut the
# continuous density w_tail (a-1) k0^(a-1) k^-a is integrated.
# Pure power law: pm_low = numeric(0), k0 = 1, w_tail = 1.
bridge_law <- function(pm_low, w_tail, k0, a, kcut = 1e6) {
  k  <- k0:kcut
  pt <- w_tail * ((k / k0)^-(a - 1) - ((k + 1) / k0)^-(a - 1))
  list(k = c(seq_along(pm_low), k), p = c(pm_low, pt), a = a, k0 = k0,
       w_tail = w_tail, kcut = kcut, cdens = w_tail * (a - 1) * k0^(a - 1))
}

# M(gamma) = E[k^gamma]; finite only for gamma < a - 1
law_M <- function(L, g) {
  if (g >= L$a - 1) return(Inf)
  K <- L$kcut + 1
  sum(L$p * L$k^g) + L$cdens * K^(g - L$a + 1) / (L$a - 1 - g)
}

# Limit of C(gamma)/|B| as |B| -> infinity at fixed x = c/|B|:
#   E[(1 - exp(-x k^gamma / M(gamma))) (1 + A k)],  gamma < a - 1.
# (For gamma >= a - 1 the importation condenses on O(1) bridges and
#  C/|B| -> 0.)
# Response exponent eta: an importation into a bridge with weight w causes
# A w^eta secondary cases (degree case: w = k, eta = 1). Requires
# eta < a - 1 (finite mean response). Default eta = 1 = original model.
law_C_per_bridge <- function(L, g, x, A, eta = 1) {
  M <- law_M(L, g)
  if (!is.finite(M)) return(0)
  body <- sum(L$p * (-expm1(-x * L$k^g / M)) * (1 + A * L$k^eta))
  K <- L$kcut + 1
  tail <- integrate(function(s) {
    k <- exp(s); L$cdens * k^(1 - L$a) * (-expm1(-x * k^g / M)) * (1 + A * k^eta)
  }, lower = log(K), upper = log(K) + 80, rel.tol = 1e-8)$value
  body + tail
}

# gamma*_inf(x): maximiser of the limit curve over [0, a - 1)
gamma_star_limit <- function(L, x, A, n_grid = 60L, eta = 1) {
  gmax <- L$a - 1 - 1e-4
  gr <- seq(0, gmax, length.out = n_grid)
  v  <- vapply(gr, function(g) law_C_per_bridge(L, g, x, A, eta), numeric(1))
  j  <- which.max(v)
  lo <- gr[max(1, j - 1)]; hi <- gr[min(n_grid, j + 1)]
  if (hi - lo < 1e-12) return(gr[j])
  optimize(function(g) law_C_per_bridge(L, g, x, A, eta), c(lo, hi), maximum = TRUE, tol = 1e-5)$maximum
}

# Small-x asymptotics of gamma*_inf (derivation in the technical report).
# For x -> 0 the burden per bridge is dominated by bridges near the
# saturation scale w_s = (E[w^g]/x)^(1/g); with the tail p(w) ~ w^-a,
#   C/n - x  ~  const * (x / E[w^g])^nu * Gamma(1 - nu),
#   nu = (a - 1 - eta)/g in (0, 1),  i.e. g in (gamma_1, gamma_2),
# where the constant does not depend on g or A. Hence
#   gamma*(x) ~ argmax_g  G(g) = -nu (ln(1/x) + ln E[w^g]) + ln Gamma(1 - nu),
# with E[w^g] exact (law_M). Leading order: gamma_2 - gamma* ~ (a-1)/ln(1/x)
# (validated numerically; next-order terms are not captured by G alone).
gamma_star_smallx <- function(L, x, eta = 1) {
  # The tail approximation behind G holds only for nu bounded away from 1
  # (gamma away from gamma_1): as nu -> 1, ln Gamma(1 - nu) -> +inf, a
  # spurious divergence because the integral is then dominated by small w.
  # The asymptotic optimum is therefore the interior local maximum of G
  # closest to gamma_2 (NA if there is none).
  g1 <- max(L$a - 1 - eta, 0) + 1e-6; g2 <- L$a - 1 - 1e-6
  G <- function(g) { nu <- (L$a - 1 - eta) / g; -nu * (log(1 / x) + log(law_M(L, g))) + lgamma(1 - nu) }
  gr <- seq(g1, g2, length.out = 300)
  v  <- vapply(gr, G, numeric(1))
  j  <- which(diff(sign(diff(v))) == -2) + 1               # interior local maxima
  if (!length(j)) return(NA_real_)
  j  <- max(j)
  optimize(G, c(gr[j - 1], gr[j + 1]), maximum = TRUE, tol = 1e-8)$maximum
}

# Large-x asymptote gamma* ~ K / x, where K > 0 is the root of
#   E[(log k - <log k>) (1 + A k) k^-K] = 0   (weights w = k)
# Derivation: with u_i = |B| pi_i and d_i = log k_i - <log k>, the
# condition dC/dgamma = 0 reads sum_i e^{-x u_i} u_i (d_i - dbar_gamma)
# (1 + A k_i) = 0; for gamma = K/x -> 0, u_i -> 1 + gamma d_i and
# e^{-x u_i} -> e^{-x} e^{-K d_i}, leaving E[d (1 + A k) e^{-K d}] = 0.
# The root is unique: the left side is A Cov(log k, k) > 0 at K = 0 and
# strictly decreasing in K. (Truncated at kcut; the neglected tail is
# O(kcut^{2-a-K} log kcut).)
asymptote_K <- function(L, A, eta = 1) {
  lk <- log(L$k); lbar <- sum(L$p * lk) / sum(L$p)
  d  <- lk - lbar
  g  <- function(K) sum(L$p * d * (1 + A * L$k^eta) * exp(-K * d))
  uniroot(g, c(1e-8, 50), tol = 1e-10)$root
}

# ── (4) Major outbreaks (above threshold) ───────────────────────────────────

# h(j, x) for vector j; geometric D summed to negligible tail
h_fun <- function(j, x, p, r, tol = 1e-12) {
  dmax <- ceiling(log(tol) / log(1 - r))
  d    <- 0:dmax
  wD   <- r * (1 - r)^d
  tauD <- 1 - (1 - p)^d
  sapply(j, function(jj) sum(wD * (1 - tauD * (1 - x))^jj))
}

major_outbreak_fixed_point <- function(pm, pf, p, r, iters = 2000L, tol = 1e-12) {
  pm1 <- excess_pmf(pm); pf1 <- excess_pmf(pf)
  jm <- seq_along(pm1) - 1; jf <- seq_along(pf1) - 1
  um <- 0; uf <- 0
  for (it in seq_len(iters)) {
    um_new <- sum(pm1 * h_fun(jm, uf, p, r))
    uf_new <- sum(pf1 * h_fun(jf, um_new, p, r))
    if (abs(um_new - um) + abs(uf_new - uf) < tol) { um <- um_new; uf <- uf_new; break }
    um <- um_new; uf <- uf_new
  }
  list(u_m = um, u_f = uf)
}
p_major_from_man <- function(k, uf, p, r) 1 - h_fun(k, uf, p, r)

# Approximate expected size of a major outbreak (number infected), used
# to set the major/minor cut-off in simulations. A person of degree k
# escapes if none of their k partnerships brings infection; each does so
# with prob tau (1 - u) (u = prob that the partner's side does not reach
# the giant outbreak). Ignores duration correlation; fine for a cut-off.
major_outbreak_size <- function(pm, pf, n_m, n_f, tau, um, uf) {
  km <- seq_along(pm) - 1; kf <- seq_along(pf) - 1
  frac_m <- 1 - sum(pm * (1 - tau * (1 - uf))^km)
  frac_f <- 1 - sum(pf * (1 - tau * (1 - um))^kf)
  n_m * frac_m + n_f * frac_f
}

p_any_major_forced <- function(kB, gamma, c, uf, p, r, w = "k") {
  hk <- h_fun(kB, uf, p, r)
  sapply(gamma, function(g) {
    q <- 1 - exp(-c * target_weights(kB, g, w))
    1 - prod(1 - q + q * hk)
  })
}
