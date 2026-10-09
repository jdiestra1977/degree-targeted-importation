###############################################################
# sim.R
# Discrete-time stochastic SIR on the bipartite network, built to match
# theory.R exactly. Each day:
#   (1) [forced runs only] each susceptible bridge i is imported with
#       prob 1 - exp(-beta_ext * alpha_t * |B| * pi_i);
#   (2) each infectious person infects each susceptible partner with
#       prob p, independently (equivalently, a susceptible with m
#       infectious partners escapes with prob (1 - p)^m = exp(-beta m));
#       (1) and (2) both use the state at the start of the day;
#   (3) every infectious person, including those infected today,
#       recovers with prob r.
# Runs continue until nobody is infectious, so final sizes are not
# truncated by a time window. The latent (E) stage of an SEIR model
# delays infections but does not change final size, so it is omitted.
# States: 0 = S, 1 = I, 2 = R.
###############################################################

# Single importation into `seed` (no forcing). Returns secondary cases;
# if they reach `stop_at` the run stops and returns stop_at (used to
# classify major outbreaks without simulating them to the end).
sim_cluster <- function(net, seed, p, r, stop_at = Inf) {
  state <- integer(net$N)
  state[seed] <- 1L
  I <- seed[runif(length(seed)) >= r]          # may recover before transmitting
  state[setdiff(seed, I)] <- 2L
  secondary <- 0L
  while (length(I)) {
    nb  <- unlist(net$adj[I], use.names = FALSE)
    nb  <- nb[state[nb] == 0L]
    new <- unique(nb[runif(length(nb)) < p])
    state[new] <- 1L
    secondary <- secondary + length(new)
    if (secondary >= stop_at) return(stop_at)
    I   <- c(I, new)
    rec <- runif(length(I)) < r
    state[I[rec]] <- 2L
    I   <- I[!rec]
  }
  secondary
}

# Continuous importation driven by alpha (daily, in [0, 1]) into bridges B
# with targeting probabilities pi (same order as B).
# Returns realised importations, secondary (network) cases, total cases,
# and major = 1 if total cases reached `stop_at` (run stopped there).
sim_forced <- function(net, B, pi, alpha, beta_ext, p, r, stop_at = Inf) {
  state <- integer(net$N)
  nB <- length(B)
  I <- integer(0)
  n_imp <- 0L; n_sec <- 0L
  t <- 0L; Tf <- length(alpha)
  while (t < Tf || length(I)) {
    t <- t + 1L
    # (1) importation
    imp <- integer(0)
    if (t <= Tf && alpha[t] > 0) {
      s <- state[B] == 0L
      if (any(s)) {
        prob <- 1 - exp(-beta_ext * alpha[t] * nB * pi[s])
        imp  <- B[s][runif(sum(s)) < prob]
      }
    }
    # (2) transmission from people infectious at the start of the day
    new <- integer(0)
    if (length(I)) {
      nb  <- unlist(net$adj[I], use.names = FALSE)
      nb  <- nb[state[nb] == 0L]
      new <- unique(nb[runif(length(nb)) < p])
    }
    new <- setdiff(new, imp)                   # imported today counts as importation
    state[c(imp, new)] <- 1L
    n_imp <- n_imp + length(imp); n_sec <- n_sec + length(new)
    if (n_imp + n_sec >= stop_at)
      return(c(importations = n_imp, secondary = n_sec, total = n_imp + n_sec, major = 1))
    I <- c(I, imp, new)
    # (3) recovery
    if (length(I)) {
      rec <- runif(length(I)) < r
      state[I[rec]] <- 2L
      I <- I[!rec]
    }
  }
  c(importations = n_imp, secondary = n_sec, total = n_imp + n_sec, major = 0)
}

# Weekly incidence -> daily alpha in [0, 1] (each week's value repeated 7 days)
load_alpha <- function(path) {
  x <- read.csv(path)
  x <- x[order(x$t), ]
  rep(x$incidence / max(x$incidence), each = 7L)
}
