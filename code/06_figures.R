###############################################################
# 06_figures.R   (interactive; RStudio)
# Publication figures, built from the saved results of scripts 03-11
# plus fast analytical curves (no epidemic simulation).
#
# Step through with Ctrl+Alt+T (run current section); each figure is
# printed to the Plots pane. With SAVE_PDFS <- TRUE each is also saved
# to figures/paper/.
#
# Figures are saved at their printed size (PRE: one column = 3.4 in,
# two columns = 7.0 in), so text appears at 8-9 pt in the paper. They
# carry no titles or subtitles; that information is in the captions.
#
#   fig2_targeting_gain      main Fig. 1 (one column, panels stacked)
#   fig_transition           main Fig. 2 (one column; survey law, Y2 and theta)
#   figS_transition          SM (all laws, Y2 and theta)
#   fig3_depletion_gamma_star main Fig. 3 (two columns)
#   figS_gen_theta_eta       main Fig. 4 (one column)
#   figS_gen_crossover       main Fig. 5 (one column)
#   fig1, fig4, fig5, figS1, figS_scaling_*, figS_gen_gamma_star_eta,
#   figS_core_calibration    Supplemental Material
###############################################################

# ---- 0. Setup ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(scales)
})
# Run from the repository root: open degree-targeted-importation.Rproj in RStudio
# (or setwd() to the folder that contains code/ and results/).
if (!file.exists("code/network.R")) stop("Working directory must be the repository root")
source("code/network.R"); source("code/theory.R")

SAVE_PDFS <- TRUE
FIG_DIR   <- "figures/paper"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

RHO  <- 1 / 14
P_BI <- 0.02
COL  <- 3.4    # PRE column width (in)
PAGE <- 7.0    # PRE two-column width (in)

theme_paper <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = "grey60"),
        strip.text = element_text(size = 8),
        plot.tag = element_text(size = 9),
        legend.position = "bottom",
        legend.key.size = unit(0.35, "cm"),
        legend.key.width = unit(0.8, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        plot.margin = margin(3, 5, 3, 3))
show <- function(p, name, w = PAGE, h = 4) {
  print(p)
  if (SAVE_PDFS) ggsave(file.path(FIG_DIR, name), p, width = w, height = h)
  invisible(p)
}
tau_scale <- function(...) scale_colour_viridis_d(name = expression(tau), end = 0.85, ...)

# Labels used in the paper for the degree laws (strips and legends)
law_lab <- function(x) {
  map <- c("power a=2.8" = "power law, a = 2.8", "power a=3.31" = "power law, a = 3.31",
           "power a=4" = "power law, a = 4", "NSFG men 25-44" = "survey, a = 3.31",
           "NSFG men, capped 25" = "survey, capped at 25",
           "NSFG 25-44" = "survey, ages 25-44", "NSFG 15-44" = "survey, ages 15-44")
  x <- as.character(x); out <- unname(map[x]); out[is.na(out)] <- x[is.na(out)]; out
}
lab_law <- as_labeller(law_lab)

# ---- 1. Load saved results ---------------------------------------------------

A3 <- read.csv("results/03_cluster_ensemble.csv")
B3 <- read.csv("results/03_forced_ensemble.csv")
G4 <- read.csv("results/04_gamma_star_summary.csv")
A5 <- read.csv("results/05_pmajor_by_degree.csv")
B5 <- read.csv("results/05_forced_major.csv")

# ---- Fig 1 (SM). Secondary cases per seeding event ---------------------------
# Points: ensemble simulation (100 networks, +-2 network-clustered SE),
# degrees present in >= 5 networks, k <= 10. Lines: s(k) = A k.
# Panel (b): probability of no secondary case, h(k, 0).

f1_dat <- A3 %>% filter(enough_networks, k <= 10)
f1a <- ggplot(f1_dat, aes(k, sim_mean, colour = factor(tau))) +
  geom_line(aes(y = theory_ens), linewidth = 0.6) +
  geom_pointrange(aes(ymin = sim_mean - 2 * sim_se_net, ymax = sim_mean + 2 * sim_se_net),
                  size = 0.15, linewidth = 0.35) +
  facet_wrap(~ network, scales = "free_y", labeller = lab_law) +
  scale_x_continuous(breaks = 1:10) + tau_scale() +
  labs(x = "Degree k of the seeded node", y = "Mean secondary cases", tag = "(a)") +
  theme_paper
f1b <- ggplot(f1_dat, aes(k, sim_p0, colour = factor(tau))) +
  geom_line(aes(y = theory_p0), linewidth = 0.6) +
  geom_point(size = 1.1) +
  facet_wrap(~ network, labeller = lab_law) +
  scale_x_continuous(breaks = 1:10) + scale_y_continuous(limits = c(0, 1)) + tau_scale() +
  labs(x = "Degree k of the seeded node", y = "P(no secondary case)", tag = "(b)") +
  guides(colour = "none") +
  theme_paper
show((f1a / f1b) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
     "fig1_secondary_cases.pdf", w = PAGE, h = 5.6)

# ---- Fig 2 (main Fig. 1). Targeting gain per seeding event -------------------
# (a) gain in secondary cases mu(gamma)/mu(0) (independent of transmissibility).
#     Red: entry nodes of the survey 25-44 networks (median over 100 networks;
#     degrees capped at 25; dotted: median bound k_max/mu(0)). Greys: n entry
#     nodes with the uncapped survey law, median over 200 samples.
# (b) gain in total cases (1 + A mu(gamma)) / (1 + A mu(0)) on the survey
#     networks for three transmissibilities; dashed: secondary-case gain.

GG <- seq(-3, 8, by = 0.05)
f2_net <- list(); f2_tot <- list()
for (s in 1:100) {
  net <- make_network("NSFG 25-44", seed = 5000L + s); B <- choose_bridges(net, P_BI, seed = 6000L + s)
  S <- net_summary(net, B)
  if (length(unique(S$kB)) < 2) next
  mu <- mu_gamma(S$kB, GG, "k"); mu0 <- mu_gamma(S$kB, 0, "k")
  f2_net[[length(f2_net) + 1]] <- data.frame(net = s, gamma = GG, gain = mu / mu0, bound = max(S$kB) / mu0)
  for (tau in c(0.3, 0.6, 0.9)) {
    A <- A_factor(tau, S$kappa_m, S$kappa_f)
    f2_tot[[length(f2_tot) + 1]] <- data.frame(net = s, gamma = GG, tau = tau,
                                               gain = (1 + A * mu) / (1 + A * mu0), sec = mu / mu0)
  }
}
f2_net <- bind_rows(f2_net) %>% group_by(gamma) %>%
  summarise(gain = median(gain), bound = median(bound), .groups = "drop") %>%
  mutate(pool = "survey network, n = 90")
f2_tot <- bind_rows(f2_tot) %>% group_by(tau, gamma) %>%
  summarise(gain = median(gain), sec = median(sec), .groups = "drop")

# Uncapped survey law for men, n entry nodes (same sampler as 07_scaling.R)
NSFG_M4 <- NSFG_PCT$`25-44`$M[2:5] / sum(NSFG_PCT$`25-44`$M[2:5])
r_nsfg_uncapped <- function(n) {
  cls <- sample.int(4, n, replace = TRUE, prob = NSFG_M4); k <- cls; t <- cls == 4
  k[t] <- floor(4 * runif(sum(t))^(-1 / LILJEROS_CUM_ALPHA[["M"]])); k
}
set.seed(606)
f2_unc <- bind_rows(lapply(c(100, 1000, 10000), function(n) {
  m <- replicate(200, { k <- r_nsfg_uncapped(n); mu_gamma(k, GG, "k") / mean(k) })
  data.frame(gamma = GG, gain = apply(m, 1, median),
             pool = sprintf("uncapped, n = %s", format(n, big.mark = ",")))
}))
pool_cols <- c("survey network, n = 90" = "#d7301f", "uncapped, n = 100" = "grey70",
               "uncapped, n = 1,000" = "grey45", "uncapped, n = 10,000" = "grey15")

f2a <- ggplot(bind_rows(f2_net, f2_unc), aes(gamma, gain, colour = pool)) +
  geom_line(linewidth = 0.7) +
  geom_hline(data = f2_net[1, ], aes(yintercept = bound), colour = "#d7301f", linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = LILJEROS_CUM_ALPHA[["M"]] - 1, linetype = "dotted", colour = "grey40") +
  scale_y_log10() + scale_colour_manual(values = pool_cols, name = NULL) +
  guides(colour = guide_legend(ncol = 2)) +
  labs(x = expression("Targeting exponent " * gamma), y = expression(mu(gamma) / mu(0)), tag = "(a)") +
  theme_paper + theme(legend.key.width = unit(0.4, "cm"))   # keeps the 2-column legend inside 3.4 in
f2b <- ggplot(f2_tot, aes(gamma, gain, colour = factor(tau))) +
  geom_line(aes(y = sec), colour = "black", linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  tau_scale() +
  labs(x = expression("Targeting exponent " * gamma), y = "Gain in total cases", tag = "(b)") +
  theme_paper
show(f2a / f2b, "fig2_targeting_gain.pdf", w = COL, h = 5.6)

# ---- Fig 3 (main). Repeated forcing with depletion ---------------------------
# (a) expected total cases vs gamma, survey 25-44 ensemble (50 networks),
#     simulation (points +-2 SE) and node-level theory (lines), tau = 0.6.
# (b) gamma* vs x = c / n. Red: survey networks (median and IQR over 100
#     networks, tau = 0.6). Black / grey: n -> infinity limit for uncapped
#     laws (code/08, 08b) and K/x (dashed), A = 2. Dotted: gamma_2 = a - 1.
# (c) burden at gamma* and at gamma = 30, relative to uniform forcing.

f3a <- ggplot(filter(B3, tau == 0.6), aes(gamma, sim_cases, colour = factor(c))) +
  geom_line(aes(y = th_cases_node), linewidth = 0.6) +
  geom_pointrange(aes(ymin = sim_cases - 2 * se_cases, ymax = sim_cases + 2 * se_cases),
                  size = 0.15, linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_y_log10() +
  scale_colour_viridis_d(name = "c", option = "C", end = 0.85) +
  labs(x = expression("Targeting exponent " * gamma), y = "Expected total cases", tag = "(a)") +
  theme_paper

# survey networks (04): convert c to x = c / n per network
G4n <- read.csv("results/04_gamma_star_by_network.csv") %>%
  filter(type == "NSFG 25-44", level == 0.6) %>% mutate(x = c / nB)
g4 <- G4n %>% mutate(xb = 10^(round(log10(x) * 10) / 10)) %>% group_by(xb) %>%
  summarise(gstar_med = median(gamma_star, na.rm = TRUE),
            gstar_q25 = quantile(gamma_star, 0.25, na.rm = TRUE),
            gstar_q75 = quantile(gamma_star, 0.75, na.rm = TRUE),
            ident = mean(!is.na(gamma_star)),
            gain_opt = median(gain_opt), gain_hub = median(gain_hub), .groups = "drop")
# n -> infinity limits for uncapped laws (08b extends 08 to smaller x)
lim_file <- if (file.exists("results/08b_log_law.csv")) "results/08b_log_law.csv" else "results/08_gamma_star_limit.csv"
LIM <- bind_rows(read.csv(lim_file), read.csv("results/08_gamma_star_limit.csv")) %>%
  distinct(dist, x, .keep_all = TRUE) %>% select(dist, a, x, gstar_inf)
Kx  <- read.csv("results/08_asymptote_K.csv") %>% filter(dist == "NSFG men 25-44") %>%
  tidyr::crossing(x = 10^seq(-0.5, 1.5, by = 0.05)) %>% mutate(g = K / x)
lim_cols <- c("power a=2.8" = "grey75", "power a=3.31" = "grey50", "power a=4" = "grey25",
              "NSFG men 25-44" = "black")
# leading order of the small-x asymptotics, gamma_2 - (a - 1)/ln(1/x), x <= 1e-2
LO <- distinct(LIM, dist, a) %>% tidyr::crossing(x = 10^seq(-8, -2, by = 0.1)) %>%
  mutate(g = (a - 1) - (a - 1) / log(1 / x))
ltys <- c("limit" = "solid", "leading order" = "dotdash", "K/x" = "dashed")
f3b <- ggplot() +
  geom_hline(data = distinct(LIM, dist, a), aes(yintercept = a - 1, colour = dist),
             linetype = "dotted", show.legend = FALSE) +
  geom_line(data = LIM, aes(x, gstar_inf, colour = dist, linetype = "limit"), linewidth = 0.6) +
  geom_line(data = LO, aes(x, g, colour = dist, linetype = "leading order"), linewidth = 0.5) +
  geom_line(data = Kx, aes(x, g, linetype = "K/x"), colour = "black") +
  geom_ribbon(data = filter(g4, ident >= 0.5), aes(xb, ymin = gstar_q25, ymax = gstar_q75),
              fill = "#d7301f", alpha = 0.2) +
  geom_line(data = filter(g4, ident >= 0.5), aes(xb, gstar_med), colour = "#d7301f", linewidth = 0.8) +
  scale_colour_manual(values = lim_cols, labels = law_lab, name = NULL) +
  scale_linetype_manual(values = ltys, breaks = names(ltys), name = NULL) +
  scale_x_log10() + coord_cartesian(ylim = c(0, 4)) +
  guides(colour = guide_legend(ncol = 2, order = 1), linetype = guide_legend(nrow = 1, order = 2)) +
  labs(x = "x = c/n, seeding events per entry node", y = expression("Optimal exponent " * gamma^"*"),
       tag = "(b)") +
  theme_paper + theme(legend.box = "vertical")
f3c <- g4 %>% select(xb, gain_opt, gain_hub) %>%
  pivot_longer(c(gain_opt, gain_hub), names_to = "what", values_to = "gain") %>%
  ggplot(aes(xb, gain, linetype = what)) +
  geom_line(linewidth = 0.7, colour = "#d7301f") + geom_hline(yintercept = 1, linetype = "dotted") +
  scale_x_log10() + scale_y_log10() +
  scale_linetype_manual(values = c(gain_opt = "solid", gain_hub = "dashed"),
                        breaks = c("gain_opt", "gain_hub"),
                        labels = c(gain_opt = expression(gamma == gamma^"*"), gain_hub = expression(gamma == 30)),
                        name = NULL) +
  labs(x = "x = c/n, seeding events per entry node", y = "Burden relative to uniform", tag = "(c)") +
  theme_paper
show(f3a / (f3b | f3c), "fig3_depletion_gamma_star.pdf", w = PAGE, h = 5.8)

# ---- Fig 4 (SM). Major outbreaks above threshold -----------------------------
# NB(1.2,0.7) networks (supercritical). (a) P(major | one seeding event into a
# class-1 node of degree k). (b) P(at least one major outbreak) under repeated
# forcing vs gamma, tau = 0.6 (points: simulation +-2 SE; lines: theory).

f4a <- ggplot(A5, aes(k, sim, colour = factor(tau))) +
  geom_line(aes(y = theory), linewidth = 0.6) +
  geom_pointrange(aes(ymin = sim - 2 * se, ymax = sim + 2 * se), size = 0.15, linewidth = 0.35) +
  scale_x_continuous(breaks = 1:10) + scale_y_continuous(limits = c(0, 1)) + tau_scale() +
  labs(x = "Degree k of the seeded node", y = "P(major outbreak)", tag = "(a)") +
  theme_paper
f4b <- ggplot(filter(B5, tau == 0.6), aes(gamma, sim, colour = factor(c))) +
  geom_line(aes(y = theory), linewidth = 0.6) +
  geom_pointrange(aes(ymin = pmax(0, sim - 2 * se), ymax = pmin(1, sim + 2 * se)), size = 0.15, linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_viridis_d(name = "c", option = "C", end = 0.85) +
  labs(x = expression("Targeting exponent " * gamma), y = "P(at least one major outbreak)", tag = "(b)") +
  theme_paper
show(f4a | f4b, "fig4_major_outbreaks.pdf", w = PAGE, h = 3)

# ---- Fig S1 (SM). Network calibration ----------------------------------------
# (a) past-year partner distribution, model ensemble (20 networks, pooled)
#     vs NSFG 2006-2008 targets, by sex. (b) kappa_1 (men) and kappa_2 (women)
#     across the 100 networks of script 03.

cal <- list()
for (age in c("25-44", "15-44")) {
  dm <- integer(0); df <- integer(0)
  for (s in 1:20) {
    net <- make_network(paste("NSFG", age), seed = 1000L + s)
    dm <- c(dm, net$deg[net$sex == "M"]); df <- c(df, net$deg[net$sex == "F"])
  }
  cats <- function(k) factor(ifelse(k >= 4, "4+", as.character(k)), levels = c("0", "1", "2", "3", "4+"))
  for (sx in c("M", "F")) {
    k <- if (sx == "M") dm else df
    tgt <- NSFG_PCT[[age]][[sx]]
    cal[[length(cal) + 1]] <- rbind(
      data.frame(age = age, sex = sx, source = "model network", cat = levels(cats(0)),
                 pct = 100 * as.numeric(table(cats(k))) / length(k)),
      data.frame(age = age, sex = sx, source = "survey", cat = levels(cats(0)),
                 pct = 100 * tgt / sum(tgt)))
  }
}
fs1a <- bind_rows(cal) %>%
  mutate(sex = recode(sex, M = "Men", F = "Women"), age = paste("Ages", age)) %>%
  ggplot(aes(cat, pct, fill = source)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  facet_grid(sex ~ age) +
  scale_fill_manual(values = c("model network" = "grey40", "survey" = "#d7301f"), name = NULL) +
  labs(x = "Opposite-sex partners in the past 12 months", y = "%", tag = "(a)") +
  theme_paper
kap <- bind_rows(read.csv("results/03_kappa_by_network_NSFG_25-44.csv"),
                 read.csv("results/03_kappa_by_network_NSFG_15-44.csv")) %>%
  pivot_longer(c(kappa_m, kappa_f), names_to = "sex", values_to = "kappa") %>%
  mutate(sex = recode(sex, kappa_m = "Men", kappa_f = "Women"),
         network = sub("NSFG ", "Ages ", network))
fs1b <- ggplot(kap, aes(network, kappa, fill = sex)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6, position = position_dodge(0.7)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_fill_manual(values = c(Men = "#2c7fb8", Women = "#d7301f"), name = NULL) +
  labs(x = NULL, y = expression("Excess degree " * kappa), tag = "(b)") +
  theme_paper
show(fs1a + fs1b + plot_layout(widths = c(2, 1)), "figS1_network_calibration.pdf", w = PAGE, h = 3.6)

###############################################################################
# SCALING / TARGETING TRANSITION  (needs results of code/07_scaling.R)
# pi_i ~ k_i^gamma over n entry nodes with tail p(k) ~ k^-a, no degree cap.
# Predicted transitions as n -> infinity: gamma_1 = a - 2 (mean degree of
# the targeted node diverges), gamma_2 = a - 1 (forcing condenses).
###############################################################################

# ---- Scaling 0. Load ---------------------------------------------------------

if (file.exists("results/07_scaling_curves.csv")) {
  SC <- read.csv("results/07_scaling_curves.csv")
  TH <- read.csv("results/07_scaling_theta.csv")
  GS <- read.csv("results/07_gamma_star_scaling.csv")
  dist_levels <- c("power a=2.8", "power a=3.31", "power a=4", "NSFG men 25-44", "NSFG men, capped 25")
  SC$dist <- factor(SC$dist, dist_levels); TH$dist <- factor(TH$dist, dist_levels)
  GS$dist <- factor(GS$dist, dist_levels)
  n_scale <- scale_colour_viridis_d(name = "n", option = "D", end = 0.9,
                                    labels = function(x) format(as.numeric(x), big.mark = ",", scientific = FALSE))
  marks <- distinct(SC, dist, a) %>% mutate(g1 = a - 2, g2 = a - 1)
} else message("Run code/07_scaling.R first")

# ---- Scaling A. Order parameter Y2 (condensation) ----------------------------
# Participation ratio Y2 = sum pi_i^2 vs gamma for growing n. Black dashed:
# limit 1 - (a - 1)/gamma above gamma_2 (0 below). Vertical lines: gamma_1
# (dotted), gamma_2 (dashed).

sA <- ggplot(SC, aes(gamma, Y2_mean, colour = factor(n))) +
  geom_line(linewidth = 0.5) +
  geom_line(data = filter(SC, dist != "NSFG men, capped 25"), aes(y = Y2_limit),
            colour = "black", linetype = "dashed", linewidth = 0.4) +
  geom_vline(data = marks, aes(xintercept = g1), linetype = "dotted", colour = "grey40") +
  geom_vline(data = marks, aes(xintercept = g2), linetype = "dashed", colour = "grey40") +
  facet_wrap(~ dist, nrow = 1, labeller = lab_law) + n_scale +
  labs(x = expression("Targeting exponent " * gamma), y = expression(Y[2])) +
  theme_paper
show(sA, "figS_scaling_Y2.pdf", w = PAGE, h = 2.5)

# ---- Scaling B. Collapse of Y2 near gamma_2 ----------------------------------
# Same data vs (gamma - gamma_2) log n; curves for different n collapse near 0
# if the width of the transition scales as 1/log n.

sB <- SC %>% mutate(x = (gamma - (a - 1)) * log(n)) %>% filter(abs(x) <= 15) %>%
  ggplot(aes(x, Y2_mean, colour = factor(n))) +
  geom_line(linewidth = 0.5) + geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ dist, nrow = 1, labeller = lab_law) + n_scale +
  labs(x = expression((gamma - gamma[2]) * ln ~ n), y = expression(Y[2])) +
  theme_paper
show(sB, "figS_scaling_Y2_collapse.pdf", w = PAGE, h = 2.5)

# ---- Scaling C. Growth exponent of the targeted mean degree ------------------
# theta(gamma) = d ln mu / d ln n, fitted over n = 10^3..10^6 (points), vs the
# piecewise-linear prediction (line).

sC <- ggplot(TH, aes(gamma)) +
  geom_line(data = filter(TH, dist != "NSFG men, capped 25"), aes(y = theta_theory),
            colour = "black", linewidth = 0.5) +
  geom_point(aes(y = theta_fit), colour = "#d7301f", size = 0.6) +
  geom_vline(data = marks, aes(xintercept = g1), linetype = "dotted", colour = "grey40") +
  geom_vline(data = marks, aes(xintercept = g2), linetype = "dashed", colour = "grey40") +
  facet_wrap(~ dist, nrow = 1, labeller = lab_law) +
  labs(x = expression("Targeting exponent " * gamma), y = expression(theta(gamma))) +
  theme_paper
show(sC, "figS_scaling_theta.pdf", w = PAGE, h = 2.3)

# ---- Fig. S (SM). Targeting transition for all laws: Y2 (a), theta (b) ------
show((sA + labs(tag = "(a)") + theme(legend.position = "right")) /
       (sC + labs(tag = "(b)")),
     "figS_transition.pdf", w = PAGE, h = 4.4)

# ---- Fig. 2 (main). Targeting transition, survey law only (one column) -------
# (a) Y2 vs gamma for n = 10^2..10^6 entry nodes with the uncapped survey law
#     (a = 3.31); black dashed: n -> infinity limit. (b) theta(gamma) fitted
#     over n = 10^3..10^6 (points) vs the prediction (line). Vertical lines:
#     gamma_1 = a - 2 (dotted), gamma_2 = a - 1 (dashed).
SURV <- "NSFG men 25-44"
mk1  <- filter(marks, dist == SURV)
t2a <- ggplot(filter(SC, dist == SURV), aes(gamma, Y2_mean, colour = factor(n))) +
  geom_line(linewidth = 0.6) +
  geom_line(aes(y = Y2_limit), colour = "black", linetype = "dashed", linewidth = 0.5) +
  geom_vline(data = mk1, aes(xintercept = g1), linetype = "dotted", colour = "grey40") +
  geom_vline(data = mk1, aes(xintercept = g2), linetype = "dashed", colour = "grey40") +
  n_scale + guides(colour = guide_legend(nrow = 2)) +
  labs(x = expression("Targeting exponent " * gamma), y = expression(Y[2]), tag = "(a)") +
  theme_paper
t2b <- ggplot(filter(TH, dist == SURV), aes(gamma)) +
  geom_line(aes(y = theta_theory), colour = "black", linewidth = 0.6) +
  geom_point(aes(y = theta_fit), colour = "#d7301f", size = 0.7) +
  geom_vline(data = mk1, aes(xintercept = g1), linetype = "dotted", colour = "grey40") +
  geom_vline(data = mk1, aes(xintercept = g2), linetype = "dashed", colour = "grey40") +
  labs(x = expression("Targeting exponent " * gamma), y = expression(theta(gamma)), tag = "(b)") +
  theme_paper
show(t2a / t2b, "fig_transition.pdf", w = COL, h = 5)

# ---- Scaling D. Targeted mean degree mu(gamma)/mu(0) vs gamma, by n ---------

sD <- ggplot(SC, aes(gamma, mu_ratio_med, colour = factor(n))) +
  geom_line(linewidth = 0.5) +
  geom_vline(data = marks, aes(xintercept = g1), linetype = "dotted", colour = "grey40") +
  geom_vline(data = marks, aes(xintercept = g2), linetype = "dashed", colour = "grey40") +
  scale_y_log10() + facet_wrap(~ dist, nrow = 1, labeller = lab_law) + n_scale +
  labs(x = expression("Targeting exponent " * gamma), y = expression(mu(gamma) / mu(0))) +
  theme_paper
show(sD, "figS_scaling_mu.pdf", w = PAGE, h = 2.5)

# ---- Scaling E. gamma*(c) vs seeding events per entry node c/n --------------
# Collapse of the curves for different n means gamma* depends on c only
# through c/n. Solid: A = 2; dashed: A = 0.5.

sE <- ggplot(filter(GS, frac_ident >= 0.5), aes(x, gstar_med, colour = factor(n), linetype = factor(A))) +
  geom_line(linewidth = 0.5) +
  scale_x_log10() + facet_wrap(~ dist, nrow = 1, labeller = lab_law) + n_scale +
  scale_linetype_manual(values = c(`0.5` = "dashed", `2` = "solid"), name = "A") +
  labs(x = "x = c/n, seeding events per entry node", y = expression(gamma^"*")) +
  theme_paper
show(sE, "figS_scaling_gamma_star.pdf", w = PAGE, h = 2.6)

# ---- Scaling F. gamma*(x): finite n, n -> infinity limit, asymptote -----------
# Needs code/08_gamma_star_limit.R. x = c/n, A = 2. Colours: finite n
# (median). Black solid: n -> infinity limit. Black dashed: K/x. Grey dashed:
# gamma_2 = a - 1.

if (file.exists("results/08_gamma_star_finite.csv")) {
  F8 <- read.csv("results/08_gamma_star_finite.csv")
  L8 <- read.csv("results/08_gamma_star_limit.csv")
  K8 <- read.csv("results/08_asymptote_K.csv")
  lv8 <- c("power a=2.8", "power a=3.31", "power a=4", "NSFG men 25-44")
  F8$dist <- factor(F8$dist, lv8); L8$dist <- factor(L8$dist, lv8); K8$dist <- factor(K8$dist, lv8)
  asym_curve <- K8 %>% tidyr::crossing(x = 10^seq(-1, 1.5, by = 0.05)) %>% mutate(g = K / x)
  sF <- ggplot(filter(F8, frac_ident >= 0.5), aes(x, gstar_med, colour = factor(n))) +
    geom_line(linewidth = 0.5) +
    geom_line(data = L8, aes(x, gstar_inf), colour = "black", linewidth = 0.7, inherit.aes = FALSE) +
    geom_line(data = asym_curve, aes(x, g), colour = "black", linetype = "dashed", inherit.aes = FALSE) +
    geom_hline(data = distinct(L8, dist, gamma2), aes(yintercept = gamma2),
               linetype = "dashed", colour = "grey45") +
    scale_x_log10() + coord_cartesian(ylim = c(0, 6)) +
    facet_wrap(~ dist, nrow = 1, labeller = lab_law) + n_scale +
    labs(x = "x = c/n, seeding events per entry node", y = expression(gamma^"*")) +
    theme_paper
  show(sF, "figS_scaling_gamma_star_limit.pdf", w = PAGE, h = 2.6)
} else message("Run code/08_gamma_star_limit.R first")

###############################################################################
# CORE-GROUP NETWORK (sex-worker core; needs code/09 and code/10b)
# Baseline: 0.05% of women are sex workers, client concentration q = 8,
# 50,000 per sex. "Large outbreak" = reaching 1,000 cases (1% of the
# population); runs were stopped there, so case counts are capped.
###############################################################################

# ---- Fig 5 (SM). Core network at realistic transmissibility ------------------
# (a) P(large outbreak) after one seeding event into a client vs another
#     entry node, by tau (mean over 10 networks, +-2 SE across networks).
# (b) the same under repeated forcing vs gamma, by tau; columns c = 10, 100.

if (file.exists("results/10b_targeting_raw.csv")) {
  S10 <- read.csv("results/10b_single_importation_raw.csv") %>%
    group_by(tau, seed_type) %>%
    summarise(p = mean(frac_reach), se = sd(frac_reach) / sqrt(n()), .groups = "drop")
  T10 <- read.csv("results/10b_targeting_raw.csv") %>%
    group_by(tau, c, gamma) %>%
    summarise(p = mean(p_reach), se = sd(p_reach) / sqrt(n()), .groups = "drop")

  f5a <- ggplot(S10, aes(tau, p, colour = seed_type)) +
    geom_line(linewidth = 0.7) +
    geom_pointrange(aes(ymin = pmax(0, p - 2 * se), ymax = pmin(1, p + 2 * se)), size = 0.2) +
    scale_colour_manual(values = c("client" = "#d7301f", "non-client bridge" = "grey40"),
                        labels = c("client" = "client of the core", "non-client bridge" = "other entry node"),
                        name = NULL) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = expression("Transmissibility " * tau),
         y = "P(outbreak reaches 1% of population)", tag = "(a)") +
    theme_paper
  f5b <- ggplot(T10, aes(gamma, p, colour = factor(tau))) +
    geom_line(linewidth = 0.6) +
    geom_pointrange(aes(ymin = pmax(0, p - 2 * se), ymax = pmin(1, p + 2 * se)), size = 0.15, linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_wrap(~ paste("c =", c)) +
    scale_y_continuous(limits = c(0, 1)) + tau_scale() +
    labs(x = expression("Targeting exponent " * gamma),
         y = "P(outbreak reaches 1% of population)", tag = "(b)") +
    theme_paper
  show(f5a + f5b + plot_layout(widths = c(1, 1.8)), "fig5_core_network.pdf", w = PAGE, h = 3.2)
} else message("Run code/10b_core_realistic_tau.R first")

# ---- Fig S (SM). Core calibration (code/09) -----------------------------------
# (a) implied partners per sex worker per year vs sex-worker prevalence
#     (median over networks), with Brewer et al.'s measured median (206) and
#     mean (694) per year. (b) share of men who are clients vs client
#     concentration q, with the GSS range (0.6-2%) shaded.

if (file.exists("results/09_core_network_check_summary.csv")) {
  C9 <- read.csv("results/09_core_network_check_summary.csv") %>% filter(fsw_frac > 0)
  fca <- C9 %>% filter(q == 8) %>%
    select(fsw_frac, n_per_sex, fsw_deg_mean, fsw_deg_median) %>%
    pivot_longer(c(fsw_deg_mean, fsw_deg_median), names_to = "stat", values_to = "deg") %>%
    mutate(stat = recode(stat, fsw_deg_mean = "mean", fsw_deg_median = "median")) %>%
    ggplot(aes(100 * fsw_frac, deg, colour = stat, shape = factor(n_per_sex))) +
    geom_hline(data = data.frame(y = c(206, 694), ref = c("observed median", "observed mean")),
               aes(yintercept = y, linetype = ref), colour = "black", inherit.aes = FALSE) +
    geom_line(aes(group = interaction(stat, n_per_sex)), linewidth = 0.5) + geom_point(size = 1.5) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c(mean = "#2c7fb8", median = "#d7301f"), name = "Model") +
    scale_linetype_manual(values = c("observed median" = "dashed", "observed mean" = "dotted"), name = NULL) +
    scale_shape_discrete(name = "People per sex") +
    guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
    labs(x = "Sex workers per 100 women", y = "Partners per sex worker per year", tag = "(a)") +
    theme_paper
  fcb <- C9 %>% filter(fsw_frac == 0.0005) %>%
    ggplot(aes(q, pct_men_clients, shape = factor(n_per_sex))) +
    geom_rect(data = data.frame(lo = 0.6, hi = 2, src = "observed range"),
              aes(xmin = 0.8, xmax = 20, ymin = lo, ymax = hi, fill = src), alpha = 0.2, inherit.aes = FALSE) +
    geom_line(aes(group = n_per_sex), linewidth = 0.5) + geom_point(size = 1.5) +
    scale_x_log10(breaks = c(1, 2, 4, 8, 16)) + scale_y_continuous(limits = c(0, 12)) +
    scale_fill_manual(values = c("observed range" = "#2c7fb8"), name = NULL) +
    scale_shape_discrete(name = "People per sex") +
    guides(fill = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
    labs(x = "Client concentration q", y = "% of men who are clients", tag = "(b)") +
    theme_paper
  show(fca + fcb, "figS_core_calibration.pdf", w = PAGE, h = 3.4)
}

###############################################################################
# GENERALIZATION AND CUT-OFF CROSSOVER (needs code/11_generalization_crossover.R)
###############################################################################

# ---- Gen A (main Fig. 4). Growth exponent theta_eta(gamma) -------------------
# Points: fitted over n = 10^3..10^6; lines: piecewise prediction with
# gamma_1 = a - 1 - eta (dotted, colour of eta) and gamma_2 = a - 1 (dashed).
if (file.exists("results/11_partA_theta.csv")) {
  TA <- read.csv("results/11_partA_theta.csv") %>%
    mutate(law = factor(law, levels = c("power a=3.31", "NSFG men 25-44")))
  gA <- ggplot(TA, aes(gamma, colour = factor(eta))) +
    geom_line(aes(y = theta_theory), linewidth = 0.6) +
    geom_point(aes(y = theta_fit), size = 0.5) +
    geom_vline(data = distinct(TA, law, eta, g1), aes(xintercept = g1, colour = factor(eta)),
               linetype = "dotted", show.legend = FALSE) +
    geom_vline(xintercept = 3.31 - 1, linetype = "dashed", colour = "grey40") +
    facet_wrap(~ law, ncol = 1, labeller = lab_law) +
    scale_colour_viridis_d(name = expression(eta), end = 0.85) +
    labs(x = expression("Targeting exponent " * gamma), y = expression(theta[eta](gamma))) +
    theme_paper
  show(gA, "figS_gen_theta_eta.pdf", w = COL, h = 4.4)
}

# ---- Gen B (SM). gamma*(x): exact limit vs small-x and large-x asymptotes -----
if (file.exists("results/11_partB_gamma_star.csv")) {
  GB <- read.csv("results/11_partB_gamma_star.csv")
  # Each asymptote is drawn only in its own regime: leading order for x <= 1e-2,
  # K/x for x >= 0.3. The a = 2.8, eta = 2 panel is undefined (eta >= a - 1).
  GBl <- GB %>% select(law, eta, x, gamma2, `exact limit` = gstar_inf,
                       `maximum of G` = gstar_smallx, `leading order` = gstar_leading, `K/x` = gstar_K) %>%
    pivot_longer(c(`exact limit`, `maximum of G`, `leading order`, `K/x`), names_to = "curve", values_to = "g") %>%
    filter(!(curve == "leading order" & x > 1e-2), !(curve == "K/x" & x < 0.3)) %>%
    mutate(curve = factor(curve, levels = c("exact limit", "maximum of G", "leading order", "K/x")),
           eta_lab = paste0("eta == ", eta), law_lab = paste0('"', law_lab(law), '"'))
  undefined <- data.frame(eta_lab = "eta == 2", law_lab = paste0('"', law_lab("power a=2.8"), '"'),
                          x = 1e-4, g = 1.6)
  gB <- ggplot(GBl, aes(x, g, colour = curve, linetype = curve)) +
    geom_line(linewidth = 0.6) +
    geom_hline(aes(yintercept = gamma2), colour = "grey50", linetype = "dotted") +
    geom_text(data = undefined, aes(x, g), label = "eta >= a - 1", parse = TRUE,
              inherit.aes = FALSE, size = 2.8, colour = "grey40") +
    facet_grid(eta_lab ~ law_lab, labeller = label_parsed) +
    scale_x_log10() + coord_cartesian(ylim = c(0, 3.2)) +
    scale_colour_manual(values = c("exact limit" = "black", "maximum of G" = "#d7301f",
                                   "leading order" = "#fc8d59", "K/x" = "#2c7fb8"), name = NULL) +
    scale_linetype_manual(values = c("exact limit" = "solid", "maximum of G" = "dashed",
                                     "leading order" = "dotted", "K/x" = "dashed"), name = NULL) +
    labs(x = "x = c/n, seeding events per entry node", y = expression(gamma[infinity]^"*")) +
    theme_paper
  show(gB, "figS_gen_gamma_star_eta.pdf", w = PAGE, h = 4.8)
}

# ---- Gen C (main Fig. 5). Cut-off crossover: capped / uncapped vs n / n_c ----
# Rows: targeted mean degree (top), participation ratio Y2 (bottom).
if (file.exists("results/11_partC_crossover.csv")) {
  GCx <- read.csv("results/11_partC_crossover.csv") %>% filter(gamma > 2.31)
  gC <- GCx %>% pivot_longer(c(mu_rel, Y2_rel), names_to = "obs", values_to = "rel") %>%
    mutate(obs = factor(recode(obs, mu_rel = '"targeted degree"', Y2_rel = "Y[2]"),
                        levels = c('"targeted degree"', "Y[2]")),
           g_lab = paste0("gamma == ", gamma)) %>%
    ggplot(aes(n_over_nc, rel, colour = factor(cap), shape = law_lab(law))) +
    geom_line(aes(group = interaction(cap, law)), linewidth = 0.45) + geom_point(size = 1.1) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    facet_grid(obs ~ g_lab, scales = "free_y", labeller = label_parsed) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_viridis_d(name = expression(K[max]), end = 0.85) +
    scale_shape_discrete(name = NULL) +
    guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
    labs(x = expression(n / n[c]), y = "Capped / uncapped") +
    theme_paper + theme(legend.box = "vertical")
  show(gC, "figS_gen_crossover.pdf", w = COL, h = 4.2)
}
