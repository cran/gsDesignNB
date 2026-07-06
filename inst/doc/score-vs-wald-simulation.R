## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(gsDesignNB)
library(data.table)
library(ggplot2)

## -----------------------------------------------------------------------------
summary_file <- system.file("extdata", "score_sweep_summary.rds",
                            package = "gsDesignNB")
if (summary_file == "" && file.exists("../inst/extdata/score_sweep_summary.rds")) {
  summary_file <- "../inst/extdata/score_sweep_summary.rds"
}

raw_file <- system.file("extdata", "score_sweep_results.rds",
                        package = "gsDesignNB")
if (raw_file == "" && file.exists("../inst/extdata/score_sweep_results.rds")) {
  raw_file <- "../inst/extdata/score_sweep_results.rds"
}

if (summary_file != "") {
  res <- readRDS(summary_file)
  using_summary_cache <- TRUE
} else if (raw_file != "") {
  res <- readRDS(raw_file)
  using_summary_cache <- FALSE
} else {
  stop("Precomputed score sweep summary not found.")
}
config <- res$config
scenarios <- as.data.table(res$scenarios)
base_grid <- as.data.table(res$base_grid)

## -----------------------------------------------------------------------------
cat(sprintf(
  "Expanded scenarios: %d | Power sims: %s | Null sims: %s | RR: %.2f | alpha: %.3f\n",
  nrow(scenarios),
  format(config$n_sims_power, big.mark = ","),
  format(config$n_sims_null, big.mark = ","),
  config$rr_power,
  config$alpha
))
cat(sprintf(
  "Cache: %s\n",
  if (using_summary_cache) "compact summary" else "full raw simulation output"
))

## -----------------------------------------------------------------------------
base_display <- base_grid[, .(
  `Control rate` = lambda1,
  `Dispersion (k)` = k,
  `Event gap (days)` = gap_days,
  `N (Wald sizing)` = n_wald,
  `N (Score sizing)` = n_score,
  `Wald - Score` = n_wald - n_score
)]

knitr::kable(
  base_display,
  caption = "Base scenario grid with sample sizes by method",
  digits = 2
)

## -----------------------------------------------------------------------------
null_dt <- as.data.table(res$null_summary)
null_long <- melt(
  null_dt,
  id.vars = c("lambda1", "k", "gap_days", "n_total", "sizing"),
  measure.vars = c("rate_wald", "rate_score"),
  variable.name = "test",
  value.name = "rejection_rate"
)
null_long[, test := fifelse(test == "rate_wald", "Wald", "Score")]

se_long <- melt(
  null_dt,
  id.vars = c("lambda1", "k", "gap_days", "sizing"),
  measure.vars = c("se_wald", "se_score"),
  variable.name = "test",
  value.name = "se"
)
se_long[, test := fifelse(test == "se_wald", "Wald", "Score")]
null_long <- merge(null_long, se_long,
                   by = c("lambda1", "k", "gap_days", "sizing", "test"))
null_long[, combo := paste0(sizing, "-sized / ", test, " test")]
null_long[, `:=`(
  above_nominal_95 = rejection_rate - 1.96 * se > config$alpha,
  below_nominal_95 = rejection_rate + 1.96 * se < config$alpha
)]

## -----------------------------------------------------------------------------
type1_summary <- null_long[, .(
  `Scenarios` = .N,
  `Minimum` = min(rejection_rate),
  `Mean` = mean(rejection_rate),
  `Maximum` = max(rejection_rate),
  `Above nominal beyond MC error` = sum(above_nominal_95),
  `Below nominal beyond MC error` = sum(below_nominal_95)
), by = .(`Sizing` = sizing, `Test` = test)]

knitr::kable(
  type1_summary[order(Sizing, Test)],
  caption = "Type I error synopsis across the scenario grid",
  digits = 4
)

## -----------------------------------------------------------------------------
null_display <- null_long[order(lambda1, k, gap_days, sizing, test),
  .(
    `Control rate` = lambda1,
    Dispersion = k,
    `Gap (days)` = gap_days,
    Sizing = sizing,
    N = n_total,
    Test = test,
    `Type I error` = round(rejection_rate, 4),
    SE = round(se, 4)
  )
]

knitr::kable(
  null_display,
  caption = sprintf(
    "Type I error rate: nominal alpha = %.3f, %s null sims/scenario",
    config$alpha,
    format(config$n_sims_null, big.mark = ",")
  )
)

## ----fig.width=8, fig.height=6------------------------------------------------
null_long[, scenario := paste0("λ₁=", lambda1, " k=", k)]

p_null <- ggplot(null_long,
  aes(x = scenario, y = rejection_rate,
      color = combo, shape = test)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = rejection_rate - 1.96 * se,
                     ymax = rejection_rate + 1.96 * se),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = config$alpha, linetype = "dashed", color = "grey40") +
  facet_wrap(~ paste0("Gap = ", gap_days, "d")) +
  labs(
    title = "Type I error: sizing method × test type",
    x = NULL, y = "Rejection rate",
    color = "Sizing / Test", shape = "Test"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_null

## -----------------------------------------------------------------------------
power_dt <- as.data.table(res$power_summary)
power_long <- melt(
  power_dt,
  id.vars = c("lambda1", "k", "gap_days", "n_total", "sizing"),
  measure.vars = c("rate_wald", "rate_score"),
  variable.name = "test",
  value.name = "power"
)
power_long[, test := fifelse(test == "rate_wald", "Wald", "Score")]

se_power <- melt(
  power_dt,
  id.vars = c("lambda1", "k", "gap_days", "sizing"),
  measure.vars = c("se_wald", "se_score"),
  variable.name = "test",
  value.name = "se"
)
se_power[, test := fifelse(test == "se_wald", "Wald", "Score")]
power_long <- merge(power_long, se_power,
                    by = c("lambda1", "k", "gap_days", "sizing", "test"))
power_long[, combo := paste0(sizing, "-sized / ", test, " test")]

## -----------------------------------------------------------------------------
power_summary <- power_long[, .(
  `Scenarios` = .N,
  `Minimum` = min(power),
  `Mean` = mean(power),
  `Maximum` = max(power),
  `Below 90%` = sum(power < config$power_target)
), by = .(`Sizing` = sizing, `Test` = test)]

knitr::kable(
  power_summary[order(Sizing, Test)],
  caption = "Power synopsis across the scenario grid",
  digits = 4
)

## -----------------------------------------------------------------------------
power_display <- power_long[order(lambda1, k, gap_days, sizing, test),
  .(
    `Control rate` = lambda1,
    Dispersion = k,
    `Gap (days)` = gap_days,
    Sizing = sizing,
    N = n_total,
    Test = test,
    Power = round(power, 4),
    SE = round(se, 4)
  )
]

knitr::kable(
  power_display,
  caption = sprintf(
    "Power: RR = %.2f, %s power sims/scenario",
    config$rr_power,
    format(config$n_sims_power, big.mark = ",")
  )
)

## ----fig.width=8, fig.height=6------------------------------------------------
power_long[, scenario := paste0("λ₁=", lambda1, " k=", k)]

p_power <- ggplot(power_long,
  aes(x = scenario, y = power,
      color = combo, shape = test)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = power - 1.96 * se,
                     ymax = power + 1.96 * se),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = config$power_target, linetype = "dashed", color = "grey40") +
  facet_wrap(~ paste0("Gap = ", gap_days, "d")) +
  labs(
    title = "Power: sizing method × test type",
    subtitle = sprintf("Target = %.0f%%, RR = %.2f",
                        100 * config$power_target, config$rr_power),
    x = NULL, y = "Power",
    color = "Sizing / Test", shape = "Test"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_power

## ----fig.width=9, fig.height=7------------------------------------------------
if (!is.null(res$z_density_null)) {
  z_density_null <- as.data.table(res$z_density_null)
} else {
  z_null <- as.data.table(res$z_sample_null)
  sc_info <- data.table(
    scenario_id = seq_len(nrow(scenarios)),
    scenarios[, .(lambda1, k, gap_days, sizing)]
  )
  z_null <- merge(z_null, sc_info, by = "scenario_id")
  z_null[, label := sprintf("l1=%.2f k=%.1f gap=%dd (%s)",
                             lambda1, k, gap_days, sizing)]

  z_null_long <- melt(
    z_null,
    id.vars = c("scenario_id", "label", "sizing"),
    measure.vars = c("z_wald", "z_score"),
    variable.name = "test",
    value.name = "z"
  )
  z_null_long[, test := fifelse(test == "z_wald", "Wald", "Score")]
  z_null_long <- z_null_long[is.finite(z)]

  z_density_null <- z_null_long[, {
    dens <- stats::density(z, from = -4, to = 4, n = 128)
    .(z = dens$x, density = dens$y)
  }, by = .(scenario_id, label, sizing, test)]
}

normal_curve <- data.table(
  z = seq(-4, 4, length.out = 128),
  density = dnorm(seq(-4, 4, length.out = 128))
)

p_z <- ggplot(z_density_null, aes(x = z, y = density, color = test)) +
  geom_line(linewidth = 0.6) +
  geom_line(data = normal_curve, aes(x = z, y = density),
            inherit.aes = FALSE, color = "black", linetype = "dashed",
            linewidth = 0.4) +
  facet_wrap(~ label, scales = "free_y") +
  labs(
    title = "Null Z-score densities: Wald vs Score",
    subtitle = "Dashed line = N(0,1) reference",
    x = "Z-statistic", y = "Density",
    color = "Test"
  ) +
  theme_minimal() +
  coord_cartesian(xlim = c(-4, 4))

p_z

## -----------------------------------------------------------------------------
null_fb <- as.data.table(res$null_summary)
fb_display <- null_fb[, .(
  `Control rate` = lambda1,
  Dispersion = k,
  `Gap (days)` = gap_days,
  Sizing = sizing,
  `Poisson (Wald)` = round(pct_fallback_poisson_wald, 1),
  `MoM (Wald)` = round(pct_fallback_mom_wald, 1),
  `Poisson (Score)` = round(pct_fallback_poisson_score, 1),
  `MoM (Score)` = round(pct_fallback_mom_score, 1)
)]

knitr::kable(
  fb_display,
  caption = "Fallback method frequency (%, null sims)",
  digits = 1
)

