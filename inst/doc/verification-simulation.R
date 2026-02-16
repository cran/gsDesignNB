## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(gsDesignNB)
library(data.table)
library(ggplot2)
library(gt)

## ----design-------------------------------------------------------------------
# Parameters
lambda1 <- 0.4
lambda2 <- 0.3
dispersion <- 0.5
power <- 0.9
alpha <- 0.025
dropout_rate <- 0.1 / 12
max_followup <- 12
trial_duration <- 24
event_gap <- 20 / 30.42 # 20 days

# Accrual targeting 90% power
# We provide relative rates (1:2) and the function scales them to achieve power
accrual_rate_rel <- c(1, 2)
accrual_duration <- c(6, 6)

design <- sample_size_nbinom(
  lambda1 = lambda1, lambda2 = lambda2, dispersion = dispersion,
  power = power,
  alpha = alpha, sided = 1,
  accrual_rate = accrual_rate_rel,
  accrual_duration = accrual_duration,
  trial_duration = trial_duration,
  dropout_rate = dropout_rate,
  max_followup = max_followup,
  event_gap = event_gap
)

# Extract calculated absolute accrual rates
accrual_rate <- design$accrual_rate

print(design)

## ----load_results-------------------------------------------------------------
# Load pre-computed simulation results
results_file <- system.file("extdata", "simulation_results.rds", package = "gsDesignNB")

if (results_file == "" && file.exists("../inst/extdata/simulation_results.rds")) {
  results_file <- "../inst/extdata/simulation_results.rds"
}

if (results_file != "") {
  sim_data <- readRDS(results_file)
  results <- sim_data$results
  design_ref <- sim_data$design
} else {
  # Fallback if data is not available (e.g. not installed yet)
  # This block allows the vignette to build without the data, but warns.
  warning("Simulation results not found. Skipping verification plots.")
  results <- NULL
  design_ref <- design
}

## ----summary_table, eval=!is.null(results)------------------------------------
# ---- Compute all metrics ----

# Number of simulations
n_sims <- sum(!is.na(results$estimate))

# Total Exposure (calendar follow-up time)
# Note: exposure is the same for both arms in the design (by symmetry)
theo_exposure <- design_ref$exposure[1]

# Check which column names are available in the results
# (Support both old and new naming conventions)
has_new_cols <- "exposure_total_control" %in% names(results)

if (has_new_cols) {
  obs_exposure_ctrl <- mean(results$exposure_total_control)
  obs_exposure_exp <- mean(results$exposure_total_experimental)
  obs_exposure_at_risk_ctrl <- mean(results$exposure_at_risk_control)
  obs_exposure_at_risk_exp <- mean(results$exposure_at_risk_experimental)
} else {
  # Legacy: old simulation used 'exposure_control' which was actually at-risk time
  obs_exposure_ctrl <- NA

  obs_exposure_exp <- NA
  obs_exposure_at_risk_ctrl <- mean(results$exposure_control)
  obs_exposure_at_risk_exp <- mean(results$exposure_experimental)
}

# Exposure at risk (time at risk excluding event gaps)
theo_exposure_at_risk_ctrl <- design_ref$exposure_at_risk_n1
theo_exposure_at_risk_exp <- design_ref$exposure_at_risk_n2

# Events by treatment group
theo_events_ctrl <- design_ref$events_n1
theo_events_exp <- design_ref$events_n2
obs_events_ctrl <- mean(results$events_control)
obs_events_exp <- mean(results$events_experimental)

# Treatment effect
true_rr <- lambda2 / lambda1
true_log_rr <- log(true_rr)
mean_log_rr <- mean(results$estimate, na.rm = TRUE)

# Variance
theo_var <- design_ref$variance
# Use median of SE^2 for robust estimate
median_se_sq <- median(results$se^2, na.rm = TRUE)
# Empirical variance of estimates
emp_var <- var(results$estimate, na.rm = TRUE)

# Power
theo_power <- design_ref$power
emp_power <- mean(results$p_value < design_ref$inputs$alpha, na.rm = TRUE)

# Sample size reproduction
z_alpha <- qnorm(1 - design_ref$inputs$alpha)
z_beta <- qnorm(design_ref$inputs$power)
n_sim_total <- design_ref$n_total
n_reproduced <- n_sim_total * (emp_var * (z_alpha + z_beta)^2) / (mean_log_rr^2)

# ---- Build summary table ----
summary_df <- data.frame(
  Metric = c(
    "Total Exposure (months) - Control",
    "Total Exposure (months) - Experimental",
    "Exposure at Risk (months) - Control",
    "Exposure at Risk (months) - Experimental",
    "Events per Subject - Control",
    "Events per Subject - Experimental",
    "Treatment Effect: log(RR)",
    "Variance of log(RR)",
    "Power",
    "Sample Size"
  ),
  Theoretical = c(
    theo_exposure,
    theo_exposure,
    theo_exposure_at_risk_ctrl,
    theo_exposure_at_risk_exp,
    theo_events_ctrl / (n_sim_total / 2),
    theo_events_exp / (n_sim_total / 2),
    true_log_rr,
    theo_var,
    theo_power,
    n_sim_total
  ),
  Simulated = c(
    obs_exposure_ctrl,
    obs_exposure_exp,
    obs_exposure_at_risk_ctrl,
    obs_exposure_at_risk_exp,
    obs_events_ctrl / (n_sim_total / 2),
    obs_events_exp / (n_sim_total / 2),
    mean_log_rr,
    median_se_sq,
    emp_power,
    n_reproduced
  ),
  stringsAsFactors = FALSE
)

summary_df$Difference <- summary_df$Simulated - summary_df$Theoretical
summary_df$Rel_Diff_Pct <- 100 * summary_df$Difference / abs(summary_df$Theoretical)

summary_df |>
  gt() |>
  tab_header(
    title = md("**Verification of sample_size_nbinom() Predictions**"),
    subtitle = paste0("Based on ", n_sims, " simulated trials")
  ) |>
  tab_row_group(
    label = md("**Sample Size**"),
    rows = Metric == "Sample Size"
  ) |>
  tab_row_group(
    label = md("**Power**"),
    rows = Metric == "Power"
  ) |>
  tab_row_group(
    label = md("**Variance**"),
    rows = Metric == "Variance of log(RR)"
  ) |>
  tab_row_group(
    label = md("**Treatment Effect**"),
    rows = Metric == "Treatment Effect: log(RR)"
  ) |>
  tab_row_group(
    label = md("**Events**"),
    rows = grepl("Events", Metric)
  ) |>
  tab_row_group(
    label = md("**Exposure**"),
    rows = grepl("Exposure", Metric)
  ) |>
  row_group_order(groups = c("**Exposure**", "**Events**", "**Treatment Effect**",
                              "**Variance**", "**Power**", "**Sample Size**")) |>
  fmt_number(columns = c(Theoretical, Simulated, Difference), decimals = 4) |>
  fmt_number(columns = Rel_Diff_Pct, decimals = 2) |>
  cols_label(
    Metric = "",
    Theoretical = "Theoretical",
    Simulated = "Simulated",
    Difference = "Difference",
    Rel_Diff_Pct = "Rel. Diff (%)"
  ) |>
  sub_missing(missing_text = "—")

## ----power_ci, eval=!is.null(results)-----------------------------------------
power_ci <- binom.test(sum(results$p_value < design_ref$inputs$alpha, na.rm = TRUE), 
   nrow(results))$conf.int
cat("95% CI for empirical power: [", round(power_ci[1], 4), ", ", round(power_ci[2], 4), "]\n", sep = "")
cat("Theoretical power:", round(theo_power, 4), "\n")

## ----z_score_dist, eval=!is.null(results)-------------------------------------
# Calculate Z-scores using estimated variance (Wald statistic)
# Z = (hat(delta) - 0) / SE_est
z_scores <- results$estimate / results$se

# Theoretical mean Z-score under the alternative
# E[Z] = log(RR) / sqrt(V_theo)
theo_se <- sqrt(theo_var)
theo_mean_z <- log(lambda2 / lambda1) / theo_se

# Critical value for 1-sided alpha (since we are looking at lower tail for RR < 1)
# However, the Z-scores here are negative (log(0.3/0.4) < 0).
# Rejection region is Z < qnorm(alpha)
crit_val <- qnorm(design_ref$inputs$alpha)

# Proportion of simulations rejecting the null
prop_reject <- mean(z_scores < crit_val, na.rm = TRUE)

# Plot
ggplot(data.frame(z = z_scores), aes(x = z)) +
  geom_density(aes(color = "Simulated Density"), linewidth = 1) +
  stat_function(
    fun = dnorm,
    args = list(mean = theo_mean_z, sd = 1),
    aes(color = "Theoretical Normal"),
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_vline(xintercept = crit_val, linetype = "dotted", color = "black") +
  annotate(
    "text",
    x = crit_val, y = 0.05,
    label = paste0("  Reject: ", round(prop_reject * 100, 1), "%"),
    hjust = 0, vjust = 0
  ) +
  labs(
    title = "Distribution of Wald Z-scores",
    subtitle = paste("Theoretical Mean Z =", round(theo_mean_z, 3)),
    x = "Z-score (Estimate / Estimated SE)",
    y = "Density",
    color = "Distribution"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

## ----log_rr_stats, eval=!is.null(results)-------------------------------------
# Theoretical mean and SD of log(RR)
# Let's also add median, skewness, and kurtosis
theo_mean_log_rr <- log(lambda2 / lambda1)
theo_sd_log_rr <- sqrt(theo_var)
emp_mean_log_rr <- mean(results$estimate, na.rm = TRUE)
emp_sd_log_rr <- sd(results$estimate, na.rm = TRUE)
emp_median_log_rr <- median(results$estimate, na.rm = TRUE)

# Skewness and Kurtosis (using trimmed data to reduce outlier influence)
# Trim extreme 1% on each tail for robust estimates
ests <- results$estimate[!is.na(results$estimate)]
q_low <- quantile(ests, 0.01)
q_high <- quantile(ests, 0.99)
ests_trimmed <- ests[ests >= q_low & ests <= q_high]
n_trimmed <- length(ests_trimmed)
m3 <- sum((ests_trimmed - mean(ests_trimmed))^3) / n_trimmed
m4 <- sum((ests_trimmed - mean(ests_trimmed))^4) / n_trimmed
s2 <- sum((ests_trimmed - mean(ests_trimmed))^2) / n_trimmed
emp_skew_log_rr <- m3 / s2^(3/2)
emp_kurt_log_rr <- m4 / s2^2

comparison_log_rr <- data.frame(
  Metric = c("Mean", "SD", "Median", "Skewness (trimmed)", "Kurtosis (trimmed)"),
  Theoretical = c(theo_mean_log_rr, theo_sd_log_rr, theo_mean_log_rr, 0, 3),
  Simulated = c(emp_mean_log_rr, emp_sd_log_rr, emp_median_log_rr, emp_skew_log_rr, emp_kurt_log_rr),
  Difference = c(
    emp_mean_log_rr - theo_mean_log_rr,
    emp_sd_log_rr - theo_sd_log_rr,
    emp_median_log_rr - theo_mean_log_rr,
    emp_skew_log_rr - 0,
    emp_kurt_log_rr - 3
  )
)
# Display table
comparison_log_rr |>
  gt() |>
  tab_header(title = md("**Comparison of log(RR) Statistics**")) |>
  fmt_number(columns = where(is.numeric), decimals = 4)

