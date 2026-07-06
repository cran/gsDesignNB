## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load-packages, message=FALSE---------------------------------------------
library(gsDesignNB)
library(data.table)
library(ggplot2)
library(MASS)

## ----load-data-1--------------------------------------------------------------
# Load the problematic dataset
# Try package location first, fall back to relative path for development
data_path <- system.file("extdata", "problematic_dataset.rds", package = "gsDesignNB")
if (data_path == "") {
  # During development, find the package root
  data_path <- file.path("..", "inst", "extdata", "problematic_dataset.rds")
}
cut_data <- readRDS(data_path)
dt <- as.data.table(cut_data)

## ----event-summary-1----------------------------------------------------------
# Overall summary
cat("Total subjects:", nrow(dt), "\n")
cat("Total events:", sum(dt$events), "\n")
cat("Mean events per subject:", round(mean(dt$events), 3), "\n")
cat("Variance of events:", round(var(dt$events), 3), "\n")
cat("Variance/Mean ratio:", round(var(dt$events) / mean(dt$events), 3), "\n")

# By treatment
dt[, .(
  n = .N,
  total_events = sum(events),
  mean_events = round(mean(events), 3),
  var_events = round(var(events), 3),
  pct_zero_events = round(100 * mean(events == 0), 1)
), by = treatment]

## ----event-distribution-1-----------------------------------------------------
# Event count distribution
event_dist <- dt[, .(count = .N), by = .(treatment, events)]
event_dist[order(treatment, events)]

## ----event-rates-1------------------------------------------------------------
dt[, event_rate := events / tte]

# Summary of event rates
dt[, .(
  n = .N,
  mean_rate = round(mean(event_rate), 4),
  median_rate = round(median(event_rate), 4),
  sd_rate = round(sd(event_rate), 4),
  pct_zero_rate = round(100 * mean(event_rate == 0), 1)
), by = treatment]

## ----violin-plot-1, fig.height=6----------------------------------------------
# Add overall group for comparison
dt_plot <- rbind(
  dt[, .(treatment = treatment, event_rate = event_rate)],
  dt[, .(treatment = "Overall (Pooled)", event_rate = event_rate)]
)

ggplot(dt_plot, aes(x = treatment, y = event_rate, fill = treatment)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
  labs(
    title = "Distribution of Patient-Level Event Rates",
    subtitle = "Dataset with Near-Zero Blinded Information",
    x = "Treatment Group",
    y = "Event Rate (events / follow-up time)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Set2")

## ----info-calc-1--------------------------------------------------------------
# Blinded information calculation
blinded_res <- calculate_blinded_info(
  cut_data,
  ratio = 1,
  lambda1_planning = 1.5/12,
  lambda2_planning = 1.0/12
)

cat("Blinded Information Results:\n")
cat("  blinded_info:", blinded_res$blinded_info, "\n")
cat("  dispersion_blinded:", blinded_res$dispersion_blinded, "\n")
cat("  lambda_blinded:", blinded_res$lambda_blinded, "\n")

# Unblinded calculation via mutze_test
mutze_res <- mutze_test(cut_data)

cat("\nUnblinded (Mutze Test) Results:\n")
cat("  method:", mutze_res$method, "\n")
cat("  SE:", round(mutze_res$se, 4), "\n")
cat("  unblinded_info (1/SE^2):", round(1/mutze_res$se^2, 2), "\n")
cat("  dispersion:", mutze_res$dispersion, "\n")

## ----glm-analysis-1-----------------------------------------------------------
df <- as.data.frame(cut_data)
df <- df[df$tte > 0, ]

# Blinded fit (no treatment effect)
fit_blinded <- suppressWarnings(MASS::glm.nb(events ~ 1 + offset(log(tte)), data = df))

cat("Blinded glm.nb fit:\n")
cat("  theta:", fit_blinded$theta, "\n")
cat("  1/theta (dispersion):", 1/fit_blinded$theta, "\n")

# Unblinded fit (with treatment effect)
fit_unblinded <- suppressWarnings(MASS::glm.nb(events ~ treatment + offset(log(tte)), data = df))

cat("\nUnblinded glm.nb fit:\n")
cat("  theta:", fit_unblinded$theta, "\n")
cat("  1/theta (dispersion):", 1/fit_unblinded$theta, "\n")

## ----mom-analysis-1-----------------------------------------------------------
# Calculate MoM estimates
mom_res <- estimate_nb_mom(df)

cat("Method of Moments (Blinded) Estimation:\n")
cat("  Lambda (Rate):", round(mom_res$lambda, 4), "\n")
cat("  Dispersion (k):", round(mom_res$dispersion, 4), "\n")

## ----load-data-2--------------------------------------------------------------
# Load the historical extreme-high dataset
data_path_2 <- system.file("extdata", "extreme_high_dataset.rds", package = "gsDesignNB")
if (data_path_2 == "") {
  data_path_2 <- file.path("..", "inst", "extdata", "extreme_high_dataset.rds")
}
cut_data_2 <- readRDS(data_path_2)
dt2 <- as.data.table(cut_data_2)

## ----event-summary-2----------------------------------------------------------
# Overall summary
cat("Total subjects:", nrow(dt2), "\n")
cat("Total events:", sum(dt2$events), "\n")
cat("Mean events per subject:", round(mean(dt2$events), 3), "\n")
cat("Variance of events:", round(var(dt2$events), 3), "\n")
cat("Variance/Mean ratio:", round(var(dt2$events) / mean(dt2$events), 3), "\n")

# By treatment
dt2[, .(
  n = .N,
  total_events = sum(events),
  mean_events = round(mean(events), 3),
  var_events = round(var(events), 3),
  pct_zero_events = round(100 * mean(events == 0), 1)
), by = treatment]

## ----info-calc-2--------------------------------------------------------------
# Blinded information calculation
blinded_res_2 <- calculate_blinded_info(
  cut_data_2,
  ratio = 1,
  lambda1_planning = 1.5/12,
  lambda2_planning = 1.0/12
)

cat("Blinded Information Results:\n")
cat("  blinded_info:", format(blinded_res_2$blinded_info, scientific = TRUE), "\n")
cat("  dispersion_blinded:", format(blinded_res_2$dispersion_blinded, scientific = TRUE), "\n")
cat("  lambda_blinded:", round(blinded_res_2$lambda_blinded, 4), "\n")

# Unblinded calculation via mutze_test
mutze_res_2 <- mutze_test(cut_data_2)

cat("\nUnblinded (Mutze Test) Results:\n")
cat("  method:", mutze_res_2$method, "\n")
cat("  SE:", round(mutze_res_2$se, 4), "\n")
cat("  unblinded_info (1/SE^2):", round(1/mutze_res_2$se^2, 2), "\n")
cat("  dispersion:", mutze_res_2$dispersion, "\n")

## ----glm-analysis-2-----------------------------------------------------------
df2 <- as.data.frame(cut_data_2)
df2 <- df2[df2$tte > 0, ]

# Blinded fit (no treatment effect)
fit_blinded_2 <- tryCatch(
  suppressWarnings(MASS::glm.nb(events ~ 1 + offset(log(tte)), data = df2)),
  error = function(e) NULL
)

if (!is.null(fit_blinded_2)) {
  cat("Blinded glm.nb fit:\n")
  cat("  theta:", format(fit_blinded_2$theta, scientific = TRUE), "\n")
  cat("  1/theta (dispersion):", format(1/fit_blinded_2$theta, scientific = TRUE), "\n")
} else {
  cat("Blinded glm.nb fit FAILED\n")
}

# Unblinded fit (with treatment effect)
fit_unblinded_2 <- tryCatch(
  suppressWarnings(MASS::glm.nb(events ~ treatment + offset(log(tte)), data = df2)),
  error = function(e) NULL
)

if (!is.null(fit_unblinded_2)) {
  cat("\nUnblinded glm.nb fit:\n")
  cat("  theta:", round(fit_unblinded_2$theta, 4), "\n")
  cat("  1/theta (dispersion):", round(1/fit_unblinded_2$theta, 4), "\n")
} else {
  cat("\nUnblinded glm.nb fit FAILED - this explains the Poisson fallback\n")
}

## ----mom-analysis-2-----------------------------------------------------------
# Calculate MoM estimates
mom_res_2 <- estimate_nb_mom(df2)

cat("Method of Moments (Blinded) Estimation:\n")
cat("  Lambda (Rate):", round(mom_res_2$lambda, 4), "\n")
cat("  Dispersion (k):", mom_res_2$dispersion, "\n")

## ----session-info-------------------------------------------------------------
sessionInfo()

