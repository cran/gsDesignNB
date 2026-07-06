## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(gsDesignNB)

## ----fig.width=8, fig.height=4------------------------------------------------
par(mfrow = c(1, 3))
k_values <- c(0, 0.5, 1)
for (k in k_values) {
  mu <- 5
  x <- 0:15
  if (k == 0) {
    probs <- dpois(x, lambda = mu)
  } else {
    size <- 1 / k
    probs <- dnbinom(x, size = size, mu = mu)
  }
  barplot(probs,
    names.arg = x, horiz = TRUE, main = paste("k =", k),
    xlab = "Probability", ylab = "Event Count", las = 1, xlim = c(0, 0.2)
  )
}

## ----gap-correction-table-----------------------------------------------------
# Show correction magnitude for various parameter combos
params <- expand.grid(
  lambda = c(0.3, 0.5, 1.0, 2.0),
  k = c(0.1, 0.3, 0.5, 1.0),
  gap = c(0.25, 0.5, 1.0)
)
params$naive_eff <- params$lambda / (1 + params$lambda * params$gap)
params$correction <- 1 - params$k * params$lambda * params$gap /
  (1 + params$lambda * params$gap)^2
params$corrected_eff <- params$naive_eff * params$correction
params$pct_reduction <- round(100 * (1 - params$correction), 2)

# Show a subset
subset_df <- params[params$pct_reduction > 0.5, ]
subset_df <- subset_df[order(-subset_df$pct_reduction), ]
knitr::kable(
  head(subset_df[, c("lambda", "k", "gap", "naive_eff", "corrected_eff", "pct_reduction")], 12),
  digits = 4, row.names = FALSE,
  col.names = c("lambda", "k", "gap", "Naive eff. rate", "Corrected eff. rate", "Reduction (%)")
)

## ----exposure-verification, fig.width=7, fig.height=5-------------------------
set.seed(42)

# Design parameters
lambda1 <- 0.5
lambda2 <- 0.3
dispersion <- 0.3
dropout_rate_val <- 0.05
max_followup_val <- 8
trial_duration <- 12

# Piecewise accrual: 5/month for 4 months, then 15/month for 4 months
accrual_rate_vec <- c(5, 15)
accrual_duration_vec <- c(4, 4)

# Theoretical calculation
design <- sample_size_nbinom(
  lambda1 = lambda1, lambda2 = lambda2, dispersion = dispersion,
  power = 0.8,
  accrual_rate = accrual_rate_vec,
  accrual_duration = accrual_duration_vec,
  trial_duration = trial_duration,
  dropout_rate = dropout_rate_val,
  max_followup = max_followup_val
)

theo_exposure <- design$exposure[1] # Same for both groups (same dropout)

# Simulate many trials and compute average exposure
enroll_rate <- data.frame(rate = design$accrual_rate, duration = accrual_duration_vec)
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(lambda1, lambda2),
  dispersion = dispersion
)
dropout <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = rep(dropout_rate_val, 2),
  duration = rep(100, 2)
)

n_sims <- 100
sim_exposures <- numeric(n_sims)

for (i in seq_len(n_sims)) {
  sim_data <- nb_sim(
    enroll_rate = enroll_rate,
    fail_rate = fail_rate,
    dropout_rate = dropout,
    max_followup = max_followup_val,
    n = design$n_total
  )
  cut <- cut_data_by_date(sim_data, cut_date = trial_duration)
  sim_exposures[i] <- mean(cut$tte_total)
}

cat("Theoretical average exposure:", round(theo_exposure, 4), "\n")
cat("Simulated average exposure (mean over", n_sims, "trials):", round(mean(sim_exposures), 4), "\n")
cat("Simulation SE:", round(sd(sim_exposures) / sqrt(n_sims), 4), "\n")
cat("Relative difference:", round(100 * (mean(sim_exposures) - theo_exposure) / theo_exposure, 2), "%\n")

# Histogram of simulated average exposures
hist(sim_exposures,
  breaks = 30, col = "steelblue", border = "white",
  main = "Distribution of simulated average exposure",
  xlab = "Average exposure per subject"
)
abline(v = theo_exposure, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Theoretical"), col = "red", lty = 2, lwd = 2)

## -----------------------------------------------------------------------------
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  alpha = 0.025,
  sided = 1,
  accrual_rate = 10, # arbitrary, just for average exposure
  accrual_duration = 12,
  trial_duration = 12
)

## -----------------------------------------------------------------------------
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  accrual_rate = c(5, 10),
  accrual_duration = c(3, 3),
  trial_duration = 12
)

## -----------------------------------------------------------------------------
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  accrual_rate = c(5, 10),
  accrual_duration = c(3, 3),
  trial_duration = 12,
  dropout_rate = 0.05,
  max_followup = 6
)

## -----------------------------------------------------------------------------
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  accrual_rate = c(5, 10),
  accrual_duration = c(3, 3),
  trial_duration = 12,
  dropout_rate = c(0.10, 0.05), # Control, Treatment
  max_followup = 6
)

## -----------------------------------------------------------------------------
# Store the result from the previous calculation
design_result <- sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  power = 0.8,
  accrual_rate = c(5, 10),
  accrual_duration = c(3, 3),
  trial_duration = 12,
  dropout_rate = 0.05,
  max_followup = 6
)

# Use the computed accrual rates to calculate power for a smaller effect size
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.4, # Smaller effect size
  dispersion = 0.1,
  power = NULL, # Request power calculation
  accrual_rate = design_result$accrual_rate, # Use computed rates
  accrual_duration = c(3, 3),
  trial_duration = 12,
  dropout_rate = 0.05,
  max_followup = 6
)

## -----------------------------------------------------------------------------
sample_size_nbinom(
  lambda1 = 0.5,
  lambda2 = 0.3,
  dispersion = 0.1,
  ratio = 2,
  accrual_rate = 10,
  accrual_duration = 12,
  trial_duration = 12
)

## -----------------------------------------------------------------------------
sample_size_nbinom(
  lambda1 = 2.0,
  lambda2 = 1.0,
  dispersion = 0.1,
  power = 0.8,
  accrual_rate = 10,
  accrual_duration = 12,
  trial_duration = 12,
  event_gap = 20 / 365.25
)

