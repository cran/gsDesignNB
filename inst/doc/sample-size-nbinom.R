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
  event_gap = 30 / 365.25
)

