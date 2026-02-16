## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(gsDesignNB)
library(data.table)
library(dplyr)

## ----params-------------------------------------------------------------------
# Randomization starts in Winter
rand_start <- as.Date("2024-01-01")

# Enrollment: 20 subjects over 6 months
enroll_rate <- data.frame(
  rate = 20 / 6,
  duration = 6
) # time unit: months?
# Note: In nb_sim_seasonal, we assumed rate units match max_followup units.
# Let's use YEARS as the time unit for consistency with the package conventions often used.
# If max_followup = 1 (year), then rates are per year.
# Enrollment duration = 0.5 years.

enroll_rate <- data.frame(
  rate = 20 / 0.5,
  duration = 0.5
)

# Seasonal Failure Rates (per year)
# Winter: High
# Spring: Low
# Summer: Low
# Fall: Medium-High

fail_rate <- data.frame(
  treatment = rep(c("Control", "Experimental"), each = 4),
  season = rep(c("Winter", "Spring", "Summer", "Fall"), 2),
  rate = c(
    # Control
    2.0, # Winter
    0.5, # Spring
    0.2, # Summer
    1.5, # Fall
    # Experimental (assume 30% reduction)
    2.0 * 0.7,
    0.5 * 0.7,
    0.2 * 0.7,
    1.5 * 0.7
  ),
  dispersion = 0.5 # Constant dispersion
)

# Dropout (5% per year)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.05, 0.05),
  duration = c(100, 100)
)

## ----sim----------------------------------------------------------------------
set.seed(123)

sim_data <- nb_sim_seasonal(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = 1, # 1 year
  randomization_start_date = rand_start,
  n = 20
)

# View structure of the output
head(sim_data)

## ----cut----------------------------------------------------------------------
# Cut data at 9 months (0.75 years)
cut_date <- 0.75

# Use a small gap (e.g., 7 days)
gap_days <- 7 / 365.25

analysis_data <- cut_data_by_date(sim_data, cut_date = cut_date, event_gap = gap_days)

head(analysis_data)

## ----summary------------------------------------------------------------------
# Summarize events by season
library(dplyr)
analysis_data %>%
  group_by(season, treatment) %>%
  summarize(
    Subjects = n_distinct(id),
    TotalExposure = sum(tte),
    TotalEvents = sum(events),
    Rate = sum(events) / sum(tte)
  )

## ----analysis-----------------------------------------------------------------
# Fit negative binomial GLM
# We include season and treatment in the model
fit <- MASS::glm.nb(
  events ~ treatment + season + offset(log(tte)),
  data = analysis_data[analysis_data$tte > 0, ]
)

# Summary of the model
summary(fit)

# Extract estimates
coef_summary <- summary(fit)$coefficients

# Seasonal effects (relative to reference season, likely Fall based on alphabetical order)
# Treatment effect (Experimental vs Control)
print(coef_summary)

# Estimated Rate Ratio (Experimental / Control)
rr <- exp(coef(fit)["treatmentExperimental"])
message("Estimated Rate Ratio (Experimental / Control): ", round(rr, 3))
message("True Design Rate Ratio: ", 0.7)

## ----information--------------------------------------------------------------
# Variance of the treatment effect coefficient
var_beta <- vcov(fit)["treatmentExperimental", "treatmentExperimental"]
message("Variance of Treatment Effect (log-scale): ", var_beta)

# Statistical Information
info <- 1 / var_beta
message("Statistical Information: ", info)

