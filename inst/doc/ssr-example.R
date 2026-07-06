## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(gsDesignNB)
library(gsDesign)
library(data.table)
library(MASS)
library(gt)

## ----initial_design-----------------------------------------------------------
# Planned parameters
lambda1_plan <- 0.1
lambda2_plan <- 0.075
k_plan <- 0.5
power_plan <- 0.9
alpha_plan <- 0.025
accrual_rate_plan <- 20
accrual_dur_plan <- 12
trial_dur_plan <- 24

# Calculate sample size
design_plan <- sample_size_nbinom(
  lambda1 = lambda1_plan,
  lambda2 = lambda2_plan,
  dispersion = k_plan,
  power = power_plan,
  alpha = alpha_plan,
  accrual_rate = accrual_rate_plan,
  accrual_duration = accrual_dur_plan,
  trial_duration = trial_dur_plan,
  max_followup = 12
)

# Convert to group sequential design
gs_plan <- design_plan |> 
   gsNBCalendar(
     k = 2,
     test.type = 4, # Non-binding futility
     alpha = alpha_plan,
     sfu = sfHSD,
     sfupar = -2, # Moderately aggressive alpha-spending
     sfl = sfHSD,
     sflpar = 1, # Pocock-like futility spending
     analysis_times = c(accrual_dur_plan - 2, trial_dur_plan)
   ) |> gsDesignNB::toInteger() # Round to integer sample size


## -----------------------------------------------------------------------------
summary(gs_plan)

## -----------------------------------------------------------------------------
gsBoundSummary(gs_plan,
    deltaname = "RR",
    logdelta = TRUE,
    Nname = "Information",
    timename = "Month",
    digits = 4,
    ddigits = 2) |> gt() |>
  tab_header(
    title = "Group Sequential Design Bounds for Negative Binomial Outcome",
    subtitle = paste0(
      "N = ", ceiling(gs_plan$n_total[gs_plan$k]),
      ", Expected events = ", round(gs_plan$nb_design$total_events, 1)
    )
  )

## ----helper_workflow----------------------------------------------------------
sim_helper <- sim_ssr_nbinom(
  n_sims = 20,
  enroll_rate = data.frame(rate = accrual_rate_plan, duration = accrual_dur_plan * 2),
  fail_rate = data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0.08, 0.06),
    dispersion = 0.65
  ),
  dropout_rate = data.frame(
    treatment = c("Control", "Experimental"),
    rate = c(0, 0),
    duration = c(100, 100)
  ),
  max_followup = trial_dur_plan,
  design = gs_plan,
  n_max = ceiling(1.4 * gs_plan$n_total[gs_plan$k]),
  metadata = list(scenario = "Lower control rate, higher dispersion"),
  seed = 1234
)

helper_summary <- summarize_ssr_sim(sim_helper, by = c("scenario", "strategy"))
helper_summary$trial_summary |>
  gt() |>
  tab_header(
    title = "High-level SSR Simulation Summary",
    subtitle = "Expected participants with events and expected events observed are included directly"
  )

## ----simulation---------------------------------------------------------------
set.seed(1234)

# Actual parameters
lambda1_true <- 0.08 # Lower event rates, same rr
lambda2_true <- 0.06
k_true <- 0.65 # Higher dispersion

# Enrollment and rates for simulation
# We simulate a larger pool to allow for potential sample size increase
enroll_rate <- data.frame(rate = accrual_rate_plan, duration = accrual_dur_plan * 2) 
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(lambda1_true, lambda2_true),
  dispersion = k_true
)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0, 0),
  duration = c(100, 100)
)

sim_data <- nb_sim(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = trial_dur_plan, 
  n = 600 
)

# Limit to planned enrollment for the initial cut
# We will "open" more enrollment if needed
sim_data_planned <- sim_data[1:ceiling(design_plan$n_total), ]

## ----interim_cut--------------------------------------------------------------
interim_time <- 10
interim_data <- cut_data_by_date(sim_data_planned, cut_date = interim_time)

# Summary
table(interim_data$treatment)
sum(interim_data$events)
mean(interim_data$tte)

## ----info_calc----------------------------------------------------------------
# Blinded SSR
blinded_res <- blinded_ssr(
  data = interim_data,
  ratio = 1,
  lambda1_planning = lambda1_plan,
  lambda2_planning = lambda2_plan,
  power = power_plan,
  alpha = alpha_plan,
  accrual_rate = accrual_rate_plan,
  accrual_duration = accrual_dur_plan,
  trial_duration = trial_dur_plan
)

# Unblinded SSR
ssr_res <- unblinded_ssr(
  data = interim_data,
  ratio = 1,
  lambda1_planning = lambda1_plan,
  lambda2_planning = lambda2_plan,
  power = power_plan,
  alpha = alpha_plan,
  accrual_rate = accrual_rate_plan,
  accrual_duration = accrual_dur_plan,
  trial_duration = trial_dur_plan
)

cat("Blinded SSR N:", ceiling(blinded_res$n_total_blinded), "\n")
cat("Unblinded SSR N:", ceiling(ssr_res$n_total_unblinded), "\n")
print(ssr_res)

## ----update_design------------------------------------------------------------
n_planned <- ceiling(design_plan$n_total)
n_estimated <- ceiling(ssr_res$n_total_unblinded)

# Ensure we don't decrease sample size
n_final <- max(n_planned, n_estimated)

if (n_final > n_planned) {
  cat("Increasing sample size from", n_planned, "to", n_final, "\n")
  
  # Calculate new durations based on constant accrual rate
  # We extend enrollment to reach the new target N
  new_accrual_dur <- n_final / accrual_rate_plan
  
  # We maintain the planned max follow-up of 12 months
  # So the trial duration extends by the same amount as the accrual duration
  new_trial_dur <- new_accrual_dur + 12 
  
  cat("New Accrual Duration:", round(new_accrual_dur, 2), "months\n")
  cat("New Trial Duration:", round(new_trial_dur, 2), "months\n")
} else {
  cat("No sample size increase needed.\n")
  new_accrual_dur <- accrual_dur_plan
  new_trial_dur <- trial_dur_plan
  n_final <- n_planned
}

## ----final_analysis-----------------------------------------------------------
# Update the dataset to include the additional patients
sim_data_final <- sim_data[1:n_final, ]

# Cut at final analysis time
final_data <- cut_data_by_date(sim_data_final, cut_date = new_trial_dur)

# Average exposure
mean(final_data$tte)

# Fit final model
fit_final <- glm.nb(events ~ treatment + offset(log(tte)), data = final_data)
summary(fit_final)

# Calculate final information
var_beta <- vcov(fit_final)[2, 2]
final_info <- 1 / var_beta
cat("Final Information:", final_info, "\n")
cat("Target Information:", ssr_res$target_info, "\n")

## ----update_bounds------------------------------------------------------------
# Information fraction at final analysis
# If final_info < target_info, fraction < 1. 
# However, for the final analysis, we typically want to spend all remaining alpha.
# We can treat the current info as the "maximum" info for the spending function.

# Example group sequential design (O'Brien-Fleming spending)
# We assume the interim was at the observed information fraction
info_fractions <- c(ssr_res$info_fraction, final_info / ssr_res$target_info)

# If final fraction < 1, we can force it to 1 to spend all alpha, 
# effectively accepting the lower power.
# Or we can use the actual information fractions in a spending function approach.

# Let's assume we want to spend all alpha at the final analysis regardless of info.
# We set the final timing to 1.
info_fractions_spending <- c(ssr_res$info_fraction, 1)

# Update bounds
gs_update <- gsDesign(
  k = 2,
  test.type = 1,
  alpha = alpha_plan,
  sfu = sfLDOF,
  usTime = info_fractions_spending,
  n.I = c(ssr_res$unblinded_info, final_info)
)

gsBoundSummary(gs_update,
    deltaname = "RR",
    logdelta = TRUE,
    Nname = "Information",
    timename = "Month",
    digits = 4,
    ddigits = 2) |> gt() |>
  tab_header(
    title = "Updated Group Sequential Design Bounds",
    subtitle = paste0(
      "Final Information = ", round(final_info, 2)
    )
  )

# Test statistic
z_score <- coef(fit_final)["treatmentExperimental"] / sqrt(var_beta)
# Note: gsDesign assumes Z > 0 for efficacy. 
# Our model gives log(RR) < 0 for efficacy.
# So Z_stat = - (log(RR)) / SE
z_stat <- -coef(fit_final)["treatmentExperimental"] / sqrt(var_beta)

cat("Z-statistic:", z_stat, "\n")
cat("Final Bound:", gs_update$upper$bound[2], "\n")

if (z_stat > gs_update$upper$bound[2]) {
  cat("Reject Null Hypothesis\n")
} else {
  cat("Fail to Reject Null Hypothesis\n")
}

