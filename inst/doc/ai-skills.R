## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(gsDesignNB)

## -----------------------------------------------------------------------------
event_gap_months <- 28 / 30.4375

design_args <- list(
  lambda1 = 0.08,
  lambda2 = 0.056,
  dispersion = 0.6,
  power = 0.80,
  alpha = 0.025,
  sided = 1,
  accrual_rate = 10,
  accrual_duration = 18,
  trial_duration = 30,
  dropout_rate = 0.01,
  max_followup = 12,
  event_gap = event_gap_months
)

## -----------------------------------------------------------------------------
wald_design <- do.call(
  sample_size_nbinom,
  c(design_args, list(test_type = "wald"))
)

score_design <- do.call(
  sample_size_nbinom,
  c(design_args, list(test_type = "score"))
)

design_comparison <- data.frame(
  test_type = c(wald_design$test_type, score_design$test_type),
  n_total = c(wald_design$n_total, score_design$n_total),
  n1 = c(wald_design$n1, score_design$n1),
  n2 = c(wald_design$n2, score_design$n2),
  total_events = round(c(wald_design$total_events, score_design$total_events), 1),
  variance_alt = round(c(wald_design$variance, score_design$variance), 4),
  variance_null = round(c(wald_design$variance_null, score_design$variance_null), 4)
)

design_comparison

## -----------------------------------------------------------------------------
analysis_times <- c(18, 24, 30)

gs_design <- gsNBCalendar(
  wald_design,
  k = 3,
  test.type = 4,
  beta = 1 - wald_design$power,
  analysis_times = analysis_times
)

data.frame(
  analysis = seq_along(gs_design$n.I),
  calendar_month = analysis_times,
  planned_information = round(gs_design$n.I, 2),
  information_fraction = round(gs_design$timing, 3)
)

## -----------------------------------------------------------------------------
set.seed(2026)

demo_enroll_rate <- data.frame(rate = 30 / 6, duration = 6)
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(design_args$lambda1, design_args$lambda2),
  dispersion = c(design_args$dispersion, design_args$dispersion)
)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(design_args$dropout_rate, design_args$dropout_rate),
  duration = c(100, 100)
)

sim_data <- nb_sim(
  enroll_rate = demo_enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = design_args$max_followup,
  n = 60,
  event_gap = design_args$event_gap
)

cut_data <- cut_data_by_date(
  sim_data,
  cut_date = 12,
  event_gap = design_args$event_gap
)

head(cut_data)

## -----------------------------------------------------------------------------
score_test <- mutze_test(cut_data, test_type = "score", sided = 1)
score_test

## ----eval=FALSE---------------------------------------------------------------
# production_enroll_rate <- data.frame(
#   rate = wald_design$accrual_rate,
#   duration = wald_design$accrual_duration
# )
# 
# set.seed(2026)
# sim_results <- sim_gs_nbinom(
#   n_sims = 10000,
#   enroll_rate = production_enroll_rate,
#   fail_rate = fail_rate,
#   dropout_rate = dropout_rate,
#   max_followup = design_args$max_followup,
#   event_gap = design_args$event_gap,
#   n_target = wald_design$n_total,
#   design = gs_design,
#   analysis_times = analysis_times,
#   test_type = "score",
#   seed = TRUE
# )
# 
# bounded <- check_gs_bound(
#   sim_results,
#   gs_design,
#   info_col = "info_unblinded_ml"
# )
# summarize_gs_sim(bounded)

