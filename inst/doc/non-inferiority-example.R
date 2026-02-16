## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(gsDesignNB)

## ----ni-example---------------------------------------------------------------
ni_design <- sample_size_nbinom(
  lambda1 = 0.1,
  lambda2 = 0.09,
  rr0 = 1.1,
  dispersion = 0.5,
  power = 0.9,
  alpha = 0.025,
  trial_duration = 2,
  accrual_duration = 1,
  accrual_rate = 100 # Dummy value to get sample size
)

print(ni_design)

## ----super-sup-example--------------------------------------------------------
ss_design <- sample_size_nbinom(
  lambda1 = 0.1,
  lambda2 = 0.02,
  rr0 = 0.5,
  dispersion = 0.5,
  power = 0.9,
  alpha = 0.025,
  trial_duration = 2,
  accrual_duration = 1, # Assumed accrual duration
  accrual_rate = 100,
  max_followup = 1
)

print(ss_design)

