## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(gsDesignNB)

## -----------------------------------------------------------------------------
# All rates and durations in years
lambda1         <- 3.0                  # control exacerbation rate (events/year)
lambda2         <- 3.0                  # experimental rate (no difference under H1)
dispersion      <- 0.625                # NB dispersion k
rr0             <- 1.15                 # non-inferiority margin
event_gap       <- 0                    # no gap specified in protocol
max_followup    <- 1                    # 1-year follow-up
dropout_retained <- 0.70                # 30% dropout/violator -> inflate by 1/0.70
alpha           <- 0.025
power           <- 0.95

# Arbitrary accrual (does not affect N when all patients get the same max FU)
accrual_rate    <- 100                  # subjects/year (placeholder)
accrual_dur     <- 1                    # 1-year enrollment (placeholder)
trial_duration  <- 2                    # years

## -----------------------------------------------------------------------------
wald_evaluable <- sample_size_nbinom(
  lambda1          = lambda1,
  lambda2          = lambda2,
  dispersion       = dispersion,
  power            = power,
  alpha            = alpha,
  rr0              = rr0,
  accrual_rate     = accrual_rate,
  accrual_duration = accrual_dur,
  trial_duration   = trial_duration,
  dropout_rate     = 0,
  max_followup     = max_followup,
  event_gap        = event_gap,
  test_type        = "wald"
)

n_inflated <- ceiling(ceiling(wald_evaluable$n_total) / dropout_retained)

## -----------------------------------------------------------------------------
wald_exposure <- sample_size_nbinom(
  lambda1          = lambda1,
  lambda2          = lambda2,
  dispersion       = dispersion,
  power            = power,
  alpha            = alpha,
  rr0              = rr0,
  accrual_rate     = accrual_rate,
  accrual_duration = accrual_dur,
  trial_duration   = trial_duration,
  dropout_rate     = -log(dropout_retained) / max_followup,
  max_followup     = max_followup,
  event_gap        = event_gap,
  test_type        = "wald"
)

data.frame(
  Method           = c("Full follow-up + inflate / 0.70",
                        "Dropout via exposure integral",
                        "Protocol target"),
  N                = c(n_inflated, ceiling(wald_exposure$n_total), 3332)
) |> knitr::kable(caption = "Sample size comparison: two dropout approaches vs protocol.")

## -----------------------------------------------------------------------------
summary(wald_evaluable)

## -----------------------------------------------------------------------------
score_evaluable <- sample_size_nbinom(
  lambda1          = lambda1,
  lambda2          = lambda2,
  dispersion       = dispersion,
  power            = power,
  alpha            = alpha,
  rr0              = rr0,
  accrual_rate     = accrual_rate,
  accrual_duration = accrual_dur,
  trial_duration   = trial_duration,
  dropout_rate     = 0,
  max_followup     = max_followup,
  event_gap        = event_gap,
  test_type        = "score"
)

score_inflated <- ceiling(ceiling(score_evaluable$n_total) / dropout_retained)

summary(score_evaluable)

## -----------------------------------------------------------------------------
comparison <- data.frame(
  Sizing        = c("Wald", "Score", "Protocol target"),
  evaluable_N   = c(ceiling(wald_evaluable$n_total), ceiling(score_evaluable$n_total), NA),
  randomized_N  = c(n_inflated, score_inflated, 3332),
  events        = round(c(wald_evaluable$total_events, score_evaluable$total_events, NA), 1)
)

knitr::kable(comparison,
  caption = "Wald vs score sizing with randomized-sample inflation."
)

## -----------------------------------------------------------------------------
wald_gap_evaluable <- sample_size_nbinom(
  lambda1          = lambda1,
  lambda2          = lambda2,
  dispersion       = dispersion,
  power            = power,
  alpha            = alpha,
  rr0              = rr0,
  accrual_rate     = accrual_rate,
  accrual_duration = accrual_dur,
  trial_duration   = trial_duration,
  dropout_rate     = 0,
  max_followup     = max_followup,
  event_gap        = 28 / 365.25,       # 28 days in years
  test_type        = "wald"
)

n_gap_inflated <- ceiling(ceiling(wald_gap_evaluable$n_total) / dropout_retained)

data.frame(
  Formula   = c("No gap (FLAME protocol)", "28-day gap (Jensen-corrected)"),
  N         = c(n_inflated, n_gap_inflated),
  delta_n   = c(0, n_gap_inflated - n_inflated)
) |> knitr::kable(
  caption = "Impact of a hypothetical 28-day event gap on the FLAME sample size."
)

## ----eval=FALSE---------------------------------------------------------------
# n_sims <- 5000
# n_per_arm <- ceiling(n_inflated / 2)
# 
# enroll_rate <- data.frame(rate = n_per_arm * 2, duration = accrual_dur)
# fail_rate <- data.frame(
#   treatment  = c("Control", "Experimental"),
#   rate       = c(lambda1, lambda2),
#   dispersion = c(dispersion, dispersion)
# )
# dropout_df <- data.frame(
#   treatment = c("Control", "Experimental"),
#   rate      = c(-log(dropout_retained), -log(dropout_retained)),
#   duration  = c(100, 100)
# )
# 
# set.seed(20260506)
# rejections <- 0L
# for (i in seq_len(n_sims)) {
#   dat <- nb_sim(
#     enroll_rate  = enroll_rate,
#     fail_rate    = fail_rate,
#     dropout_rate = dropout_df,
#     max_followup = max_followup
#   )
#   dat <- cut_data_by_date(dat, trial_duration)
#   tst <- mutze_test(dat, test_type = "wald")
#   # NI z-statistic: shift from H0: beta=0 to H0: rr >= rr0
#   z_ni <- (log(rr0) - tst$estimate) / tst$se
#   if (z_ni >= qnorm(1 - alpha)) rejections <- rejections + 1L
# }
# 
# empirical_power <- rejections / n_sims
# 
# data.frame(
#   Metric = c("Target power", "Empirical power", "Replicates"),
#   Value  = c(sprintf("%.1f%%", power * 100),
#              sprintf("%.1f%%", empirical_power * 100),
#              n_sims)
# ) |> knitr::kable(caption = "Fixed-design simulation verification.")

## -----------------------------------------------------------------------------
# IMPACT: FF/UMEC/VI vs UMEC/VI
lambda_triple <- 0.80                     # FF/UMEC/VI rate (events/year)
lambda_umec   <- lambda_triple / (1 - 0.15)  # UMEC/VI rate (~0.941)

impact_wald <- sample_size_nbinom(
  lambda1          = lambda_umec,           # control = UMEC/VI
  lambda2          = lambda_triple,         # experimental = FF/UMEC/VI
  dispersion       = 0.75,
  power            = 0.90,
  alpha            = 0.005,                 # two-sided 1% → one-sided 0.5%
  sided            = 1,
  ratio            = 2,                     # 2:1 (experimental:control)
  accrual_rate     = 100,
  accrual_duration = 1,
  trial_duration   = 2,
  max_followup     = 1,                     # 52-week follow-up
  event_gap        = 28 / 365.25,           # 28-day gap
  test_type        = "wald"
)

data.frame(
  Metric = c("N control (UMEC/VI)", "N experimental (FF/UMEC/VI)",
             "N total (this comparison)", "Protocol target (comparison)"),
  Value  = c(ceiling(impact_wald$n1), ceiling(impact_wald$n2),
             ceiling(impact_wald$n_total), "2,000 + 4,000 = 6,000")
) |> knitr::kable(caption = "IMPACT: FF/UMEC/VI vs UMEC/VI.")

## -----------------------------------------------------------------------------
# IMPACT: FF/UMEC/VI vs FF/VI
lambda_ffvi <- lambda_triple / (1 - 0.12)  # FF/VI rate (~0.909)

impact_ffvi <- sample_size_nbinom(
  lambda1          = lambda_ffvi,
  lambda2          = lambda_triple,
  dispersion       = 0.75,
  power            = 0.90,
  alpha            = 0.005,
  sided            = 1,
  ratio            = 1,                     # 1:1 (4,000 : 4,000)
  accrual_rate     = 100,
  accrual_duration = 1,
  trial_duration   = 2,
  max_followup     = 1,
  event_gap        = 28 / 365.25,
  test_type        = "wald"
)

data.frame(
  Metric = c("N per arm (FF/UMEC/VI vs FF/VI)", "N total", "Protocol target"),
  Value  = c(ceiling(impact_ffvi$n1), ceiling(impact_ffvi$n_total), "4,000 + 4,000 = 8,000")
) |> knitr::kable(caption = "IMPACT: FF/UMEC/VI vs FF/VI (12% reduction).")

## ----message=FALSE, warning=FALSE, fig.width=6, fig.height=4------------------
library(ggplot2)

k_vals <- seq(0.3, 1.5, by = 0.05)
n_by_k <- vapply(k_vals, function(k) {
  res <- sample_size_nbinom(
    lambda1 = lambda_umec, lambda2 = lambda_triple, dispersion = k,
    power = 0.90, alpha = 0.005, sided = 1, ratio = 2,
    accrual_rate = 100, accrual_duration = 1, trial_duration = 2,
    max_followup = 1, event_gap = 28 / 365.25, test_type = "wald"
  )
  ceiling(res$n_total)
}, numeric(1))

sens_df <- data.frame(k = k_vals, N = n_by_k)

p <- ggplot(sens_df, aes(x = k, y = N)) +
  geom_line(linewidth = 0.8, colour = "#2166AC") +
  geom_point(
    data = sens_df[sens_df$k == 0.75, ],
    colour = "#B2182B", size = 3
  ) +
  annotate("text", x = 0.75, y = sens_df$N[sens_df$k == 0.75],
           label = "Protocol assumption", vjust = -1.2,
           colour = "#B2182B", size = 3.5) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Dispersion parameter k",
    y = "Total sample size N",
    title = "IMPACT: sample size sensitivity to dispersion"
  ) +
  theme_minimal(base_size = 12)

p

## -----------------------------------------------------------------------------
impact_nogap <- sample_size_nbinom(
  lambda1 = lambda_umec, lambda2 = lambda_triple, dispersion = 0.75,
  power = 0.90, alpha = 0.005, sided = 1, ratio = 2,
  accrual_rate = 100, accrual_duration = 1, trial_duration = 2,
  max_followup = 1, event_gap = 0, test_type = "wald"
)

data.frame(
  Design      = c("No event gap", "28-day gap (Jensen-corrected)"),
  N           = c(ceiling(impact_nogap$n_total), ceiling(impact_wald$n_total)),
  Delta       = c(0, ceiling(impact_wald$n_total) - ceiling(impact_nogap$n_total))
) |> knitr::kable(caption = "Impact of the 28-day event gap on IMPACT sample size.")

## -----------------------------------------------------------------------------
# OPERA: ocrelizumab vs Rebif
opera_k_vals <- c(0.5, 0.75, 1.0, 1.5)
opera_results <- lapply(opera_k_vals, function(k) {
  res <- sample_size_nbinom(
    lambda1          = 0.33,                 # Rebif ARR (events/year)
    lambda2          = 0.165,                # ocrelizumab ARR
    dispersion       = k,
    power            = 0.84,
    alpha            = 0.025,                # two-sided 5% → one-sided 2.5%
    sided            = 1,
    ratio            = 1,
    accrual_rate     = 100,
    accrual_duration = 1,
    trial_duration   = 3,
    max_followup     = 96 / 52,              # 96 weeks ≈ 1.846 years
    dropout_rate     = -log(0.80) / (96/52), # 20% dropout over 96 weeks
    test_type        = "wald"
  )
  data.frame(k = k, N_per_arm = ceiling(res$n1), N_total = ceiling(res$n_total))
})

opera_df <- do.call(rbind, opera_results)
opera_df$Protocol <- 400

knitr::kable(opera_df,
  caption = "OPERA sample size for varying dispersion k (84% power, RR = 0.50)."
)

## -----------------------------------------------------------------------------
friede_scenarios <- expand.grid(
  theta = c(0.7, 0.8),
  power = c(0.80, 0.90)
)
lambda_overall <- 1.5
phi <- 0.5

friede_results <- mapply(function(theta, pwr) {
  lam0 <- 2 * lambda_overall / (1 + theta)
  lam1 <- theta * lam0
  res <- sample_size_nbinom(
    lambda1          = lam0,
    lambda2          = lam1,
    dispersion       = phi,
    power            = pwr,
    alpha            = 0.025,
    sided            = 1,
    ratio            = 1,
    accrual_rate     = 100,
    accrual_duration = 1,
    trial_duration   = 2,
    max_followup     = 1,               # ti = 1 year (equal FU for all)
    dropout_rate     = 0,
    test_type        = "wald"
  )
  ceiling(res$n1)
}, friede_scenarios$theta, friede_scenarios$power)

published_n0 <- c(147, 370, 196, 496)

knitr::kable(
  data.frame(
    theta     = friede_scenarios$theta,
    power     = friede_scenarios$power,
    published = published_n0,
    gsDesignNB = friede_results,
    difference = friede_results - published_n0,
    within_one = ifelse(abs(friede_results - published_n0) <= 1, "Yes", "No")
  ),
  caption = "Friede & Schmidli (2010) Table 1: fixed-design validation."
)

## -----------------------------------------------------------------------------
friede_ni <- sample_size_nbinom(
  lambda1          = 1.16,                # active control rate
  lambda2          = 1.16,                # experimental (no difference)
  dispersion       = 0.46,
  power            = 0.80,
  alpha            = 0.025,
  sided            = 1,
  ratio            = 1,
  rr0              = 1.15,               # NI margin
  accrual_rate     = 100,
  accrual_duration = 1,
  trial_duration   = 2,
  max_followup     = 1,
  dropout_rate     = 0,
  test_type        = "wald"
)

summary(friede_ni)

## -----------------------------------------------------------------------------
# Mutze Table 4: HF scenarios, κ=2 (k=0.5), θ=0.70
# Their λ₁=λ₂ under H₀ gives the "overall" rate; under H₁ λ₂=θ·λ₁
# We verify a subset of fixed-design n₁ values
mutze_hf <- expand.grid(
  lambda_annual = c(0.08, 0.10, 0.12, 0.14),
  theta         = c(0.70, 0.80),
  kappa         = 2,
  power         = c(0.80, 0.90)
)

mutze_hf_results <- mapply(function(lam, theta, kappa, pwr) {
  k <- 1 / kappa
  # Rates are annualized; follow-up is variable (15-mo enrollment, 48-mo study)
  # Convert to monthly for enrollment/trial_duration consistency
  lam_mo <- lam / 12
  res <- sample_size_nbinom(
    lambda1          = lam_mo,              # control rate (monthly)
    lambda2          = theta * lam_mo,      # experimental rate
    dispersion       = k,
    power            = pwr,
    alpha            = 0.025,
    sided            = 1,
    ratio            = 1,
    accrual_rate     = 100,
    accrual_duration = 15,                  # 15 months enrollment
    trial_duration   = 48,                  # 48 months study
    dropout_rate     = 0,
    test_type        = "wald"
  )
  ceiling(res$n1)
}, mutze_hf$lambda_annual, mutze_hf$theta, mutze_hf$kappa, mutze_hf$power)

mutze_hf$k <- 1 / mutze_hf$kappa
mutze_hf$n1_gsDesignNB <- mutze_hf_results

knitr::kable(
  mutze_hf[, c("lambda_annual", "theta", "k", "power", "n1_gsDesignNB")],
  caption = "Mutze et al. (2019): HF fixed-design n₁ (κ=2, k=0.5)."
)

## -----------------------------------------------------------------------------
# Mutze Table 5: MS scenarios, κ=2 (k=0.5), fixed 6-month FU
mutze_ms <- expand.grid(
  lambda_annual = c(6, 8, 10),   # annualized monthly CUAL rates
  theta         = c(0.50, 0.70),
  kappa         = 2,
  power         = c(0.80, 0.90)
)

mutze_ms_results <- mapply(function(lam, theta, kappa, pwr) {
  k <- 1 / kappa
  lam_mo <- lam / 12          # monthly rate
  res <- sample_size_nbinom(
    lambda1          = lam_mo,
    lambda2          = theta * lam_mo,
    dispersion       = k,
    power            = pwr,
    alpha            = 0.025,
    sided            = 1,
    ratio            = 1,
    accrual_rate     = 100,
    accrual_duration = 1,
    trial_duration   = 7,       # 6-month FU + enrollment
    max_followup     = 6,       # fixed 6 months
    dropout_rate     = 0,
    test_type        = "wald"
  )
  ceiling(res$n1)
}, mutze_ms$lambda_annual, mutze_ms$theta, mutze_ms$kappa, mutze_ms$power)

mutze_ms$k <- 1 / mutze_ms$kappa
mutze_ms$n1_gsDesignNB <- mutze_ms_results

knitr::kable(
  mutze_ms[, c("lambda_annual", "theta", "k", "power", "n1_gsDesignNB")],
  caption = "Mutze et al. (2019): MS fixed-design n₁ (κ=2, k=0.5)."
)

## -----------------------------------------------------------------------------
# Estimate rates from CHARM-Preserved
rate_placebo   <- 547 / 4374     # events per person-year
rate_treatment <- 392 / 4425
rr_observed    <- rate_treatment / rate_placebo

charm <- sample_size_nbinom(
  lambda1          = rate_placebo,
  lambda2          = rate_treatment,
  dispersion       = 0.5,           # typical HF dispersion
  power            = 0.90,
  alpha            = 0.025,
  sided            = 1,
  ratio            = 1,
  accrual_rate     = 100,
  accrual_duration = 12,
  trial_duration   = 48,            # ~4-year study
  max_followup     = 36,            # ~3-year median FU
  dropout_rate     = 0,
  test_type        = "wald"
)

data.frame(
  Metric = c("Control rate (events/person-year)", "Rate ratio",
             "Required N (90% power)", "CHARM-Preserved actual N"),
  Value  = c(round(rate_placebo, 3), round(rr_observed, 3),
             ceiling(charm$n_total), 3023)
) |> knitr::kable(caption = "Design from CHARM-Preserved observed data.")

## -----------------------------------------------------------------------------
bold_wald <- sample_size_nbinom(
  lambda1          = 1.39,                # monthly CUAL rate, placebo
  lambda2          = 0.42,                # siponimod 2mg
  dispersion       = 0.5,                # assumed
  power            = 0.80,
  alpha            = 0.025,
  sided            = 1,
  ratio            = 1,
  accrual_rate     = 100,
  accrual_duration = 1,
  trial_duration   = 4,
  max_followup     = 3,                  # 3 months
  dropout_rate     = 0,
  test_type        = "wald"
)

bold_score <- sample_size_nbinom(
  lambda1          = 1.39,
  lambda2          = 0.42,
  dispersion       = 0.5,
  power            = 0.80,
  alpha            = 0.025,
  sided            = 1,
  ratio            = 1,
  accrual_rate     = 100,
  accrual_duration = 1,
  trial_duration   = 4,
  max_followup     = 3,
  dropout_rate     = 0,
  test_type        = "score"
)

data.frame(
  Test  = c("Wald", "Score"),
  N_per_arm = c(ceiling(bold_wald$n1), ceiling(bold_score$n1)),
  N_total   = c(ceiling(bold_wald$n_total), ceiling(bold_score$n_total))
) |> knitr::kable(
  caption = "BOLD trial sizing: Wald vs Score (RR = 0.30, high rate)."
)

## -----------------------------------------------------------------------------
paradigm_nb <- sample_size_nbinom(
  lambda1          = 0.3,                  # HF hospitalization rate, monthly
  lambda2          = 0.3 * 0.80,           # 20% reduction
  dispersion       = 0.5,
  power            = 0.90,
  alpha            = 0.025,
  sided            = 1,
  ratio            = 1,
  accrual_rate     = 362,                  # ~7,980 over 22 months
  accrual_duration = 22,
  trial_duration   = 43,                   # 22 enrollment + 21 minimum FU
  dropout_rate     = 0,
  test_type        = "wald"
)

data.frame(
  Design      = c("NB recurrent events (hypothetical)",
                   "Cox time-to-first (protocol)"),
  N           = c(ceiling(paradigm_nb$n_total), 7980),
  Endpoint    = c("Recurrent HF hospitalizations", "CV death or first HF hosp.")
) |> knitr::kable(
  caption = "PARADIGM-HF: NB recurrent-event sizing vs actual protocol."
)

