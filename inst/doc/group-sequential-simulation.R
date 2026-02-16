## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(gsDesign)
library(gsDesignNB)
library(data.table)
library(ggplot2)
library(gt)

## ----sample-size--------------------------------------------------------------
# Sample size calculation
# Enrollment: constant rate over 12 months
# Trial duration: 24 months
event_gap_val <- 20 / 30.4375 # Minimum gap between events is 20 days (approx)

nb_ss <- sample_size_nbinom(
  lambda1 = 1.5 / 12, # Control event rate (per month)
  lambda2 = 1.0 / 12, # Experimental event rate (per month)
  dispersion = 0.5, # Overdispersion parameter
  power = 0.9, # 90% power
  alpha = 0.025, # One-sided alpha
  accrual_rate = 1, # This will be scaled to achieve the target power
  accrual_duration = 12, # 12 months enrollment
  trial_duration = 24, # 24 months trial
  max_followup = 12, # 12 months of follow-up per patient
  dropout_rate = -log(0.95) / 12, # 5% dropout rate at 1 year
  event_gap = event_gap_val
)

# Print key results
message("Fixed design")
nb_ss

## ----gs-design----------------------------------------------------------------
# Analysis times (in months)
analysis_times <- c(10, 18, 24)

# Create group sequential design with integer sample sizes
gs_nb <- gsNBCalendar(
  x = nb_ss, # Input fixed design for negative binomial
  k = 3, # 3 analyses
  test.type = 4, # Two-sided asymmetric, non-binding futility
  sfu = sfLinear, # Linear spending function for upper bound
  sfupar = c(.5, .5), # Identity function
  sfl = sfHSD, # HSD spending for lower bound
  sflpar = -8, # Conservative futility bound
  usTime = c(.1, .18, 1), # Upper spending timing
  lsTime = NULL, # Spending based on information
  analysis_times = analysis_times # Calendar times in months
) |> gsDesignNB::toInteger() # Round to integer sample size

## -----------------------------------------------------------------------------
summary(gs_nb)

## -----------------------------------------------------------------------------
gs_nb |>
  gsDesign::gsBoundSummary(
    deltaname = "RR",
    logdelta = TRUE,
    Nname = "Information",
    timename = "Month",
    digits = 4,
    ddigits = 2
  ) |>
  gt() |>
  tab_header(
    title = "Group Sequential Design Bounds for Negative Binomial Outcome",
    subtitle = paste0(
      "N = ", ceiling(gs_nb$n_total[gs_nb$k]),
      ", Expected events = ", round(gs_nb$nb_design$total_events, 1)
    )
  )

## ----load-results-------------------------------------------------------------
# Load pre-computed simulation results
results_file <- system.file("extdata", "gs_simulation_results.rds", package = "gsDesignNB")

if (results_file == "" && file.exists("../inst/extdata/gs_simulation_results.rds")) {
  results_file <- "../inst/extdata/gs_simulation_results.rds"
}

if (results_file != "") {
  sim_data <- readRDS(results_file)
  results <- sim_data$results
  summary_gs <- sim_data$summary_gs
  n_sims <- sim_data$n_sims
  params <- sim_data$params
} else {
  # Fallback if data is not available (e.g. not installed yet)
  warning("Simulation results not found. Skipping simulation analysis.")
  results <- NULL
  summary_gs <- NULL
  n_sims <- 0
}

## ----summary-table, eval=!is.null(results), results='asis'--------------------
# Helper function for trimmed mean (to handle outliers in blinded info)
trimmed_mean <- function(x, trim = 0.01) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x, trim = trim)
}

# Create comprehensive theoretical vs simulation comparison table for each analysis
dt <- data.table::as.data.table(results)

# Function to compute comparison for a specific analysis
get_analysis_comparison <- function(analysis_num) {
  sub_dt <- dt[analysis == analysis_num]
  
  # Get theoretical values from gs_nb design
  theo_n <- gs_nb$n_total[analysis_num]
  theo_n1 <- gs_nb$n1[analysis_num]
  theo_n2 <- gs_nb$n2[analysis_num]
  theo_exposure <- gs_nb$exposure[analysis_num]
  theo_exposure_at_risk1 <- gs_nb$exposure_at_risk1[analysis_num]
  theo_exposure_at_risk2 <- gs_nb$exposure_at_risk2[analysis_num]
  theo_events1 <- gs_nb$events1[analysis_num]
  theo_events2 <- gs_nb$events2[analysis_num]
  theo_events <- gs_nb$events[analysis_num]
  theo_variance <- gs_nb$variance[analysis_num]
  theo_info <- gs_nb$n.I[analysis_num]
  
  # Observed values (using trimmed means for robustness)
  obs_n <- mean(sub_dt$n_enrolled)
  obs_n1 <- mean(sub_dt$n_ctrl)
  obs_n2 <- mean(sub_dt$n_exp)
  obs_exposure1 <- mean(sub_dt$exposure_total_ctrl)
  obs_exposure2 <- mean(sub_dt$exposure_total_exp)
  obs_exposure_at_risk1 <- mean(sub_dt$exposure_at_risk_ctrl)
  obs_exposure_at_risk2 <- mean(sub_dt$exposure_at_risk_exp)
  obs_events1 <- mean(sub_dt$events_ctrl)
  obs_events2 <- mean(sub_dt$events_exp)
  obs_events <- mean(sub_dt$events_total)
  obs_variance <- median(sub_dt$se^2, na.rm = TRUE)
  obs_info_blinded <- trimmed_mean(sub_dt$blinded_info, trim = 0.01)
  obs_info_unblinded <- trimmed_mean(sub_dt$unblinded_info, trim = 0.01)
  
  # Build comparison data frame
  data.frame(
    Metric = c(
      "N Enrolled",
      "N Control",
      "N Experimental",
      "Total Exposure - Control",
      "Total Exposure - Experimental",
      "Exposure at Risk - Control",
      "Exposure at Risk - Experimental",
      "Events - Control",
      "Events - Experimental",
      "Events - Total",
      "Variance of log(RR)",
      "Information (Blinded)",
      "Information (Unblinded)"
    ),
    Theoretical = c(
      theo_n, theo_n1, theo_n2,
      theo_n1 * theo_exposure, theo_n2 * theo_exposure,
      theo_n1 * theo_exposure_at_risk1, theo_n2 * theo_exposure_at_risk2,
      theo_events1, theo_events2, theo_events,
      theo_variance, theo_info, theo_info
    ),
    Simulated = c(
      obs_n, obs_n1, obs_n2,
      obs_exposure1, obs_exposure2,
      obs_exposure_at_risk1, obs_exposure_at_risk2,
      obs_events1, obs_events2, obs_events,
      obs_variance, obs_info_blinded, obs_info_unblinded
    ),
    stringsAsFactors = FALSE
  )
}

# Generate comparison table for each analysis
for (k in 1:3) {
  cat(sprintf("\n### Analysis %d (Month %d)\n\n", k, params$analysis_times[k]))
  
  comparison_k <- get_analysis_comparison(k)
  comparison_k$Difference <- comparison_k$Simulated - comparison_k$Theoretical
  comparison_k$Rel_Diff_Pct <- 100 * comparison_k$Difference / abs(comparison_k$Theoretical)
  
  print(
    comparison_k |>
      gt() |>
      tab_header(
        title = sprintf("Analysis %d: Theoretical vs Simulated", k),
        subtitle = sprintf("Calendar time = %d months", params$analysis_times[k])
      ) |>
      fmt_number(columns = c(Theoretical, Simulated, Difference), decimals = 2) |>
      fmt_number(columns = Rel_Diff_Pct, decimals = 1) |>
      cols_label(
        Metric = "",
        Theoretical = "Theoretical",
        Simulated = "Simulated",
        Difference = "Difference",
        Rel_Diff_Pct = "Rel. Diff (%)"
      ) |>
      tab_row_group(label = md("**Information**"), rows = grepl("Information|Variance", Metric)) |>
      tab_row_group(label = md("**Events**"), rows = grepl("Events", Metric)) |>
      tab_row_group(label = md("**Exposure**"), rows = grepl("Exposure", Metric)) |>
      tab_row_group(label = md("**Sample Size**"), rows = grepl("^N ", Metric)) |>
      row_group_order(groups = c("**Sample Size**", "**Exposure**", "**Events**", "**Information**"))
  )
  
  # Add boundary crossing summary
  sub_dt <- dt[analysis == k]
  cat(sprintf("\n**Boundary Crossing:**\n"))
  cat(sprintf("- Efficacy (upper): %.1f%% (n=%d)\n", 
              mean(sub_dt$cross_upper) * 100, sum(sub_dt$cross_upper)))
  cat(sprintf("- Futility (lower): %.1f%% (n=%d)\n", 
              mean(sub_dt$cross_lower) * 100, sum(sub_dt$cross_lower)))
  cat(sprintf("- Cumulative Efficacy: %.1f%%\n\n", 
              sum(dt[analysis <= k]$cross_upper) / n_sims * 100))
}

## ----overall-summary, eval=!is.null(results)----------------------------------
cat("=== Overall Operating Characteristics ===\n")
cat(sprintf("Number of simulations: %d\n", n_sims))
cat(sprintf("Overall Power (P[reject H0]): %.1f%% (SE: %.1f%%)\n", 
            summary_gs$power * 100, 
            sqrt(summary_gs$power * (1 - summary_gs$power) / n_sims) * 100))
cat(sprintf("Futility Stopping Rate: %.1f%%\n", summary_gs$futility * 100))
cat(sprintf("Design Power (target): %.1f%%\n", (1 - gs_nb$beta) * 100))

## ----power-comparison, eval=!is.null(results)---------------------------------
# Create comparison table
crossing_summary <- data.frame(
  Analysis = 1:3,
  Analysis_Time = params$analysis_times,
  Sim_Power = summary_gs$analysis_summary$prob_cross_upper,
  Sim_Cum_Power = summary_gs$analysis_summary$cum_prob_upper,
  Design_Cum_Power = cumsum(gs_nb$upper$prob[, 2])
)

crossing_summary |>
  gt() |>
  tab_header(
    title = "Power Comparison: Simulation vs Design",
    subtitle = sprintf("Based on %d simulated trials", n_sims)
  ) |>
  cols_label(
    Analysis = "Analysis",
    Analysis_Time = "Month",
    Sim_Power = "Incremental Power (Sim)",
    Sim_Cum_Power = "Cumulative Power (Sim)",
    Design_Cum_Power = "Cumulative Power (Design)"
  ) |>
  fmt_percent(columns = c(Sim_Power, Sim_Cum_Power, Design_Cum_Power), decimals = 1)

## ----plot-z-stats, fig.width=7, fig.height=5, fig.alt="Z-statistics across analyses with group sequential boundaries", eval=!is.null(results)----
# Prepare data for plotting
plot_data <- results
plot_data$z_flipped <- -plot_data$z_stat # Flip for efficacy direction

# Boundary data
bounds_df <- data.frame(
  analysis = 1:gs_nb$k,
  upper = gs_nb$upper$bound,
  lower = gs_nb$lower$bound
)

ggplot(plot_data, aes(x = factor(analysis), y = z_flipped)) +
  geom_violin(fill = "steelblue", alpha = 0.5, color = "steelblue") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  # Draw bounds as lines connecting analyses
  geom_line(
    data = bounds_df, aes(x = analysis, y = upper, group = 1),
    linetype = "dashed", color = "darkgreen", linewidth = 1
  ) +
  geom_line(
    data = bounds_df, aes(x = analysis, y = lower, group = 1),
    linetype = "dashed", color = "darkred", linewidth = 1
  ) +
  # Draw points for bounds
  geom_point(data = bounds_df, aes(x = analysis, y = upper), color = "darkgreen") +
  geom_point(data = bounds_df, aes(x = analysis, y = lower), color = "darkred") +
  geom_hline(yintercept = 0, color = "gray50") +
  labs(
    title = "Simulated Z-Statistics by Analysis",
    subtitle = "Green dashed = efficacy bound, Red dashed = futility bound",
    x = "Analysis",
    y = "Z-statistic (positive = favors experimental)"
  ) +
  theme_minimal() +
  ylim(c(-4, 6))

