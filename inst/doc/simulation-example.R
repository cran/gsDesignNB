## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(gsDesignNB)
library(data.table)
library(ggplot2)
library(lt)
library(scales)

## -----------------------------------------------------------------------------
# Enrollment: rate of 4 patients/month for 5 months -> ~20 patients
enroll_rate <- data.frame(
  rate = 4,
  duration = 5
)

# Failure rates (events per unit time)
fail_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.5, 0.3) # events per year (assuming time unit is year, adjust enroll duration if needed)
)

# Let's ensure time units are consistent.
# If fail_rate is per year, then durations should be in years.
# 5 months = 5/12 years.
enroll_rate <- data.frame(
  rate = 20 / (5 / 12), # 20 patients over 5/12 years
  duration = 5 / 12
)

# Dropout rates (per year)
dropout_rate <- data.frame(
  treatment = c("Control", "Experimental"),
  rate = c(0.1, 0.05),
  duration = c(100, 100) # constant rate for long duration
)

# Maximum follow-up per patient (years)
max_followup <- 2

## -----------------------------------------------------------------------------
set.seed(123)

sim_data <- nb_sim(
  enroll_rate = enroll_rate,
  fail_rate = fail_rate,
  dropout_rate = dropout_rate,
  max_followup = max_followup,
  n = 20
)

head(sim_data)

## -----------------------------------------------------------------------------
sim_dt <- as.data.table(sim_data)
sim_dt[, censor_followup := ifelse(event == 0, tte, 0)]
summary_stats <- sim_dt[
  ,
  .(
    n_subjects = uniqueN(id),
    total_events = sum(event == 1),
    total_followup = sum(censor_followup),
    observed_rate = sum(event == 1) / sum(censor_followup)
  ),
  by = treatment
]
summary_stats |>
  lt() |>
  lt_header(title = "Summary Statistics by Treatment") |>
  lt_label(
    treatment = "Treatment",
    n_subjects = "N",
    total_events = "Events",
    total_followup = "Follow-up",
    observed_rate = "Rate"
  ) |>
  lt_format(columns = "total_followup", decimals = 2) |>
  lt_format(columns = "observed_rate", decimals = 3)

## -----------------------------------------------------------------------------
head(sim_data, 10)

## ----fig.width=7, fig.height=5, fig.alt="Patient timelines with events (dots), 20-day gaps (gray segments), and censoring (X)"----
sim_plot <- as.data.frame(sim_data)
names(sim_plot) <- make.names(names(sim_plot), unique = TRUE)
events_df <- sim_plot[sim_plot$event == 1, ]
censor_df <- sim_plot[sim_plot$event == 0, ]

# Define a 20-day gap for visualization
gap_duration <- 20 / 365.25

# Create segments for the gap after each event
events_df$gap_start <- events_df$tte
events_df$gap_end <- events_df$tte + gap_duration

ggplot(sim_plot, aes(x = tte, y = factor(id), color = treatment)) +
  geom_line(aes(group = id), color = "gray80") +
  # Add gap segments
  geom_segment(
    data = events_df,
    aes(x = gap_start, xend = gap_end, y = factor(id), yend = factor(id)),
    color = "gray50", linewidth = 2, alpha = 0.7
  ) +
  geom_point(data = events_df, shape = 19, size = 2) +
  geom_point(data = censor_df, shape = 4, size = 3) +
  labs(
    title = "Patient Timelines",
    x = "Time from Randomization (Years)",
    y = "Patient ID",
    caption = "Dots = Events, Gray Bars = 20-day Gap, X = Censoring/Dropout"
  ) +
  theme_minimal()

## -----------------------------------------------------------------------------
cut_summary <- cut_data_by_date(sim_data, cut_date = 1.5)
head(cut_summary)

## -----------------------------------------------------------------------------
mutze_res <- mutze_test(cut_summary)
mutze_res

## -----------------------------------------------------------------------------
# Target 15 total events
target_events <- 15
analysis_date <- get_analysis_date(sim_data, planned_events = target_events)

print(paste("Calendar date for", target_events, "events:", round(analysis_date, 3)))

# Cut data at this date
cut_events <- cut_data_by_date(sim_data, cut_date = analysis_date)

# Verify event count
sum(cut_events$events)

## -----------------------------------------------------------------------------
# Define failure rates with dispersion
fail_rate_nb <- data.frame(
  treatment = "Control",
  rate = 10, # Mean event rate
  dispersion = 2 # Variance = mean + 2 * mean^2
)

enroll_rate_nb <- data.frame(
  rate = 100,
  duration = 1
)

set.seed(1)
# Simulate 10000 subjects for a quick local review build
sim_nb <- nb_sim(
  enroll_rate = enroll_rate_nb,
  fail_rate = fail_rate_nb,
  max_followup = 1,
  n = 10000,
  block = "Control" # Assign all to Control for simplicity
)

## -----------------------------------------------------------------------------
# Count events per subject
counts_nb <- as.data.table(sim_nb)[, .(events = sum(event)), by = id]

m <- mean(counts_nb$events)
v <- var(counts_nb$events)
k_mom <- (v - m) / (m^2)

# Theoretical values
mu_true <- 10
k_true <- 2
v_true <- mu_true + k_true * mu_true^2

# Also estimate using GLM
# We use MASS::glm.nb to fit the negative binomial model
# We suppress warnings because fitting intercept-only models on simulated data
# can occasionally produce convergence warnings despite valid estimates.
k_glm <- tryCatch(
  {
    fit <- suppressWarnings(MASS::glm.nb(events ~ 1, data = counts_nb))
    1 / fit$theta
  },
  error = function(e) NA
)

print(paste("True Mean:", mu_true, "| Observed Mean:", signif(m, 4)))
print(paste("True Variance:", v_true, "| Observed Variance:", signif(v, 4)))
print(paste("True Dispersion:", k_true))
print(paste("Estimated Dispersion (MoM):", signif(k_mom, 4)))
print(paste("Estimated Dispersion (GLM):", signif(k_glm, 4)))

## -----------------------------------------------------------------------------
# Calculate observed proportions
obs_dist <- counts_nb[, .N, by = events]
obs_dist[, prop := N / sum(N)]
obs_dist[, type := "Observed"]

# Calculate theoretical probabilities
# Parameters: mu = 10, size = 1/k = 1/2 = 0.5
mu <- 10
k <- 2
size <- 1/k
max_events <- max(obs_dist$events)
theo_probs <- dnbinom(0:max_events, size = size, mu = mu)
theo_dist <- data.table(
  events = 0:max_events,
  N = NA,
  prop = theo_probs,
  type = "Theoretical"
)

# Combine for plotting
plot_data <- rbind(obs_dist, theo_dist)

# Filter for visualization (limit x-axis to 50)
plot_data <- plot_data[events <= 50]

# Plot
ggplot(plot_data, aes(x = events, y = prop, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  scale_y_log10(labels = scales::label_number()) +
  labs(
    title = "Observed vs. Theoretical Negative Binomial Distribution",
    subtitle = paste("Mean =", mu, ", Dispersion =", k, "(Log Scale)"),
    x = "Number of Events",
    y = "Proportion",
    fill = "Distribution"
  ) +
  theme_minimal()

