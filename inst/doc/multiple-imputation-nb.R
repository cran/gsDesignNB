## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 7,
  fig.height = 4
)

## ----simulate-----------------------------------------------------------------
set.seed(2025)

n_subj   <- 90L   # 30 per arm
n_visit  <- 4L
k_true   <- 0.5   # dispersion
lambda   <- c(placebo = 1.8, low = 1.1, high = 0.7)   # mean counts/visit

trt_labels <- rep(c("placebo", "low", "high"), each = n_subj / 3L)
strat      <- sample(c(0L, 1L), n_subj, replace = TRUE)
baseline   <- rnbinom(n_subj, mu = 2.5, size = 1 / k_true)
id         <- seq_len(n_subj)

# Build long-format data (all visits complete at first)
long_data <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
  mu_i <- lambda[trt_labels[i]] * exp(0.15 * strat[i] + 0.05 * baseline[i])
  data.frame(
    id       = i,
    visit    = seq_len(n_visit),
    trt      = trt_labels[i],
    strat1   = strat[i],
    baseline = baseline[i],
    count    = rnbinom(n_visit, mu = mu_i, size = 1 / k_true),
    stringsAsFactors = FALSE
  )
}))

# Sort
long_data <- long_data[order(long_data$id, long_data$visit), ]
rownames(long_data) <- NULL
str(long_data)

## ----introduce_missing--------------------------------------------------------
long_data$miss_flag <- NA_character_

# --- MAR: ~15% of subjects drop out from visit 3 onward ---
mar_subjects <- sample(id, size = round(0.15 * n_subj))
long_data$miss_flag[
  long_data$id %in% mar_subjects & long_data$visit >= 3
] <- "MAR"

# --- MNAR: ~10% of active-arm subjects withdraw at visit 2 ---
active_ids  <- id[trt_labels != "placebo"]
mnar_subjects <- sample(active_ids, size = round(0.10 * n_subj))
long_data$miss_flag[
  long_data$id %in% mnar_subjects & long_data$visit >= 2
] <- "MNAR"

# --- Composite ICE: ~5% of all subjects, post-event visits ---
comp_subjects <- sample(setdiff(id, c(mar_subjects, mnar_subjects)),
                        size = round(0.05 * n_subj))
comp_ice_visit <- sample(2L:4L, length(comp_subjects), replace = TRUE)
for (ci in seq_along(comp_subjects)) {
  long_data$miss_flag[
    long_data$id == comp_subjects[ci] &
    long_data$visit >= comp_ice_visit[ci]
  ] <- "Comp"
}

# Set the outcome to NA for flagged rows (simulate non-collection)
long_data$count[!is.na(long_data$miss_flag)] <- NA_integer_

cat("Missing rows by mechanism:\n")
print(table(long_data$miss_flag, useNA = "ifany"))

## ----run_mi, eval=requireNamespace("glmmTMB", quietly = TRUE)-----------------
# library(gsDesignNB)
# 
# # The formula mirrors the fixed-effects structure from the SAS code.
# # A random intercept per subject captures between-subject heterogeneity.
# mi_formula <- count ~ baseline + strat1 + trt + visit + (1 | id)
# 
# imp_result <- impute_nb(
#   data            = long_data,
#   formula         = mi_formula,
#   outcome_col     = "count",
#   miss_flag_col   = "miss_flag",
#   baseline_col    = "baseline",
#   trt_col         = "trt",
#   reference_trt   = "placebo",
#   subject_col     = "id",
#   strata_cols     = c("trt", "strat1"),
#   mar_values      = "MAR",
#   mnar_value      = "MNAR",
#   composite_value = "Comp",
#   n_imp           = 5L,
#   n_boot          = 1L,
#   seed            = 42L
# )
# 
# # Dimensions: n_subj * n_visit * n_imp rows
# dim(imp_result)
# names(imp_result)

## ----run_mi_fallback, echo=FALSE, eval=!requireNamespace("glmmTMB", quietly = TRUE)----
message(
  "NOTE: 'glmmTMB' is not installed in this build environment.\n",
  "Install it with install.packages('glmmTMB') to run this vignette."
)

## ----inspect, eval=requireNamespace("glmmTMB", quietly = TRUE)----------------
# # Show imputed rows only, first 2 imputations
# imp_missing <- imp_result[
#   !is.na(imp_result$miss_flag) & imp_result$imputation <= 2,
#   c("id", "visit", "trt", "miss_flag", "baseline", "count", "imputation",
#     "imputed_value")
# ]
# head(imp_missing[order(imp_missing$miss_flag, imp_missing$id,
#                         imp_missing$imputation), ], 20)

## ----means_table, eval=requireNamespace("glmmTMB", quietly = TRUE)------------
# # Mean imputed value per mechanism and treatment arm
# agg <- aggregate(
#   imputed_value ~ miss_flag + trt,
#   data   = imp_result[!is.na(imp_result$miss_flag), ],
#   FUN    = mean
# )
# agg <- agg[order(agg$miss_flag, agg$trt), ]
# print(agg, row.names = FALSE)

## ----pooling, eval=requireNamespace("glmmTMB", quietly = TRUE)----------------
# # Per-imputation mean at visit 4
# v4 <- imp_result[imp_result$visit == 4, ]
# per_imp <- lapply(split(v4, v4$imputation), function(d) {
#   tapply(d$imputed_value, d$trt, mean, na.rm = TRUE)
# })
# 
# # Point estimate: average over imputations
# Q_bar <- Reduce("+", per_imp) / length(per_imp)
# 
# # Within-imputation variance (using single-imputation SE ~ mean/n, illustrative)
# # In practice use regression SE from each imputed dataset
# cat("Pooled mean imputed count at visit 4 by treatment arm:\n")
# print(round(Q_bar, 2))

## ----boot_mi, eval=requireNamespace("glmmTMB", quietly = TRUE)----------------
# imp_boot <- impute_nb(
#   data            = long_data,
#   formula         = mi_formula,
#   outcome_col     = "count",
#   miss_flag_col   = "miss_flag",
#   baseline_col    = "baseline",
#   trt_col         = "trt",
#   reference_trt   = "placebo",
#   subject_col     = "id",
#   strata_cols     = c("trt", "strat1"),
#   n_imp           = 5L,
#   n_boot          = 3L,   # use larger values in practice (e.g., 100–500)
#   seed            = 99L
# )
# 
# # Each replicate × imputation combination is an independent completed dataset
# cat("Total rows:", nrow(imp_boot),
#     "= n_subj * n_visit * n_boot * n_imp =",
#     n_subj * n_visit * 3L * 5L, "\n")
# 
# # Variance of the arm-level mean at visit 4 across all completed datasets
# v4b <- imp_boot[imp_boot$visit == 4, ]
# dataset_means <- tapply(
#   v4b$imputed_value,
#   list(paste(v4b$replicate, v4b$imputation, sep = "_"), v4b$trt),
#   mean
# )
# cat("\nVariance of dataset-level means across boot×imp datasets:\n")
# print(round(apply(dataset_means, 2, var), 4))

## ----step1, eval=requireNamespace("glmmTMB", quietly = TRUE)------------------
# obs_data <- long_data[!is.na(long_data$count), ]
# fits <- fit_nb_glmm(
#   data    = obs_data,
#   formula = count ~ baseline + strat1 + trt + visit + (1 | id)
# )
# cat("Estimated dispersion k:", round(fits[["1"]]$k, 3), "\n")

## ----step2, eval=requireNamespace("glmmTMB", quietly = TRUE)------------------
# imp_mar <- impute_nb_mar(
#   data          = long_data,
#   fits          = fits,
#   outcome_col   = "count",
#   miss_flag_col = "miss_flag",
#   mar_values    = c("MAR", "MNAR"),  # reference-arm MNAR treated as MAR
#   n_imp         = 3L
# )
# cat("MAR imputed rows:", sum(!is.na(imp_mar$imputed_value[
#   imp_mar$imputation == 1 & !is.na(imp_mar$miss_flag)
# ])), "\n")

## ----step3, eval=requireNamespace("glmmTMB", quietly = TRUE)------------------
# imp_mnar <- impute_nb_mnar_ref(
#   data          = long_data,
#   fits          = fits,
#   outcome_col   = "count",
#   miss_flag_col = "miss_flag",
#   mnar_value    = "MNAR",
#   trt_col       = "trt",
#   reference_trt = "placebo",
#   n_imp         = 3L
# )

## ----step4--------------------------------------------------------------------
library(gsDesignNB)

# Composite fill can be applied to any data frame, no glmmTMB needed
df_comp_example <- data.frame(
  count         = c(3L,   NA,     NA,      5L),
  imputed_value = c(3L,    7L,    NA,      5L),
  miss_flag     = c(NA,  "MAR", "Comp",    NA),
  baseline      = c(4L,    4L,    4L,      6L)
)
impute_nb_composite(
  df_comp_example,
  outcome_col     = "count",
  miss_flag_col   = "miss_flag",
  composite_value = "Comp",
  baseline_col    = "baseline"
)

## ----session_info-------------------------------------------------------------
sessionInfo()

