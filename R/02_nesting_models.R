# ==============================================================================
# 02_nesting_models.R
#
# Updated nesting population models for the 2026 SSLL Technical Memorandum.
# Runs, in order:
#   1) loggerhead trend + current abundance (1986-2025)
#   2) leatherback monthly imputation (through season 2025)
#   3) leatherback trend + current abundance, median/low/high branches
#
# The validated Martin et al. (2020) singleUQ model is written below directly so
# the team-facing pipeline does not depend on separate validation scripts.
# ==============================================================================

source("R/00_config.R")

single_uq_model <- "
model{
  A[1] <- 0

  for(j in 2:n.timeseries){
    A[j] ~ dnorm(a_mean,1/(a_sd^2))
  }

  U ~ dnorm(u_mean,1/(u_sd^2))

  tauQ ~ dgamma(q_alpha,q_beta)
  Q <- 1/tauQ

  X0 ~ dnorm(x0_mean,1/(x0_sd^2))
  predX[1] <- X0 + U
  X[1] <- predX[1]

  for(j in 1:n.timeseries){
    tauR[j] ~ dgamma(r_alpha,r_beta)
    R[j] <- 1/tauR[j]
    predY[j,1] <- Z[j,1]*X[1] + A[j]
    Y[j,1] ~ dnorm(predY[j,1],tauR[j])
  }

  for(tt in 2:n.yrs){
    predX[tt] <- X[tt-1] + U
    X[tt] ~ dnorm(predX[tt],tauQ)

    for(j in 1:n.timeseries){
      predY[j,tt] <- Z[j,1]*X[tt] + A[j]
      Y[j,tt] ~ dnorm(predY[j,tt],tauR[j])
    }
  }
}
"

model_dir <- cfg$output_dirs$validation
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(single_uq_model, file.path(model_dir, "singleUQ_validation.txt"))


# ==============================================================================
# LOGGERHEAD UPDATED TREND AND ABUNDANCE
# ==============================================================================

if (!requireNamespace("jagsUI", quietly = TRUE)) {
  stop("Package 'jagsUI' is required.")
}

out_dir <- cfg$output_dirs$update
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nest_file <- "outputs/nesting/cc_annual_nesters.csv"
model_file <- file.path(
  cfg$output_dirs$validation,
  "singleUQ_validation.txt"
)

if (!file.exists(nest_file)) {
  stop("Missing ", nest_file, ". Run 01_prepare_data.R first.")
}

if (!file.exists(model_file)) {
  stop("Missing validated singleUQ model: ", model_file)
}


# ------------------------------------------------------------------------------
# 1. READ PREPARED LOGGERHEAD NESTERS
# ------------------------------------------------------------------------------

cc <- read.csv(
  nest_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

beaches <- c(
  "Inakahama",
  "Maehama",
  "Yotsusehama"
)

required <- c("Year", beaches)

if (!all(required %in% names(cc))) {
  stop(
    "Loggerhead nesting file missing: ",
    paste(setdiff(required, names(cc)), collapse = ", ")
  )
}

cc <- cc[
  cc$Year >= 1986 &
    cc$Year <= 2025,
  required
]

cc <- cc[order(cc$Year), ]

if (!identical(cc$Year, 1986:2025)) {
  stop("Loggerhead model must contain the complete 1986-2025 calendar sequence.")
}

if (any(as.matrix(cc[beaches]) <= 0, na.rm = TRUE)) {
  stop("Non-positive annual-nester value cannot be logged.")
}

all_missing <- cc$Year[
  rowSums(!is.na(cc[beaches])) == 0
]


# ------------------------------------------------------------------------------
# 2. EXACT VALIDATED MARTIN singleUQ INPUT
# ------------------------------------------------------------------------------

data_mat <- t(
  log(
    as.matrix(
      cc[beaches]
    )
  )
)

n.yrs <- ncol(data_mat)
n.timeseries <- nrow(data_mat)

# Exact Martin dimension-preservation construction.
Y <- rbind(
  data_mat,
  NA
)

whichPop <- rep(
  1,
  n.timeseries
)

n.states <- max(whichPop)

Z <- matrix(
  0,
  n.timeseries + 1,
  n.states + 1
)

Z[n.timeseries + 1, ] <- NA
Z[, n.states + 1] <- NA

for (i in seq_along(whichPop)) {
  Z[i, whichPop[i]] <- 1
}


# ------------------------------------------------------------------------------
# 3. EXACT VALIDATED PRIORS
# ------------------------------------------------------------------------------

jags_data <- list(
  
  Y = Y,
  
  n.yrs = n.yrs,
  n.timeseries = n.timeseries,
  
  Z = Z,
  
  a_mean = 0,
  a_sd = 4,
  
  u_mean = 0,
  u_sd = 0.5,
  
  q_alpha = 0.01,
  q_beta = 0.01,
  
  r_alpha = 0.01,
  r_beta = 0.01,
  
  x0_mean = data_mat[1, 1],
  x0_sd = 10
)


# ------------------------------------------------------------------------------
# 4. EXACT VALIDATED MCMC SETTINGS
# ------------------------------------------------------------------------------

n.samples <- 10000
mcmc.chains <- 2
mcmc.thin <- 50
mcmc.burn <- 5000

samples2Save <-
  (mcmc.burn + n.samples) *
  mcmc.thin


cat("\n")
cat("============================================================\n")
cat("LOGGERHEAD 2026 NESTING UPDATE\n")
cat("============================================================\n")
cat("Years: 1986-2025\n")
cat("Clutch frequency: 3 nests/female\n")
cat(
  "All-beach missing years: ",
  if (length(all_missing))
    paste(all_missing, collapse = ", ")
  else
    "none",
  "\n",
  sep = ""
)
cat("Running validated Martin singleUQ model...\n\n")


set.seed(132)

fit <- jagsUI::jags(
  
  data = jags_data,
  
  inits = NULL,
  
  parameters.to.save = c(
    "A",
    "U",
    "Q",
    "R",
    "X0",
    "X"
  ),
  
  model.file = model_file,
  
  n.chains = mcmc.chains,
  
  n.burnin =
    mcmc.burn *
    mcmc.thin,
  
  n.thin =
    mcmc.thin,
  
  n.iter =
    samples2Save,
  
  DIC = TRUE,
  
  parallel = TRUE,
  
  verbose = FALSE
)


# ------------------------------------------------------------------------------
# 5. SAVE FIT
# ------------------------------------------------------------------------------

fit_file <- file.path(
  out_dir,
  "update_trend_loggerhead_1986_2025.rds"
)

saveRDS(
  fit,
  fit_file
)


# ------------------------------------------------------------------------------
# 6. DIAGNOSTICS
# ------------------------------------------------------------------------------

s <- as.data.frame(
  fit$summary
)

finite_rhat <- s$Rhat[
  is.finite(s$Rhat)
]

max_rhat <- if (length(finite_rhat)) {
  max(finite_rhat)
} else {
  NA_real_
}

neff_name <- intersect(
  c("n.eff", "n_eff"),
  names(s)
)

if (length(neff_name)) {
  
  # Exclude deterministic/fixed nodes from effective sample-size diagnostics.
  # In the Martin singleUQ model A[1] is fixed at 0 and therefore has
  # sd = 0, Rhat = NA, and n.eff = 1; it is not an estimated parameter.
  stochastic_neff <- s[[neff_name[1]]][
    is.finite(s[[neff_name[1]]]) &
      is.finite(s$sd) &
      s$sd > 0
  ]
  
  min_neff <- if (length(stochastic_neff)) {
    min(stochastic_neff)
  } else {
    NA_real_
  }
  
} else {
  
  min_neff <- NA_real_
}


# ------------------------------------------------------------------------------
# 7. EXACT MARTIN POSTERIOR THINNING
# ------------------------------------------------------------------------------

U <- as.numeric(
  fit$sims.list$U
)

Q <- as.numeric(
  fit$sims.list$Q
)

X <- as.matrix(
  fit$sims.list$X
)

A <- as.matrix(
  fit$sims.list$A
)

nsim <- 10000L

X_len <- nrow(X)

if (X_len < nsim) {
  stop("Fewer than 10,000 posterior X draws available.")
}

thin_id <- floor(
  seq(
    from = 1,
    to = X_len,
    length.out = nsim
  )
)

thin_id <- unique(
  pmax(
    1L,
    pmin(
      X_len,
      thin_id
    )
  )
)

if (length(thin_id) != nsim) {
  
  step <- X_len / nsim
  
  thin_id <- as.integer(
    seq(
      from = 1,
      by = step,
      length.out = nsim
    )
  )
}


# ------------------------------------------------------------------------------
# 8. REGIONAL ANNUAL NESTERS
# ------------------------------------------------------------------------------

regional_nesters <- function(year_index) {
  
  x <- X[
    thin_id,
    year_index
  ]
  
  total <- exp(x)
  
  if (n.timeseries > 1L) {
    
    for (j in 2:n.timeseries) {
      
      total <-
        total +
        exp(
          x +
            A[thin_id, j]
        )
    }
  }
  
  total
}


fy <- ncol(X)

N_fym0 <- regional_nesters(fy)
N_fym1 <- regional_nesters(fy - 1L)
N_fym2 <- regional_nesters(fy - 2L)
N_fym3 <- regional_nesters(fy - 3L)


# ------------------------------------------------------------------------------
# 9. CURRENT TOTAL NESTER ABUNDANCE
#
# Exact Martin curr.abund.fn:
# final four annual nester estimates * RI / 4
# Loggerhead RI = 3.3
# ------------------------------------------------------------------------------

RI <- cfg$biology$cc$RI

current_abundance <-
  rowSums(
    cbind(
      N_fym0,
      N_fym1,
      N_fym2,
      N_fym3
    )
  ) *
  RI / 4


q3 <- function(x) {
  
  unname(
    quantile(
      x,
      probs = c(
        0.025,
        0.5,
        0.975
      ),
      na.rm = TRUE
    )
  )
}


r_q <- q3(U[thin_id])

annual_pct_q <-
  (
    exp(r_q) -
      1
  ) *
  100

abundance_q <-
  q3(current_abundance)

final_year_q <-
  q3(N_fym0)


# ------------------------------------------------------------------------------
# 10. SAVE POSTERIOR INPUT FOR LATER HISTORICAL-TAKE/PVA STEPS
# ------------------------------------------------------------------------------

posts <- data.frame(
  
  U = U[thin_id],
  
  Q = Q[thin_id],
  
  N_fym0 = N_fym0,
  N_fym1 = N_fym1,
  N_fym2 = N_fym2,
  N_fym3 = N_fym3
)

write.csv(
  posts,
  file.path(
    out_dir,
    "update_posteriors_loggerhead_1986_2025.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 11. SUMMARY TABLE
# ------------------------------------------------------------------------------

summary_out <- data.frame(
  
  species = "Loggerhead",
  
  start_year = 1986,
  end_year = 2025,
  
  r_median = r_q[2],
  r_L95 = r_q[1],
  r_U95 = r_q[3],
  
  annual_percent_median =
    annual_pct_q[2],
  
  annual_percent_L95 =
    annual_pct_q[1],
  
  annual_percent_U95 =
    annual_pct_q[3],
  
  final_year_nesters_median =
    final_year_q[2],
  
  final_year_nesters_L95 =
    final_year_q[1],
  
  final_year_nesters_U95 =
    final_year_q[3],
  
  current_abundance_median =
    abundance_q[2],
  
  current_abundance_L95 =
    abundance_q[1],
  
  current_abundance_U95 =
    abundance_q[3],
  
  max_Rhat =
    max_rhat,
  
  min_n_eff =
    min_neff
)

write.csv(
  summary_out,
  file.path(
    out_dir,
    "update_loggerhead_trend_abundance_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 12. REPORT
# ------------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("LOGGERHEAD UPDATE RESULTS\n")
cat("============================================================\n")

cat(
  sprintf(
    "r: %.4f (95%% CI %.4f to %.4f)\n",
    r_q[2],
    r_q[1],
    r_q[3]
  )
)

cat(
  sprintf(
    "annual %%: %.2f%% (95%% CI %.2f%% to %.2f%%)\n",
    annual_pct_q[2],
    annual_pct_q[1],
    annual_pct_q[3]
  )
)

cat(
  sprintf(
    "2025 annual nesters: %.0f (95%% CI %.0f to %.0f)\n",
    final_year_q[2],
    final_year_q[1],
    final_year_q[3]
  )
)

cat(
  sprintf(
    "current total abundance: %.0f (95%% CI %.0f to %.0f)\n",
    abundance_q[2],
    abundance_q[1],
    abundance_q[3]
  )
)

cat(
  sprintf(
    "max Rhat: %.4f\n",
    max_rhat
  )
)

cat(
  sprintf(
    "min n_eff: %.0f\n",
    min_neff
  )
)

cat("\nCurrent abundance uses fitted years 2022-2025 and RI = 3.3.\n")
cat("This is NESTING-ONLY abundance; historical SSLL ANE has not yet been added back.\n")
cat("============================================================\n")


# ==============================================================================
# LEATHERBACK UPDATED MONTHLY IMPUTATION
# ==============================================================================

if (!requireNamespace("jagsUI", quietly = TRUE)) {
  stop("Package 'jagsUI' is required.")
}

out_dir <- cfg$output_dirs$update
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- "outputs/nesting/dc_monthly_season_grid.csv"

if (!file.exists(input_file)) {
  stop("Missing ", input_file, ". Run 01_prepare_data.R first.")
}

dc <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c(
  "Season",
  "JM_Nests",
  "W_Nests"
)

if (!all(required %in% names(dc))) {
  stop(
    "Missing columns: ",
    paste(setdiff(required, names(dc)), collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 1. UPDATED ANALYSIS WINDOW
# ------------------------------------------------------------------------------

dc <- dc[
  dc$Season >= 2001L &
    dc$Season <= 2025L,
]

dc <- dc[
  order(dc$Season, dc$Seq_month),
]

seasons <- 2001:2025
n_years <- length(seasons)

season_counts <- table(dc$Season)

if (
  length(season_counts) != n_years ||
  any(season_counts != 12)
) {
  stop("Expected exactly 12 monthly rows for every season 2001-2025.")
}

# Martin W series begins in season 2006.
dc$W_Nests[
  dc$Season < 2006L
] <- NA_real_


# ------------------------------------------------------------------------------
# 2. CHECK OBSERVED VALUES
#
# Exact Martin model logs observed nest counts.
# Do not silently alter zeros/non-positive values.
# ------------------------------------------------------------------------------

bad <- dc[
  (!is.na(dc$JM_Nests) & dc$JM_Nests <= 0) |
    (!is.na(dc$W_Nests)  & dc$W_Nests <= 0),
]

if (nrow(bad)) {
  
  print(bad)
  
  stop(
    "Observed zero/non-positive leatherback nest counts occur in the ",
    "updated analysis window. Martin's model logs observed counts directly."
  )
}


# ------------------------------------------------------------------------------
# 3. EXACT MARTIN DATA ORIENTATION
# ------------------------------------------------------------------------------

jm_mat <- matrix(
  log(dc$JM_Nests),
  nrow = n_years,
  ncol = 12,
  byrow = TRUE
)

w_mat <- matrix(
  log(dc$W_Nests),
  nrow = n_years,
  ncol = 12,
  byrow = TRUE
)

y_dc <- cbind(
  as.vector(t(jm_mat)),
  as.vector(t(w_mat))
)

dc_jags_data <- list(
  
  y = y_dc,
  
  m = rep(
    1:12,
    times = n_years
  ),
  
  n.steps = nrow(y_dc),
  
  n.months = 12,
  
  pi = pi,
  
  period = c(
    12,
    6
  ),
  
  n.timeseries = 2,
  
  n.years = n_years
)


# ------------------------------------------------------------------------------
# 4. EXACT VALIDATED MARTIN IMPUTATION MODEL
# ------------------------------------------------------------------------------

imputation_model <- "
model{

    for(j in 1:n.timeseries) {

       predX0[j] ~ dnorm(5, 0.1)

       predX[1,j] <- c[j, m[1]] + predX0[j]

       X[1,j] ~ dnorm(
          predX[1,j],
          tau.X[j]
       )

       y[1,j] ~ dnorm(
          X[1,j],
          tau.y[j]
       )

       for(t in 2:n.steps) {

           predX[t,j] <-
               c[j,m[t]] +
               X[t-1,j]

           X[t,j] ~ dnorm(
               predX[t,j],
               tau.X[j]
           )

           y[t,j] ~ dnorm(
               X[t,j],
               tau.y[j]
           )
       }

       for(yy in 1:n.years) {

           for(mm in 1:12) {

               tmp2[yy,mm,j] <-
                   exp(
                       X[
                           (yy*12 - mm + 1),
                           j
                       ]
                   )
           }

           N[yy,j] <-
               log(
                   sum(
                       tmp2[yy,,j]
                   )
               )
       }
    }

    for(j in 1:n.timeseries) {

        for(k in 1:n.months) {

            c.const[j,k] <-
                2 * pi * k /
                period[j]

            c[j,k] <-
                beta.cos[j] *
                cos(c.const[j,k]) +

                beta.sin[j] *
                sin(c.const[j,k])
        }

        sigma.y[j] ~ dgamma(2, 0.5)

        tau.y[j] <-
            1 /
            (
                sigma.y[j] *
                sigma.y[j]
            )

        beta.cos[j] ~ dnorm(0, 1)

        beta.sin[j] ~ dnorm(0, 1)

        sigma.X[j] ~ dgamma(2, 0.5)

        tau.X[j] <-
            1 /
            (
                sigma.X[j] *
                sigma.X[j]
            )
    }
}
"

model_file <- file.path(
  out_dir,
  "model_norm_norm_Four_imputation_update.txt"
)

writeLines(
  imputation_model,
  model_file
)


# ------------------------------------------------------------------------------
# 5. RUN
#
# Martin baseline:
#   100000 iterations
#   50000 burn-in
#   thin 5
#
# Updated run:
#   SAME MODEL
#   longer burn-in and sampling path only
# ------------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("LEATHERBACK 2026 IMPUTATION UPDATE\n")
cat("============================================================\n")
cat("Seasons: 2001-2025\n")
cat("Monthly rows: ", nrow(y_dc), "\n", sep = "")
cat(
  "JM observed months: ",
  sum(!is.na(y_dc[, 1])),
  " / ",
  nrow(y_dc),
  "\n",
  sep = ""
)
cat(
  "W observed months: ",
  sum(!is.na(y_dc[, 2])),
  " / ",
  nrow(y_dc),
  "\n",
  sep = ""
)
cat("\n")
cat("Exact validated Martin imputation model\n")
cat("MCMC tuning only:\n")
cat("  5 chains\n")
cat("  1000000 iterations\n")
cat("  500000 burn-in\n")
cat("  thin 20\n\n")


# Explicit chain-specific RNG states make the parallel JAGS run reproducible.
# These affect only random-number initialization, not the model or priors.
dc_inits <- lapply(
  seq_len(5),
  function(i) {
    list(
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = 132L + i - 1L
    )
  }
)

fit <- jagsUI::jags(
  
  data = dc_jags_data,
  
  inits = dc_inits,
  
  parameters.to.save = c(
    "c",
    "beta.cos",
    "beta.sin",
    "sigma.X",
    "sigma.y",
    "N",
    "y",
    "X"
  ),
  
  model.file = model_file,
  
  n.chains = 5,
  
  n.iter = 1000000,
  
  n.burnin = 500000,
  
  n.thin = 20,
  
  DIC = TRUE,
  
  parallel = TRUE,
  
  verbose = FALSE
)


saveRDS(
  fit,
  file.path(
    out_dir,
    "update_leatherback_imputation_fit.rds"
  )
)


# ------------------------------------------------------------------------------
# 6. CONVERGENCE
# ------------------------------------------------------------------------------

s <- as.data.frame(
  fit$summary
)

s$parameter <- rownames(s)

finite <- s[
  is.finite(s$Rhat),
]

finite <- finite[
  order(
    finite$Rhat,
    decreasing = TRUE
  ),
]

max_rhat <- max(
  finite$Rhat
)

cat("\n")
cat("============================================================\n")
cat("LEATHERBACK IMPUTATION CONVERGENCE\n")
cat("============================================================\n")

cat(
  sprintf(
    "Maximum finite Rhat: %.4f\n\n",
    max_rhat
  )
)

cat("10 highest Rhat values:\n")

print(
  finite[
    seq_len(
      min(
        10,
        nrow(finite)
      )
    ),
    c(
      "parameter",
      "mean",
      "sd",
      "Rhat"
    )
  ],
  row.names = FALSE
)


# Do not propagate an unconverged imputation into the trend/PVA pipeline.
if (max_rhat > 1.05) {
  stop(
    sprintf(
      paste0(
        "Leatherback imputation did not meet the convergence criterion ",
        "(maximum Rhat = %.4f > 1.05). ",
        "The fit has been saved for inspection, but annual imputed nests ",
        "and downstream trend estimates were not updated."
      ),
      max_rhat
    )
  )
}

# ------------------------------------------------------------------------------
# 7. ANNUAL IMPUTED NESTS
# ------------------------------------------------------------------------------

log_low <- as.matrix(
  fit$q2.5$N
)

log_med <- as.matrix(
  fit$q50$N
)

log_high <- as.matrix(
  fit$q97.5$N
)

if (
  !identical(
    dim(log_med),
    c(n_years, 2L)
  )
) {
  stop("Unexpected dimensions for annual N output.")
}


annual <- data.frame(
  
  Season = seasons,
  
  JM_low =
    exp(log_low[, 1]),
  
  JM_median =
    exp(log_med[, 1]),
  
  JM_high =
    exp(log_high[, 1]),
  
  W_low =
    exp(log_low[, 2]),
  
  W_median =
    exp(log_med[, 2]),
  
  W_high =
    exp(log_high[, 2])
)


write.csv(
  annual,
  file.path(
    out_dir,
    "update_leatherback_annual_imputed_nests_2001_2025.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8. REPORT
# ------------------------------------------------------------------------------

cat("\nFinal four seasons — median annual nests:\n")

print(
  annual[
    annual$Season >= 2021,
    c(
      "Season",
      "JM_median",
      "W_median"
    )
  ],
  row.names = FALSE
)

cat("\n")

if (max_rhat <= 1.05) {
  
  cat("CONVERGENCE STATUS: PASS\n")
  cat("Proceed to leatherback trend + abundance.\n")
  
} else if (max_rhat <= 1.10) {
  
  cat("CONVERGENCE STATUS: REVIEW\n")
  cat("Inspect the parameters above before proceeding.\n")
  
} else {
  
  cat("CONVERGENCE STATUS: FAIL\n")
  cat("Do not run the updated leatherback trend yet.\n")
}

cat("============================================================\n")


# ==============================================================================
# LEATHERBACK UPDATED TREND AND ABUNDANCE
# ==============================================================================

if (!requireNamespace("jagsUI", quietly = TRUE)) {
  stop("Package 'jagsUI' is required.")
}

out_dir <- cfg$output_dirs$update
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(
  out_dir,
  "update_leatherback_annual_imputed_nests_2001_2025.csv"
)

model_file <- file.path(
  cfg$output_dirs$validation,
  "singleUQ_validation.txt"
)

if (!file.exists(input_file)) {
  stop("Missing leatherback imputation output created earlier in 02_nesting_models.R: ", input_file)
}

if (!file.exists(model_file)) {
  stop("Missing validated singleUQ model: ", model_file)
}


# ------------------------------------------------------------------------------
# 1. READ IMPUTED ANNUAL NESTS
# ------------------------------------------------------------------------------

d <- read.csv(
  input_file,
  stringsAsFactors = FALSE
)

if (!identical(d$Season, 2001:2025)) {
  stop("Expected leatherback seasons 2001-2025.")
}

CF <- cfg$biology$dc$CF
RI <- cfg$biology$dc$RI


make_branch <- function(branch) {
  
  jm_col <- paste0("JM_", branch)
  w_col  <- paste0("W_", branch)
  
  out <- data.frame(
    Season = d$Season,
    JM = d[[jm_col]] / CF,
    W  = d[[w_col]] / CF
  )
  
  # Exact Martin site window.
  out$W[out$Season < 2006L] <- NA_real_
  
  out
}


dc_median <- make_branch("median")
dc_low    <- make_branch("low")
dc_high   <- make_branch("high")


# ------------------------------------------------------------------------------
# 2. HELPERS
# ------------------------------------------------------------------------------

q3 <- function(x) {
  unname(
    quantile(
      x,
      c(0.025, 0.5, 0.975),
      na.rm = TRUE
    )
  )
}


max_rhat <- function(fit) {
  
  s <- as.data.frame(fit$summary)
  
  z <- s$Rhat[
    is.finite(s$Rhat)
  ]
  
  if (length(z)) max(z) else NA_real_
}


# ------------------------------------------------------------------------------
# 3. EXACT VALIDATED singleUQ FIT
# ------------------------------------------------------------------------------

run_trend <- function(dat, label) {
  
  years <- dat$Season
  
  mat <- as.matrix(
    dat[, c("JM", "W")]
  )
  
  if (any(mat <= 0, na.rm = TRUE)) {
    stop(label, ": non-positive annual-nester value.")
  }
  
  data_mat <- t(log(mat))
  
  n.yrs <- ncol(data_mat)
  n.timeseries <- nrow(data_mat)
  
  # Exact Martin dimension-preservation construction.
  Y <- rbind(
    data_mat,
    NA
  )
  
  whichPop <- rep(
    1,
    n.timeseries
  )
  
  n.states <- max(whichPop)
  
  Z <- matrix(
    0,
    n.timeseries + 1,
    n.states + 1
  )
  
  Z[n.timeseries + 1, ] <- NA
  Z[, n.states + 1] <- NA
  
  for (i in seq_along(whichPop)) {
    Z[i, whichPop[i]] <- 1
  }
  
  
  jags_data <- list(
    
    Y = Y,
    
    n.yrs = n.yrs,
    n.timeseries = n.timeseries,
    
    Z = Z,
    
    a_mean = 0,
    a_sd = 4,
    
    u_mean = 0,
    u_sd = 0.5,
    
    q_alpha = 0.01,
    q_beta = 0.01,
    
    r_alpha = 0.01,
    r_beta = 0.01,
    
    x0_mean = data_mat[1, 1],
    x0_sd = 10
  )
  
  
  # Exact validated Martin trend settings.
  n.samples <- 10000
  mcmc.chains <- 2
  mcmc.thin <- 50
  mcmc.burn <- 5000
  
  samples2Save <-
    (mcmc.burn + n.samples) *
    mcmc.thin
  
  
  cat("\nRunning leatherback ", label, " trend...\n", sep = "")
  
  set.seed(132)
  
  fit <- jagsUI::jags(
    
    data = jags_data,
    
    inits = NULL,
    
    parameters.to.save = c(
      "A",
      "U",
      "Q",
      "R",
      "X0",
      "X"
    ),
    
    model.file = model_file,
    
    n.chains = mcmc.chains,
    
    n.burnin =
      mcmc.burn *
      mcmc.thin,
    
    n.thin =
      mcmc.thin,
    
    n.iter =
      samples2Save,
    
    DIC = TRUE,
    
    parallel = TRUE,
    
    verbose = FALSE
  )
  
  
  saveRDS(
    fit,
    file.path(
      out_dir,
      paste0(
        "update_trend_leatherback_",
        label,
        "_2001_2025.rds"
      )
    )
  )
  
  
  # --------------------------------------------------------------------------
  # Posterior thinning — exact validated approach
  # --------------------------------------------------------------------------
  
  X <- as.matrix(
    fit$sims.list$X
  )
  
  A <- as.matrix(
    fit$sims.list$A
  )
  
  U <- as.numeric(
    fit$sims.list$U
  )
  
  Q <- as.numeric(
    fit$sims.list$Q
  )
  
  nsim <- 10000L
  X_len <- nrow(X)
  
  if (X_len < nsim) {
    stop(label, ": fewer than 10,000 posterior draws.")
  }
  
  thin_id <- floor(
    seq(
      1,
      X_len,
      length.out = nsim
    )
  )
  
  thin_id <- unique(
    pmax(
      1L,
      pmin(X_len, thin_id)
    )
  )
  
  if (length(thin_id) != nsim) {
    
    step <- X_len / nsim
    
    thin_id <- as.integer(
      seq(
        1,
        by = step,
        length.out = nsim
      )
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Regional annual nesters
  # --------------------------------------------------------------------------
  
  regional_nesters <- function(year_index) {
    
    x <- X[
      thin_id,
      year_index
    ]
    
    total <- exp(x)
    
    for (j in 2:n.timeseries) {
      
      total <-
        total +
        exp(
          x +
            A[thin_id, j]
        )
    }
    
    total
  }
  
  
  fy <- ncol(X)
  
  N_fym0 <- regional_nesters(fy)
  N_fym1 <- regional_nesters(fy - 1L)
  N_fym2 <- regional_nesters(fy - 2L)
  N_fym3 <- regional_nesters(fy - 3L)
  
  
  posts <- data.frame(
    
    U = U[thin_id],
    
    Q = Q[thin_id],
    
    N_fym0 = N_fym0,
    N_fym1 = N_fym1,
    N_fym2 = N_fym2,
    N_fym3 = N_fym3
  )
  
  
  write.csv(
    posts,
    file.path(
      out_dir,
      paste0(
        "update_posteriors_leatherback_",
        label,
        "_2001_2025.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # --------------------------------------------------------------------------
  # Current abundance — exact Martin calculation
  # --------------------------------------------------------------------------
  
  current_abundance <-
    rowSums(
      posts[
        c(
          "N_fym0",
          "N_fym1",
          "N_fym2",
          "N_fym3"
        )
      ]
    ) *
    RI / 4
  
  
  r_q <- q3(
    U[thin_id]
  )
  
  pct_q <-
    (
      exp(r_q) -
        1
    ) *
    100
  
  final_q <- q3(
    N_fym0
  )
  
  abundance_q <- q3(
    current_abundance
  )
  
  
  data.frame(
    
    branch = toupper(label),
    
    start_season = 2001,
    end_season = 2025,
    
    r_median = r_q[2],
    r_L95 = r_q[1],
    r_U95 = r_q[3],
    
    annual_percent_median = pct_q[2],
    annual_percent_L95 = pct_q[1],
    annual_percent_U95 = pct_q[3],
    
    final_annual_nesters_median = final_q[2],
    final_annual_nesters_L95 = final_q[1],
    final_annual_nesters_U95 = final_q[3],
    
    current_abundance_median = abundance_q[2],
    current_abundance_L95 = abundance_q[1],
    current_abundance_U95 = abundance_q[3],
    
    max_Rhat = max_rhat(fit)
  )
}


# ------------------------------------------------------------------------------
# 4. RUN ALL THREE MARTIN IMPUTATION BRANCHES
# ------------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("LEATHERBACK 2026 TREND + ABUNDANCE UPDATE\n")
cat("============================================================\n")
cat("Seasons: 2001-2025\n")
cat("Clutch frequency: ", CF, "\n", sep = "")
cat("Remigration interval: ", RI, "\n", sep = "")
cat("Branches: MEDIAN, LOW, HIGH\n")
cat("Using validated Martin singleUQ model unchanged.\n")
cat("============================================================\n")


med <- run_trend(
  dc_median,
  "MEDIAN"
)

low <- run_trend(
  dc_low,
  "LOW"
)

high <- run_trend(
  dc_high,
  "HIGH"
)


results <- rbind(
  med,
  low,
  high
)


write.csv(
  results,
  file.path(
    out_dir,
    "update_leatherback_trend_abundance_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 5. REPORT
# ------------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("LEATHERBACK UPDATE RESULTS\n")
cat("============================================================\n")

for (i in seq_len(nrow(results))) {
  
  z <- results[i, ]
  
  cat("\n", z$branch, "\n", sep = "")
  
  cat(
    sprintf(
      "  r: %.4f (95%% CI %.4f to %.4f)\n",
      z$r_median,
      z$r_L95,
      z$r_U95
    )
  )
  
  cat(
    sprintf(
      "  annual %%: %.2f%% (95%% CI %.2f%% to %.2f%%)\n",
      z$annual_percent_median,
      z$annual_percent_L95,
      z$annual_percent_U95
    )
  )
  
  cat(
    sprintf(
      "  2025 annual nesters: %.0f (95%% CI %.0f to %.0f)\n",
      z$final_annual_nesters_median,
      z$final_annual_nesters_L95,
      z$final_annual_nesters_U95
    )
  )
  
  cat(
    sprintf(
      "  current total abundance: %.0f (95%% CI %.0f to %.0f)\n",
      z$current_abundance_median,
      z$current_abundance_L95,
      z$current_abundance_U95
    )
  )
  
  cat(
    sprintf(
      "  max Rhat: %.4f\n",
      z$max_Rhat
    )
  )
}

cat("\n")
cat("Current abundance uses fitted seasons 2021-2025 and RI = 3.06.\n")
cat("These are NESTING-ONLY estimates; updated historical SSLL ANE has not yet been added back.\n")
cat("============================================================\n")
