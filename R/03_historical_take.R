# ==============================================================================
# 03_historical_take.R
#
# Historical SSLL take for the 2026 update.
#   1) Calculate historical Adult Nester Equivalents (ANE).
#   2) Add historical SSLL ANE back to nesting series.
#   3) Rerun the frozen Martin singleUQ model to obtain no-historical-take
#      population posteriors used by the PVA.
# ==============================================================================

source("R/00_config.R")

in_dir  <- "outputs/input_audit"
out_dir <- cfg$output_dirs$update

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

fishery_start <- 2005L
historical_RI <- 3L
age_grid <- seq(0, 100, by = 0.01)

trend_end <- c(
  cc = 2025L,
  dc = 2025L
)

species_name <- c(
  cc = "Loggerhead",
  dc = "Leatherback"
)


# ------------------------------------------------------------------------------
# 1. READ CURRENT COMBINED PIRO INPUT
# ------------------------------------------------------------------------------

interaction_file <- file.path(
  in_dir,
  "ssll_interactions_prepared.csv"
)

coverage_file <- file.path(
  in_dir,
  "ssll_annual_counts.csv"
)

if (!file.exists(interaction_file)) {
  stop("Missing: ", interaction_file)
}

if (!file.exists(coverage_file)) {
  stop("Missing: ", coverage_file)
}

ssll <- read.csv(
  interaction_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

coverage <- read.csv(
  coverage_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

needed_ssll <- c(
  "year",
  "species",
  "scl_cm",
  "m_mean",
  "source_row"
)

if (!all(needed_ssll %in% names(ssll))) {
  stop(
    "Combined SSLL file missing: ",
    paste(
      setdiff(needed_ssll, names(ssll)),
      collapse = ", "
    )
  )
}

needed_coverage <- c(
  "year",
  "species",
  "n",
  "status"
)

if (!all(needed_coverage %in% names(coverage))) {
  stop(
    "Coverage file missing: ",
    paste(
      setdiff(needed_coverage, names(coverage)),
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# 2. STORAGE
# ------------------------------------------------------------------------------

annual_results <- list()
summary_results <- list()


# ------------------------------------------------------------------------------
# 3. RUN MARTIN HISTORICAL ANE BY SPECIES
# ------------------------------------------------------------------------------

for (spp in c("cc", "dc")) {
  
  spp_name <- species_name[[spp]]
  b <- cfg$biology[[spp]]
  end_year <- trend_end[[spp]]
  
  
  # --------------------------------------------------------------------------
  # Coverage status
  # --------------------------------------------------------------------------
  
  cov_spp <- coverage[
    coverage$species == spp_name &
      coverage$year >= fishery_start &
      coverage$year <= end_year,
    ,
    drop = FALSE
  ]
  
  pending_years <- cov_spp$year[
    toupper(cov_spp$status) == "PENDING"
  ]
  
  confirmed_zero_years <- cov_spp$year[
    tolower(cov_spp$status) == "complete" &
      cov_spp$n == 0
  ]
  
  
  # --------------------------------------------------------------------------
  # Select current PIRO records
  # --------------------------------------------------------------------------
  
  d <- ssll[
    ssll$species == spp_name &
      ssll$year >= fishery_start &
      ssll$year <= end_year,
    ,
    drop = FALSE
  ]
  
  n_before_pending_exclusion <- nrow(d)
  
  # Pending years are incomplete.
  # Exclude their currently available partial records entirely.
  if (length(pending_years) > 0L) {
    
    d <- d[
      !d$year %in% pending_years,
      ,
      drop = FALSE
    ]
  }
  
  n_pending_records_excluded <-
    n_before_pending_exclusion - nrow(d)
  
  if (nrow(d) == 0L) {
    stop(
      spp_name,
      ": no complete-year interaction records remain."
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Standardize names for Martin ANE calculation
  # --------------------------------------------------------------------------
  
  d <- data.frame(
    Year = d$year,
    SCL_cm = d$scl_cm,
    M_mean = d$m_mean,
    source_row = d$source_row,
    stringsAsFactors = FALSE
  )
  
  
  # --------------------------------------------------------------------------
  # Missing SCL and mortality — Martin method
  # --------------------------------------------------------------------------
  
  d$length_imputed <- is.na(d$SCL_cm)
  d$mortality_imputed <- is.na(d$M_mean)
  
  median_length <- median(
    d$SCL_cm,
    na.rm = TRUE
  )
  
  mean_mortality <- mean(
    d$M_mean,
    na.rm = TRUE
  )
  
  d$Len_used_cm <- d$SCL_cm
  d$M_used <- d$M_mean
  
  d$Len_used_cm[d$length_imputed] <-
    median_length
  
  d$M_used[d$mortality_imputed] <-
    mean_mortality
  
  if (anyNA(d$Len_used_cm)) {
    stop(spp_name, ": missing SCL remains.")
  }
  
  if (anyNA(d$M_used)) {
    stop(spp_name, ": missing mortality remains.")
  }
  
  if (any(d$M_used < 0 | d$M_used > 1)) {
    stop(spp_name, ": mortality outside 0-1.")
  }
  
  
  # --------------------------------------------------------------------------
  # Back-calculate age from SCL — Martin method
  # --------------------------------------------------------------------------
  
  if (spp == "cc") {
    
    predicted_length <-
      b$Linf -
      (b$Linf - b$Lknot) *
      exp(-b$K * age_grid)
    
  } else {
    
    predicted_length <-
      b$Linf *
      (
        1 -
          exp(
            -b$K *
              (age_grid - b$tknot)
          )
      )
  }
  
  d$Age <- vapply(
    d$Len_used_cm,
    function(x) {
      
      age_grid[
        which.min(
          abs(x - predicted_length)
        )
      ]
    },
    numeric(1)
  )
  
  d$Stage <- ifelse(
    d$Age > b$Amat,
    "A",
    "J"
  )
  
  d$years_to_maturity <-
    b$Amat - d$Age
  
  d$years_to_first_nesting <-
    pmax(
      0,
      round(d$years_to_maturity)
    )
  
  d$first_nesting_year <-
    d$Year +
    d$years_to_first_nesting
  
  
  # --------------------------------------------------------------------------
  # Historical ANE — Martin method
  # --------------------------------------------------------------------------
  
  events <- lapply(
    seq_len(nrow(d)),
    function(i) {
      
      first_year <-
        d$first_nesting_year[i]
      
      # Interaction cannot affect fitted nesting series
      # if first nesting occurs after its endpoint.
      if (first_year > end_year) {
        return(NULL)
      }
      
      nesting_years <- seq(
        from = first_year,
        to = end_year,
        by = historical_RI
      )
      
      elapsed <-
        nesting_years -
        d$Year[i]
      
      first_ane <-
        (b$Pj ^
           d$years_to_first_nesting[i]) *
        b$PF *
        d$M_used[i]
      
      ane <- rep(
        first_ane,
        length(nesting_years)
      )
      
      later <-
        seq_along(nesting_years) > 1L
      
      ane[later] <-
        (b$Pa ^ elapsed[later]) *
        first_ane
      
      data.frame(
        species = spp_name,
        source_row = d$source_row[i],
        interaction_year = d$Year[i],
        nesting_year = nesting_years,
        occasion = seq_along(nesting_years),
        elapsed_years = elapsed,
        ANE = ane,
        stringsAsFactors = FALSE
      )
    }
  )
  
  events <- Filter(
    Negate(is.null),
    events
  )
  
  if (length(events) > 0L) {
    
    event_table <- do.call(
      rbind,
      events
    )
    
    totals <- aggregate(
      ANE ~ nesting_year,
      data = event_table,
      FUN = sum
    )
    
  } else {
    
    event_table <- data.frame(
      species = character(),
      source_row = integer(),
      interaction_year = integer(),
      nesting_year = integer(),
      occasion = integer(),
      elapsed_years = numeric(),
      ANE = numeric()
    )
    
    totals <- data.frame(
      nesting_year = integer(),
      ANE = numeric()
    )
  }
  
  
  # --------------------------------------------------------------------------
  # Complete annual ANE series
  # --------------------------------------------------------------------------
  
  annual <- merge(
    data.frame(
      nesting_year =
        seq.int(
          fishery_start,
          end_year
        )
    ),
    totals,
    by = "nesting_year",
    all.x = TRUE,
    sort = TRUE
  )
  
  annual$ANE[
    is.na(annual$ANE)
  ] <- 0
  
  annual$species <- spp_name
  
  annual <- annual[
    c(
      "species",
      "nesting_year",
      "ANE"
    )
  ]
  
  
  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------
  
  write.csv(
    d,
    file.path(
      out_dir,
      paste0(
        spp,
        "_historical_ane_interactions_update.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    event_table,
    file.path(
      out_dir,
      paste0(
        spp,
        "_historical_ane_events_update.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    annual,
    file.path(
      out_dir,
      paste0(
        spp,
        "_historical_ane_annual_update.csv"
      )
    ),
    row.names = FALSE
  )
  
  annual_results[[spp]] <- annual
  
  
  # --------------------------------------------------------------------------
  # Summary
  # --------------------------------------------------------------------------
  
  summary_results[[spp]] <- data.frame(
    
    species = spp_name,
    
    interactions_used = nrow(d),
    
    pending_records_excluded =
      n_pending_records_excluded,
    
    interaction_start =
      min(d$Year),
    
    interaction_end =
      max(d$Year),
    
    trend_end =
      end_year,
    
    lengths_imputed =
      sum(d$length_imputed),
    
    median_SCL_cm =
      median_length,
    
    mortalities_imputed =
      sum(d$mortality_imputed),
    
    mean_mortality =
      mean_mortality,
    
    juveniles =
      sum(d$Stage == "J"),
    
    adults =
      sum(d$Stage == "A"),
    
    ANE_through_trend_end =
      sum(annual$ANE),
    
    pending_years =
      if (length(pending_years)) {
        paste(
          pending_years,
          collapse = ","
        )
      } else {
        "none"
      },
    
    confirmed_zero_years =
      if (length(confirmed_zero_years)) {
        paste(
          confirmed_zero_years,
          collapse = ","
        )
      } else {
        "none"
      },
    
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# 4. COMBINED OUTPUTS
# ------------------------------------------------------------------------------

summary_table <- do.call(
  rbind,
  summary_results
)

annual_table <- do.call(
  rbind,
  annual_results
)

write.csv(
  summary_table,
  file.path(
    out_dir,
    "historical_ane_update_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  annual_table,
  file.path(
    out_dir,
    "historical_ane_update_both_species.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 5. REPORT
# ------------------------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("UPDATED HISTORICAL SSLL ANE — MARTIN METHOD\n")
cat("============================================================\n")
cat("Input: ssll_interactions_prepared.csv\n")
cat("PIRO interaction start:", fishery_start, "\n")
cat("Historical nesting interval:", historical_RI, "years\n")
cat("Loggerhead trend endpoint:", trend_end["cc"], "\n")
cat("Leatherback trend endpoint:", trend_end["dc"], "\n\n")

for (spp in c("cc", "dc")) {
  
  x <- summary_results[[spp]]
  
  cat(toupper(x$species), "\n")
  
  cat(
    "  Interactions used:",
    x$interactions_used,
    "\n"
  )
  
  cat(
    "  Pending records excluded:",
    x$pending_records_excluded,
    "\n"
  )
  
  cat(
    "  Interaction years used:",
    x$interaction_start,
    "-",
    x$interaction_end,
    "\n"
  )
  
  cat(
    "  Missing SCL imputed:",
    x$lengths_imputed,
    "\n"
  )
  
  cat(
    "  Median SCL used:",
    round(x$median_SCL_cm, 2),
    "cm\n"
  )
  
  cat(
    "  Missing mortality imputed:",
    x$mortalities_imputed,
    "\n"
  )
  
  cat(
    "  Mean mortality:",
    round(x$mean_mortality, 4),
    "\n"
  )
  
  cat(
    "  Juveniles:",
    x$juveniles,
    "\n"
  )
  
  cat(
    "  Adults:",
    x$adults,
    "\n"
  )
  
  cat(
    "  Total ANE through trend endpoint:",
    round(
      x$ANE_through_trend_end,
      4
    ),
    "\n"
  )
  
  cat(
    "  Pending/unavailable years:",
    x$pending_years,
    "\n"
  )
  
  cat(
    "  Confirmed zero-interaction years:",
    x$confirmed_zero_years,
    "\n\n"
  )
}

cat("IMPORTANT:\n")
cat("  PENDING years are excluded, not treated as zero.\n")
cat("  Complete zero-interaction years remain valid observed zeros.\n")
cat("  No future take has been applied.\n")
cat("============================================================\n")


# ==============================================================================
# HISTORICAL SSLL TAKE EFFECT ON NESTING TRAJECTORY
# ==============================================================================

if (!requireNamespace("jagsUI", quietly = TRUE)) stop("Package 'jagsUI' is required.")

nest_dir <- cfg$output_dirs$nesting
out_dir <- cfg$output_dirs$update
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

q3 <- function(x) unname(quantile(x, c(.025, .5, .975), na.rm = TRUE))

max_finite_rhat <- function(fit) {
  s <- as.data.frame(fit$summary)
  if (!"Rhat" %in% names(s)) return(NA_real_)
  z <- s$Rhat[is.finite(s$Rhat)]
  if (length(z)) max(z) else NA_real_
}

min_finite_neff <- function(fit) {
  s <- as.data.frame(fit$summary)
  nm <- intersect(c("n.eff", "n_eff"), names(s))
  if (!length(nm)) return(NA_real_)
  z <- s[[nm[1]]][is.finite(s[[nm[1]]])]
  if (length(z)) min(z) else NA_real_
}

# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------

cc_file <- file.path(nest_dir, "cc_annual_nesters.csv")
dc_file <- file.path(out_dir, "update_leatherback_annual_imputed_nests_2001_2025.csv")
cc_ane_file <- file.path(out_dir, "cc_historical_ane_annual_update.csv")
dc_ane_file <- file.path(out_dir, "dc_historical_ane_annual_update.csv")

cc_obs_fit_file <- file.path(out_dir, "update_trend_loggerhead_1986_2025.rds")
dc_med_obs_fit_file <- file.path(out_dir, "update_trend_leatherback_MEDIAN_2001_2025.rds")
dc_low_obs_fit_file <- file.path(out_dir, "update_trend_leatherback_LOW_2001_2025.rds")
dc_high_obs_fit_file <- file.path(out_dir, "update_trend_leatherback_HIGH_2001_2025.rds")

required_files <- c(
  cc_file, dc_file, cc_ane_file, dc_ane_file,
  cc_obs_fit_file, dc_med_obs_fit_file,
  dc_low_obs_fit_file, dc_high_obs_fit_file
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required file(s):\n", paste(missing_files, collapse = "\n"))
}

# ------------------------------------------------------------------------------
# Observed loggerhead nesting
# ------------------------------------------------------------------------------

cc_obs <- read.csv(cc_file, check.names = FALSE)
cc_needed <- c("Year", "Inakahama", "Maehama", "Yotsusehama")

if (!all(cc_needed %in% names(cc_obs))) {
  stop("Loggerhead nesting file missing: ",
       paste(setdiff(cc_needed, names(cc_obs)), collapse = ", "))
}

cc_obs <- cc_obs[
  cc_obs$Year >= 1986 & cc_obs$Year <= 2025,
  cc_needed,
  drop = FALSE
]

# ------------------------------------------------------------------------------
# Observed leatherback nesting branches
# ------------------------------------------------------------------------------

dc_imp <- read.csv(dc_file, check.names = FALSE)

dc_needed <- c(
  "Season",
  "JM_low", "JM_median", "JM_high",
  "W_low", "W_median", "W_high"
)

if (!all(dc_needed %in% names(dc_imp))) {
  stop("Leatherback annual imputation file missing: ",
       paste(setdiff(dc_needed, names(dc_imp)), collapse = ", "))
}

dc_imp <- dc_imp[dc_imp$Season >= 2001 & dc_imp$Season <= 2025, , drop = FALSE]
dc_cf <- cfg$biology$dc$CF

make_dc_branch <- function(branch) {
  out <- data.frame(
    Season = dc_imp$Season,
    JM = dc_imp[[paste0("JM_", branch)]] / dc_cf,
    W = dc_imp[[paste0("W_", branch)]] / dc_cf
  )
  out$W[out$Season < 2006L] <- NA_real_
  out
}

dc_obs_med <- make_dc_branch("median")
dc_obs_low <- make_dc_branch("low")
dc_obs_high <- make_dc_branch("high")

# ------------------------------------------------------------------------------
# Historical ANE
# ------------------------------------------------------------------------------

cc_ane <- read.csv(cc_ane_file, check.names = FALSE)
dc_ane <- read.csv(dc_ane_file, check.names = FALSE)

ane_needed <- c("nesting_year", "ANE")
if (!all(ane_needed %in% names(cc_ane))) stop("Loggerhead ANE file has unexpected columns.")
if (!all(ane_needed %in% names(dc_ane))) stop("Leatherback ANE file has unexpected columns.")

# ------------------------------------------------------------------------------
# Loggerhead: add historical ANE back using observed annual beach proportions
# ------------------------------------------------------------------------------

cc_no_take <- cc_obs
cc_allocation <- data.frame(
  Year = integer(),
  Target_Year = integer(),
  ANE = numeric(),
  allocated = logical(),
  note = character(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(cc_ane))) {
  yr  <- cc_ane$nesting_year[i]
  ane <- cc_ane$ANE[i]
  if (!is.finite(ane) || ane == 0) next
  
  # Re-route 2020 missing year ANE to 2021
  target_yr <- if (yr == 2020L) 2021L else yr
  
  row_id <- match(target_yr, cc_no_take$Year)
  if (is.na(row_id)) next
  
  # Base observed beach values from cc_obs to calculate true site proportions
  beach_values <- as.numeric(
    cc_obs[row_id, c("Inakahama", "Maehama", "Yotsusehama")]
  )
  beach_total <- sum(beach_values, na.rm = TRUE)
  
  # Skip if target year also lacks valid observed nesting data
  if (!is.finite(beach_total) || beach_total <= 0) {
    cc_allocation <- rbind(
      cc_allocation,
      data.frame(
        Year = yr,
        Target_Year = target_yr,
        ANE = ane,
        allocated = FALSE,
        note = "Target year missing beach counts",
        stringsAsFactors = FALSE
      )
    )
    next
  }
  
  beach_prop <- beach_values / beach_total
  observed   <- !is.na(beach_values)
  
  # Pull current values from cc_no_take to accumulate multiple ANE additions (e.g. 2020 + 2021)
  current_vals <- as.numeric(
    cc_no_take[row_id, c("Inakahama", "Maehama", "Yotsusehama")]
  )
  
  adjusted <- current_vals
  adjusted[observed] <- current_vals[observed] + (ane * beach_prop[observed])
  
  cc_no_take[row_id, c("Inakahama", "Maehama", "Yotsusehama")] <- adjusted
  
  note_str <- if (target_yr != yr) paste0("Allocated to ", target_yr) else "Direct allocation"
  
  cc_allocation <- rbind(
    cc_allocation,
    data.frame(
      Year = yr,
      Target_Year = target_yr,
      ANE = ane,
      allocated = TRUE,
      note = note_str,
      stringsAsFactors = FALSE
    )
  )
}

# ------------------------------------------------------------------------------
# Leatherback: use MEDIAN JM/W proportions for all three branches
# ------------------------------------------------------------------------------

dc_no_take_med <- dc_obs_med
dc_no_take_low <- dc_obs_low
dc_no_take_high <- dc_obs_high

median_total <- rowSums(dc_obs_med[c("JM", "W")], na.rm = TRUE)
prop_JM <- dc_obs_med$JM / median_total

dc_allocation <- data.frame(
  Season = integer(), ANE = numeric(),
  prop_JM = numeric(), prop_W = numeric()
)

for (i in seq_len(nrow(dc_ane))) {
  yr <- dc_ane$nesting_year[i]
  ane <- dc_ane$ANE[i]
  if (!is.finite(ane) || ane == 0) next
  
  row_id <- match(yr, dc_obs_med$Season)
  if (is.na(row_id)) next
  
  p_jm <- prop_JM[row_id]
  if (!is.finite(p_jm)) {
    stop("Cannot calculate leatherback median site proportion for season ", yr)
  }
  p_w <- 1 - p_jm
  
  for (obj_name in c("dc_no_take_med", "dc_no_take_low", "dc_no_take_high")) {
    obj <- get(obj_name)
    if (!is.na(obj$JM[row_id])) obj$JM[row_id] <- obj$JM[row_id] + ane * p_jm
    if (!is.na(obj$W[row_id])) obj$W[row_id] <- obj$W[row_id] + ane * p_w
    assign(obj_name, obj)
  }
  
  dc_allocation <- rbind(
    dc_allocation,
    data.frame(Season = yr, ANE = ane, prop_JM = p_jm, prop_W = p_w)
  )
}

write.csv(dc_no_take_med,
          file.path(out_dir, "update_leatherback_MEDIAN_NO_SSLL_TAKE_nesting.csv"),
          row.names = FALSE)
write.csv(dc_no_take_low,
          file.path(out_dir, "update_leatherback_LOW_NO_SSLL_TAKE_nesting.csv"),
          row.names = FALSE)
write.csv(dc_no_take_high,
          file.path(out_dir, "update_leatherback_HIGH_NO_SSLL_TAKE_nesting.csv"),
          row.names = FALSE)
write.csv(dc_allocation,
          file.path(out_dir, "update_leatherback_historical_ANE_allocation.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Frozen Martin singleUQ model
# ------------------------------------------------------------------------------

single_uq_model <- "
model
{
   A[1] <- 0;

   for(j in 2:n.timeseries) {
      A[j] ~ dnorm(a_mean,1/(a_sd^2));
   }

   U ~ dnorm(u_mean,1/(u_sd^2));

   tauQ ~ dgamma(q_alpha,q_beta);
   Q <- 1/tauQ;

   X0 ~ dnorm(x0_mean,1/(x0_sd^2));
   predX[1] <- X0 + U;
   X[1] <- predX[1];

   for(j in 1:n.timeseries) {
      tauR[j] ~ dgamma(r_alpha,r_beta);
      R[j] <- 1/tauR[j];
      predY[j,1] <- Z[j,1] * X[1] + A[j];
      Y[j,1] ~ dnorm(predY[j,1],tauR[j]);
   }

   for(tt in 2:n.yrs) {
      predX[tt] <- X[tt-1] + U;
      X[tt] ~ dnorm(predX[tt],tauQ);

      for(j in 1:n.timeseries) {
         predY[j,tt] <- Z[j,1] * X[tt] + A[j];
         Y[j,tt] ~ dnorm(predY[j,tt],tauR[j]);
      }
   }
}
"

single_uq_file <- file.path(out_dir, "singleUQ_update_historical_take.txt")
writeLines(single_uq_model, single_uq_file)

run_single_uq <- function(thedata, year_col, data_cols, label) {
  years <- as.integer(thedata[[year_col]])
  mat <- as.matrix(thedata[data_cols])
  
  if (any(mat <= 0, na.rm = TRUE)) {
    stop(label, ": non-positive annual-nester value cannot be logged.")
  }
  
  data_mat <- t(log(mat))
  n.yrs <- ncol(data_mat)
  n.timeseries <- nrow(data_mat)
  
  Y <- rbind(data_mat, NA)
  
  whichPop <- rep(1, n.timeseries)
  n.states <- max(whichPop)
  Z <- matrix(0, n.timeseries + 1, n.states + 1)
  Z[n.timeseries + 1, ] <- NA
  Z[, n.states + 1] <- NA
  for (i in seq_along(whichPop)) Z[i, whichPop[i]] <- 1
  
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
  
  n.samples <- 10000
  mcmc.chains <- 2
  mcmc.thin <- 50
  mcmc.burn <- 5000
  samples2Save <- (mcmc.burn + n.samples) * mcmc.thin
  
  cat("\nRunning ", label, "...\n", sep = "")
  cat("  Years: ", min(years), "-", max(years), "\n", sep = "")
  
  set.seed(132)
  
  fit <- jagsUI::jags(
    data = jags_data,
    inits = NULL,
    parameters.to.save = c("A", "U", "Q", "R", "X0", "X"),
    model.file = single_uq_file,
    n.chains = mcmc.chains,
    n.burnin = mcmc.burn * mcmc.thin,
    n.thin = mcmc.thin,
    n.iter = samples2Save,
    DIC = TRUE,
    parallel = TRUE,
    verbose = FALSE
  )
  
  saveRDS(fit, file.path(out_dir, paste0(label, ".rds")))
  
  list(
    label = label,
    fit = fit,
    years = years,
    n.timeseries = n.timeseries,
    max_rhat = max_finite_rhat(fit),
    min_neff = min_finite_neff(fit)
  )
}

# ------------------------------------------------------------------------------
# Run NO-SSLL-TAKE fits
# ------------------------------------------------------------------------------

cc_nt_fit <- run_single_uq(
  cc_no_take, "Year",
  c("Inakahama", "Maehama", "Yotsusehama"),
  "update_trend_loggerhead_NO_SSLL_TAKE_1986_2025"
)

dc_med_nt_fit <- run_single_uq(
  dc_no_take_med, "Season", c("JM", "W"),
  "update_trend_leatherback_MEDIAN_NO_SSLL_TAKE_2001_2025"
)

dc_low_nt_fit <- run_single_uq(
  dc_no_take_low, "Season", c("JM", "W"),
  "update_trend_leatherback_LOW_NO_SSLL_TAKE_2001_2025"
)

dc_high_nt_fit <- run_single_uq(
  dc_no_take_high, "Season", c("JM", "W"),
  "update_trend_leatherback_HIGH_NO_SSLL_TAKE_2001_2025"
)

# ------------------------------------------------------------------------------
# Load frozen OBSERVED fits from scripts 10 and 12
# ------------------------------------------------------------------------------

make_existing_fit <- function(file, years, n.timeseries, label) {
  fit <- readRDS(file)
  list(
    label = label,
    fit = fit,
    years = years,
    n.timeseries = n.timeseries,
    max_rhat = max_finite_rhat(fit),
    min_neff = min_finite_neff(fit)
  )
}

cc_obs_fit <- make_existing_fit(
  cc_obs_fit_file, 1986:2025, 3, "loggerhead_OBSERVED"
)
dc_med_obs_fit <- make_existing_fit(
  dc_med_obs_fit_file, 2001:2025, 2, "leatherback_MEDIAN_OBSERVED"
)
dc_low_obs_fit <- make_existing_fit(
  dc_low_obs_fit_file, 2001:2025, 2, "leatherback_LOW_OBSERVED"
)
dc_high_obs_fit <- make_existing_fit(
  dc_high_obs_fit_file, 2001:2025, 2, "leatherback_HIGH_OBSERVED"
)

# ------------------------------------------------------------------------------
# Martin trend/current-abundance summaries
# ------------------------------------------------------------------------------

summarise_martin_fit <- function(res, RI) {
  fit <- res$fit
  fy <- length(res$years)
  nsim <- 10000
  
  U <- as.numeric(fit$sims.list$U)
  trend_q <- q3(U)
  X <- as.matrix(fit$sims.list$X)
  A <- as.matrix(fit$sims.list$A)
  
  X_len <- nrow(X)
  if (X_len < nsim) stop(res$label, ": fewer than 10000 posterior X draws.")
  
  X_thin <- floor(seq(1, X_len, length.out = nsim))
  X_thin <- unique(pmax(1L, pmin(X_len, X_thin)))
  
  if (length(X_thin) != nsim) {
    X_thin <- as.integer(seq(1, by = X_len / nsim, length.out = nsim))
  }
  
  regional_nesters <- function(year_index) {
    x <- X[X_thin, year_index]
    total <- exp(x)
    if (res$n.timeseries > 1L) {
      for (i in 2:res$n.timeseries) {
        total <- total + exp(x + A[X_thin, i])
      }
    }
    total
  }
  
  N_fym0 <- regional_nesters(fy)
  N_fym1 <- regional_nesters(fy - 1L)
  N_fym2 <- regional_nesters(fy - 2L)
  N_fym3 <- regional_nesters(fy - 3L)
  
  posts <- data.frame(
    U = U[X_thin],
    Q = as.numeric(fit$sims.list$Q)[X_thin],
    N_fym0 = N_fym0,
    N_fym1 = N_fym1,
    N_fym2 = N_fym2,
    N_fym3 = N_fym3
  )
  
  total_draws <- rowSums(
    posts[c("N_fym0", "N_fym1", "N_fym2", "N_fym3")]
  ) * (RI / 4)
  
  abundance <- rbind(
    `N-3` = q3(N_fym3),
    `N-2` = q3(N_fym2),
    `N-1` = q3(N_fym1),
    `N0` = q3(N_fym0),
    Sum = q3(total_draws)
  )
  colnames(abundance) <- c("L95", "Median", "U95")
  
  list(
    r = c(L95 = trend_q[1], Median = trend_q[2], U95 = trend_q[3]),
    abundance = abundance,
    posts = posts,
    max_rhat = res$max_rhat,
    min_neff = res$min_neff
  )
}

cc_obs_sum <- summarise_martin_fit(cc_obs_fit, cfg$biology$cc$RI)
cc_nt_sum <- summarise_martin_fit(cc_nt_fit, cfg$biology$cc$RI)

dc_med_obs_sum <- summarise_martin_fit(dc_med_obs_fit, cfg$biology$dc$RI)
dc_med_nt_sum <- summarise_martin_fit(dc_med_nt_fit, cfg$biology$dc$RI)

dc_low_obs_sum <- summarise_martin_fit(dc_low_obs_fit, cfg$biology$dc$RI)
dc_low_nt_sum <- summarise_martin_fit(dc_low_nt_fit, cfg$biology$dc$RI)

dc_high_obs_sum <- summarise_martin_fit(dc_high_obs_fit, cfg$biology$dc$RI)
dc_high_nt_sum <- summarise_martin_fit(dc_high_nt_fit, cfg$biology$dc$RI)

# ------------------------------------------------------------------------------
# Save posterior outputs
# ------------------------------------------------------------------------------

write.csv(cc_nt_sum$posts,
          file.path(out_dir, "update_posteriors_loggerhead_NO_SSLL_TAKE.csv"),
          row.names = FALSE)
write.csv(dc_med_nt_sum$posts,
          file.path(out_dir, "update_posteriors_leatherback_MEDIAN_NO_SSLL_TAKE.csv"),
          row.names = FALSE)
write.csv(dc_low_nt_sum$posts,
          file.path(out_dir, "update_posteriors_leatherback_LOW_NO_SSLL_TAKE.csv"),
          row.names = FALSE)
write.csv(dc_high_nt_sum$posts,
          file.path(out_dir, "update_posteriors_leatherback_HIGH_NO_SSLL_TAKE.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Comparison table
# ------------------------------------------------------------------------------

make_comparison <- function(species, branch, observed, no_take, provisional = FALSE) {
  data.frame(
    species = species,
    branch = branch,
    observed_r = observed$r["Median"],
    no_take_r = no_take$r["Median"],
    delta_r = no_take$r["Median"] - observed$r["Median"],
    observed_r_L95 = observed$r["L95"],
    observed_r_U95 = observed$r["U95"],
    no_take_r_L95 = no_take$r["L95"],
    no_take_r_U95 = no_take$r["U95"],
    observed_abundance = observed$abundance["Sum", "Median"],
    no_take_abundance = no_take$abundance["Sum", "Median"],
    delta_abundance =
      no_take$abundance["Sum", "Median"] - observed$abundance["Sum", "Median"],
    no_take_abundance_L95 = no_take$abundance["Sum", "L95"],
    no_take_abundance_U95 = no_take$abundance["Sum", "U95"],
    no_take_max_Rhat = no_take$max_rhat,
    provisional = provisional
  )
}

comparison <- rbind(
  make_comparison("Loggerhead", "MAIN", cc_obs_sum, cc_nt_sum, TRUE),
  make_comparison("Leatherback", "MEDIAN", dc_med_obs_sum, dc_med_nt_sum),
  make_comparison("Leatherback", "LOW", dc_low_obs_sum, dc_low_nt_sum),
  make_comparison("Leatherback", "HIGH", dc_high_obs_sum, dc_high_nt_sum)
)

write.csv(
  comparison,
  file.path(out_dir, "update_historical_SSLL_take_effect_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Report
# ------------------------------------------------------------------------------

cc_allocated <- sum(cc_allocation$ANE[cc_allocation$allocated], na.rm = TRUE)
cc_unallocated <- sum(cc_allocation$ANE[!cc_allocation$allocated], na.rm = TRUE)
dc_allocated <- sum(dc_allocation$ANE, na.rm = TRUE)

print_result <- function(label, observed, no_take) {
  cat(label, "\n")
  cat("  OBSERVED r:      ", sprintf("%.4f", observed$r["Median"]), "\n", sep = "")
  cat("  NO-SSLL-TAKE r:  ", sprintf("%.4f", no_take$r["Median"]),
      " (95% CI ", sprintf("%.4f", no_take$r["L95"]), " to ",
      sprintf("%.4f", no_take$r["U95"]), ")\n", sep = "")
  cat("  Delta r:         ",
      sprintf("%+.4f", no_take$r["Median"] - observed$r["Median"]), "\n", sep = "")
  cat("  OBSERVED abundance:     ",
      round(observed$abundance["Sum", "Median"]), "\n", sep = "")
  cat("  NO-SSLL-TAKE abundance: ",
      round(no_take$abundance["Sum", "Median"]),
      " (95% CI ", round(no_take$abundance["Sum", "L95"]), " to ",
      round(no_take$abundance["Sum", "U95"]), ")\n", sep = "")
  cat("  Delta abundance:        ",
      sprintf("%+.1f",
              no_take$abundance["Sum", "Median"] -
                observed$abundance["Sum", "Median"]), "\n", sep = "")
  cat("  NO-SSLL-TAKE max Rhat:  ",
      sprintf("%.4f", no_take$max_rhat), "\n\n", sep = "")
}

cat("\n============================================================\n")
cat("UPDATED HISTORICAL SSLL TAKE EFFECT — MARTIN METHOD\n")
cat("============================================================\n\n")

cat("LOGGERHEAD ANE ALLOCATION\n")
cat("  ANE allocated: ", round(cc_allocated, 4), "\n", sep = "")
cat("  ANE unallocated in fully missing nesting years: ",
    round(cc_unallocated, 4), "\n\n", sep = "")

cat("LEATHERBACK ANE ALLOCATION\n")
cat("  ANE allocated: ", round(dc_allocated, 4), "\n", sep = "")
cat("  Median JM/W proportions used for MEDIAN, LOW and HIGH.\n\n")

print_result("LOGGERHEAD", cc_obs_sum, cc_nt_sum)
print_result("LEATHERBACK MEDIAN", dc_med_obs_sum, dc_med_nt_sum)
print_result("LEATHERBACK LOW", dc_low_obs_sum, dc_low_nt_sum)
print_result("LEATHERBACK HIGH", dc_high_obs_sum, dc_high_nt_sum)

cat("IMPORTANT:\n")
cat("  Pending years were not treated as zero.\n")
cat("  No future take was applied.\n")
cat("============================================================\n")
