# ==============================================================================
# 06_pva.R
#
# Updated SSLL future take / no-take PVA
# Martin et al. (2020) projection framework with:
#   - updated no-historical-take population posteriors
#   - updated annual-CMP ATL distributions
#   - updated take-demographic fits
#
# Primary projection:
#   10,000 simulations x 100 years
#   stochastic take demographics
#   standard mortality branch (grim_reaper = FALSE)
#
# Loggerhead RI sensitivity:
#   The Martin et al. (2020) loggerhead RI distribution is a Normal
#   (mean = 3.3 y, sd = 2.3 y) truncated only at RI > 0. Because juvenile
#   ANE contains 1 / RI, rare near-zero RI draws create an extreme ANE tail.
#   The primary analysis remains Martin-faithful (RI > 0). A secondary
#   sensitivity constrains loggerhead RI to >= 1 year while holding the exact
#   primary ATL matrix fixed. This sensitivity does not replace the primary.
# ==============================================================================

source(file.path("R", "00_config.R"))

for (pkg in c("mvtnorm", "truncnorm")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required.")
  }
}

out <- cfg$output_dirs$update
dir.create(out, recursive = TRUE, showWarnings = FALSE)

n_sim <- 10000L
n_y   <- 100L
seed  <- 132L

thresholds <- c(0.50, 0.25, 0.125)
horizons   <- c(5L, 10L, 25L, 50L, 100L)


# ==============================================================================
# 1. LOAD UPDATED NO-HISTORICAL-TAKE POSTERIORS
# ==============================================================================

post_files <- c(
  CC = "update_posteriors_loggerhead_NO_SSLL_TAKE.csv",
  DC_MEDIAN = "update_posteriors_leatherback_MEDIAN_NO_SSLL_TAKE.csv",
  DC_LOW = "update_posteriors_leatherback_LOW_NO_SSLL_TAKE.csv",
  DC_HIGH = "update_posteriors_leatherback_HIGH_NO_SSLL_TAKE.csv"
)

post <- lapply(post_files, function(f) {
  p <- file.path(out, f)
  if (!file.exists(p)) stop("Missing 03_historical_take.R output: ", p)
  
  x <- read.csv(p)
  needed <- c("U", "Q", "N_fym0")
  
  if (!all(needed %in% names(x))) {
    stop(f, " is missing: ",
         paste(setdiff(needed, names(x)), collapse = ", "))
  }
  
  x <- x[, needed]
  names(x)[3] <- "N0"
  
  if (nrow(x) < n_sim) {
    stop(f, " contains fewer than ", n_sim, " posterior draws.")
  }
  
  x[seq_len(n_sim), ]
})


# ==============================================================================
# 2. LOAD UPDATED ATL PMFs
# ==============================================================================

read_pmf <- function(file) {
  p <- file.path(out, file)
  if (!file.exists(p)) stop("Missing 04_atl.R output: ", p)
  
  x <- read.csv(p)
  if (!all(c("ATL", "probability") %in% names(x))) {
    stop(file, " must contain ATL and probability.")
  }
  
  x <- x[is.finite(x$ATL) & is.finite(x$probability), ]
  x$probability <- x$probability / sum(x$probability)
  x
}

atl_pmf <- list(
  CC = read_pmf("atl_loggerhead_update_CMP.csv"),
  DC = read_pmf("atl_leatherback_update_CMP.csv")
)


# ==============================================================================
# 3. LOAD UPDATED TAKE-DEMOGRAPHIC PARAMETERS
# ==============================================================================

td_file <- file.path(out, "update_take_demographics_parameters.rds")
if (!file.exists(td_file)) stop("Missing 05_take_demographics.R output: ", td_file)

TD <- readRDS(td_file)

if (!all(c("CC", "DC") %in% names(TD))) {
  stop("Take-demographic RDS must contain CC and DC.")
}

check_td <- function(x, spp) {
  needed <- c("mu0", "beta0", "beta1", "cov")
  if (!all(needed %in% names(x))) {
    stop(spp, " take-demographic object missing: ",
         paste(setdiff(needed, names(x)), collapse = ", "))
  }
}

check_td(TD$CC, "CC")
check_td(TD$DC, "DC")


# ==============================================================================
# 4. BIOLOGICAL PARAMETERS — MARTIN 2020
# ==============================================================================

bio <- list(
  
  CC = list(
    VBGF = list(
      model = "Lknot",
      Linf = 80.4473850,
      K = 0.1396317,
      Lknot = 4.7363329,
      Amat = 26.4950786
    ),
    Pj_mean = 0.80,
    Pj_sd   = 0.031,
    PF      = 0.65,
    RI_type = "normal",
    RI_mean = 3.3,
    RI_sd   = 2.3
  ),
  
  DC = list(
    VBGF = list(
      model = "tknot",
      Linf = 142.7,
      K = 0.2262,
      tknot = -0.17,
      Amat = 16.1
    ),
    Pj_mean = 0.81,
    Pj_sd   = 0.030,
    PF      = 0.73,
    RI_type = "CMP",
    RI_lambda = 17.583686,
    RI_nu     = 2.364034
  )
)


# ==============================================================================
# 5. HELPERS
# ==============================================================================

# Positive leatherback remigration-interval CMP.
# Martin redraws zero values; conditioning the same CMP on RI > 0 is equivalent.
make_positive_cmp <- function(lambda, nu, max_k = 50L) {
  k <- 1:max_k
  lp <- k * log(lambda) - nu * lgamma(k + 1)
  p <- exp(lp - max(lp))
  p <- p / sum(p)
  
  function(n) sample(k, n, replace = TRUE, prob = p)
}

DC_RI <- make_positive_cmp(
  bio$DC$RI_lambda,
  bio$DC$RI_nu
)


# Length-at-age grids.
make_growth_grid <- function(b) {
  age <- seq(0, 100, by = 0.01)
  
  len <- if (b$VBGF$model == "tknot") {
    with(
      b$VBGF,
      Linf * (1 - exp(-K * (age - tknot)))
    )
  } else {
    # Biologically coherent Martin loggerhead VBGF.
    with(
      b$VBGF,
      Linf - (Linf - Lknot) * exp(-K * age)
    )
  }
  
  list(age = age, length = len)
}

growth <- list(
  CC = make_growth_grid(bio$CC),
  DC = make_growth_grid(bio$DC)
)


# Nearest age on Martin's 0.01-y growth grid.
age_from_length <- function(x, g) {
  
  j <- findInterval(x, g$length)
  
  j1 <- pmax(1L, pmin(length(g$length), j))
  j2 <- pmax(1L, pmin(length(g$length), j + 1L))
  
  use2 <- abs(x - g$length[j2]) < abs(x - g$length[j1])
  idx <- ifelse(use2, j2, j1)
  
  g$age[idx]
}


# ==============================================================================
# 6. FUTURE TAKE
# ==============================================================================

simulate_future_take <- function(
    spp,
    pmf,
    td,
    ATL_fixed = NULL,
    RI_min = 0
) {
  
  b <- bio[[spp]]
  g <- growth[[spp]]
  
  if (is.null(ATL_fixed)) {
    ATL <- matrix(
      sample(
        pmf$ATL,
        n_y * n_sim,
        replace = TRUE,
        prob = pmf$probability
      ),
      nrow = n_y,
      ncol = n_sim
    )
  } else {
    ATL <- ATL_fixed

    if (!all(dim(ATL) == c(n_y, n_sim))) {
      stop(
        "ATL_fixed must have dimensions ",
        n_y, " x ", n_sim, "."
      )
    }
  }

  if (spp != "CC" && RI_min != 0) {
    stop("RI_min sensitivity is implemented for loggerheads only.")
  }
  
  ANE <- matrix(0, n_y, n_sim)
  
  cat("\nSimulating ", spp, " future take...\n", sep = "")
  
  for (s in seq_len(n_sim)) {
    
    atl <- ATL[, s]
    n_total <- sum(atl)
    
    if (n_total > 0) {
      
      yr <- rep(seq_len(n_y), atl)
      
      # Martin: mean log-length depends on annual ATL.
      mu_length <- td$beta0 + td$beta1 * atl
      
      z <- mvtnorm::rmvnorm(
        n_total,
        mean = c(0, 0),
        sigma = td$cov
      )
      
      length_cm <- exp(mu_length[yr] + z[, 1])
      mortality <- plogis(td$mu0 + z[, 2])
      
      age <- age_from_length(length_cm, g)
      years_to_maturity <- b$VBGF$Amat - age
      stage <- ifelse(age > b$VBGF$Amat, "A", "J")
      
      RI <- if (b$RI_type == "normal") {
        truncnorm::rtruncnorm(
          n_total,
          a = RI_min,
          mean = b$RI_mean,
          sd = b$RI_sd
        )
      } else {
        DC_RI(n_total)
      }
      
      Pj <- rnorm(n_total, b$Pj_mean, b$Pj_sd)
      
      ANEj <- (Pj ^ years_to_maturity) * (1 / RI)
      stage_ANE <- ifelse(stage == "A", 1, ANEj)
      
      sex  <- rbinom(n_total, 1, b$PF)
      dead <- rbinom(n_total, 1, mortality)
      
      rANE <- stage_ANE * sex * dead
      
      sums <- rowsum(
        rANE,
        group = yr,
        reorder = FALSE
      )
      
      ANE[as.integer(rownames(sums)), s] <- sums[, 1]
    }
    
    if (s %% 1000L == 0L) {
      cat("  ", s, "/", n_sim, "\n", sep = "")
    }
  }
  
  list(ATL = ATL, ANE = ANE)
}


# ==============================================================================
# 7. PVA PROJECTION
# ==============================================================================

project_branch <- function(trend, take_ANE, dynUQ, label) {
  
  N_take    <- matrix(NA_real_, n_y, n_sim)
  N_no_take <- matrix(NA_real_, n_y, n_sim)
  
  cat("\nProjecting ", label, "...\n", sep = "")
  
  for (s in seq_len(n_sim)) {
    
    # Martin:
    #   CC = static U/Q
    #   DC = dynamic annual U/Q
    idx <- if (dynUQ) {
      c(
        s,
        sample(
          seq_len(nrow(trend)),
          n_y - 1L,
          replace = FALSE
        )
      )
    } else {
      rep(s, n_y)
    }
    
    lambda <- exp(trend$U[idx])
    Qsd <- sqrt(trend$Q[idx])
    
    # Independent process-error draws, as in Martin 2020.
    N_no_take[1, s] <- rnorm(
      1,
      trend$N0[s] * lambda[1],
      Qsd[1]
    )
    
    N_take[1, s] <- rnorm(
      1,
      (trend$N0[s] - take_ANE[1, s]) * lambda[1],
      Qsd[1]
    )
    
    for (y in 2:n_y) {
      
      N_no_take[y, s] <- rnorm(
        1,
        N_no_take[y - 1, s] * lambda[y],
        Qsd[y]
      )
      
      N_take[y, s] <- rnorm(
        1,
        (N_take[y - 1, s] - take_ANE[y, s]) * lambda[y],
        Qsd[y]
      )
    }
    
    if (s %% 1000L == 0L) {
      cat("  ", s, "/", n_sim, "\n", sep = "")
    }
  }
  
  list(
    N_take = N_take,
    N_no_take = N_no_take,
    N0 = trend$N0
  )
}


# ==============================================================================
# 8. RUN FUTURE TAKE
# ==============================================================================

set.seed(seed)

# Separate realizations, matching Martin's separate PVA runs.
CC_take     <- simulate_future_take("CC", atl_pmf$CC, TD$CC)

DC_take_med <- simulate_future_take("DC", atl_pmf$DC, TD$DC)
DC_take_low <- simulate_future_take("DC", atl_pmf$DC, TD$DC)
DC_take_hi  <- simulate_future_take("DC", atl_pmf$DC, TD$DC)


# ==============================================================================
# 9. RUN PVA
# ==============================================================================

PVA_CC <- project_branch(
  post$CC,
  CC_take$ANE,
  dynUQ = FALSE,
  "LOGGERHEAD"
)

PVA_DC_MEDIAN <- project_branch(
  post$DC_MEDIAN,
  DC_take_med$ANE,
  dynUQ = TRUE,
  "LEATHERBACK MEDIAN"
)

PVA_DC_LOW <- project_branch(
  post$DC_LOW,
  DC_take_low$ANE,
  dynUQ = TRUE,
  "LEATHERBACK LOW"
)

PVA_DC_HIGH <- project_branch(
  post$DC_HIGH,
  DC_take_hi$ANE,
  dynUQ = TRUE,
  "LEATHERBACK HIGH"
)


# ==============================================================================
# 10. MARTIN-STYLE EXTINCTION HANDLING
# ==============================================================================

fix_extinction <- function(x) {
  
  for (s in seq_len(ncol(x))) {
    neg <- which(x[, s] < 0)
    
    if (length(neg)) {
      x[min(neg):nrow(x), s] <- 0
    }
  }
  
  x
}


# ==============================================================================
# 11. THRESHOLD SUMMARIES
# ==============================================================================

summarise_thresholds <- function(PVA, species, branch) {
  
  mats <- list(
    TAKE = fix_extinction(PVA$N_take),
    NO_TAKE = fix_extinction(PVA$N_no_take)
  )
  
  out <- list()
  
  for (scenario in names(mats)) {
    
    X <- mats[[scenario]]
    
    for (p in thresholds) {
      
      limit <- PVA$N0 * p
      
      first_year <- vapply(
        seq_len(n_sim),
        function(s) {
          z <- which(X[, s] < limit[s])
          if (length(z)) min(z) else NA_integer_
        },
        integer(1)
      )
      
      yrs <- first_year[!is.na(first_year)]
      
      for (h in horizons) {
        
        crossed <- vapply(
          seq_len(n_sim),
          function(s) any(X[seq_len(h), s] < limit[s]),
          logical(1)
        )
        
        out[[length(out) + 1L]] <- data.frame(
          species = species,
          branch = branch,
          threshold = p,
          scenario = scenario,
          horizon = h,
          probability_below = mean(crossed),
          mean_year = if (h == 100L && length(yrs)) mean(yrs) else NA,
          median_year = if (h == 100L && length(yrs)) median(yrs) else NA,
          L95_year = if (h == 100L && length(yrs))
            unname(quantile(yrs, 0.025)) else NA,
          U95_year = if (h == 100L && length(yrs))
            unname(quantile(yrs, 0.975)) else NA
        )
      }
    }
  }
  
  do.call(rbind, out)
}


threshold_results <- rbind(
  summarise_thresholds(PVA_CC, "Loggerhead", "PRIMARY"),
  summarise_thresholds(PVA_DC_MEDIAN, "Leatherback", "MEDIAN"),
  summarise_thresholds(PVA_DC_LOW, "Leatherback", "LOW"),
  summarise_thresholds(PVA_DC_HIGH, "Leatherback", "HIGH")
)


# TAKE - NO_TAKE probability difference.
threshold_differences <- merge(
  subset(
    threshold_results,
    scenario == "TAKE",
    select = c(
      species, branch, threshold,
      horizon, probability_below
    )
  ),
  subset(
    threshold_results,
    scenario == "NO_TAKE",
    select = c(
      species, branch, threshold,
      horizon, probability_below
    )
  ),
  by = c("species", "branch", "threshold", "horizon"),
  suffixes = c("_take", "_no_take")
)

threshold_differences$difference <-
  threshold_differences$probability_below_take -
  threshold_differences$probability_below_no_take


# ==============================================================================
# 12. TRAJECTORY SUMMARIES
# ==============================================================================

summarise_trajectory <- function(PVA, species, branch) {
  
  mats <- list(
    TAKE = fix_extinction(PVA$N_take),
    NO_TAKE = fix_extinction(PVA$N_no_take)
  )
  
  do.call(
    rbind,
    lapply(names(mats), function(scenario) {
      
      X <- mats[[scenario]]
      
      do.call(
        rbind,
        lapply(horizons, function(h) {
          
          q <- quantile(
            X[h, ],
            c(0.025, 0.5, 0.975),
            na.rm = TRUE
          )
          
          data.frame(
            species = species,
            branch = branch,
            scenario = scenario,
            horizon = h,
            L95 = unname(q[1]),
            median = unname(q[2]),
            U95 = unname(q[3])
          )
        })
      )
    })
  )
}

trajectory_summary <- rbind(
  summarise_trajectory(PVA_CC, "Loggerhead", "PRIMARY"),
  summarise_trajectory(PVA_DC_MEDIAN, "Leatherback", "MEDIAN"),
  summarise_trajectory(PVA_DC_LOW, "Leatherback", "LOW"),
  summarise_trajectory(PVA_DC_HIGH, "Leatherback", "HIGH")
)


# ==============================================================================
# 13. FUTURE-TAKE QA
# ==============================================================================

take_qa <- function(x, species, branch) {
  
  atl <- as.vector(x$ATL)
  ane <- as.vector(x$ANE)
  
  data.frame(
    species = species,
    branch = branch,
    ATL_mean = mean(atl),
    ATL_sd = sd(atl),
    ATL_median = median(atl),
    ATL_L95 = unname(quantile(atl, 0.025)),
    ATL_U95 = unname(quantile(atl, 0.975)),
    ANE_mean = mean(ane),
    ANE_sd = sd(ane),
    ANE_median = median(ane),
    ANE_L95 = unname(quantile(ane, 0.025)),
    ANE_U95 = unname(quantile(ane, 0.975))
  )
}

future_take_QA <- rbind(
  take_qa(CC_take, "Loggerhead", "PRIMARY"),
  take_qa(DC_take_med, "Leatherback", "MEDIAN"),
  take_qa(DC_take_low, "Leatherback", "LOW"),
  take_qa(DC_take_hi, "Leatherback", "HIGH")
)


# ==============================================================================
# 14. LOGGERHEAD RI >= 1 YEAR SENSITIVITY
# ==============================================================================
#
# Primary Martin-faithful analysis above:
#   RI ~ Normal(3.3, 2.3), truncated at RI > 0.
#
# Sensitivity:
#   Same model, but loggerhead RI is truncated at RI >= 1 year.
#
# The exact primary loggerhead ATL matrix is reused so the sensitivity is not
# confounded by a different set of annual interaction counts.
#
# TAKE and NO_TAKE process errors remain independently drawn, as in the primary
# Martin projection framework. Therefore very small differences between primary
# and sensitivity threshold probabilities should not be interpreted beyond
# Monte Carlo precision.
# ==============================================================================

set.seed(seed)

CC_take_RImin1 <- simulate_future_take(
  "CC",
  atl_pmf$CC,
  TD$CC,
  ATL_fixed = CC_take$ATL,
  RI_min = 1
)

PVA_CC_RImin1 <- project_branch(
  post$CC,
  CC_take_RImin1$ANE,
  dynUQ = FALSE,
  "LOGGERHEAD RI >= 1"
)


# Annual realized ANE comparison.
ane_sensitivity_summary <- function(x, scenario) {

  ane <- as.vector(x$ANE)

  data.frame(
    scenario = scenario,
    mean = mean(ane),
    sd = sd(ane),
    median = median(ane),
    L95 = unname(quantile(ane, 0.025)),
    U95 = unname(quantile(ane, 0.975)),
    q99 = unname(quantile(ane, 0.99)),
    q999 = unname(quantile(ane, 0.999)),
    maximum = max(ane)
  )
}

loggerhead_RI_ANE_sensitivity <- rbind(
  ane_sensitivity_summary(
    CC_take,
    "PRIMARY_RI_GT_0"
  ),
  ane_sensitivity_summary(
    CC_take_RImin1,
    "SENSITIVITY_RI_GE_1"
  )
)


# Threshold summaries under the RI >= 1 sensitivity.
loggerhead_RI_threshold_results <- summarise_thresholds(
  PVA_CC_RImin1,
  "Loggerhead",
  "RI_MIN_1_SENSITIVITY"
)

loggerhead_RI_threshold_differences <- merge(
  subset(
    loggerhead_RI_threshold_results,
    scenario == "TAKE",
    select = c(
      species, branch, threshold,
      horizon, probability_below
    )
  ),
  subset(
    loggerhead_RI_threshold_results,
    scenario == "NO_TAKE",
    select = c(
      species, branch, threshold,
      horizon, probability_below
    )
  ),
  by = c("species", "branch", "threshold", "horizon"),
  suffixes = c("_take", "_no_take")
)

loggerhead_RI_threshold_differences$difference <-
  loggerhead_RI_threshold_differences$probability_below_take -
  loggerhead_RI_threshold_differences$probability_below_no_take


# Direct comparison of the population-level take effect.
primary_CC_difference <- subset(
  threshold_differences,
  species == "Loggerhead" & branch == "PRIMARY",
  select = c(threshold, horizon, difference)
)

names(primary_CC_difference)[3] <- "primary_RI_GT_0"

RImin1_CC_difference <- loggerhead_RI_threshold_differences[
  ,
  c("threshold", "horizon", "difference")
]

names(RImin1_CC_difference)[3] <- "sensitivity_RI_GE_1"

loggerhead_RI_effect_comparison <- merge(
  primary_CC_difference,
  RImin1_CC_difference,
  by = c("threshold", "horizon")
)

loggerhead_RI_effect_comparison$change_due_to_sensitivity <-
  loggerhead_RI_effect_comparison$sensitivity_RI_GE_1 -
  loggerhead_RI_effect_comparison$primary_RI_GT_0


# ==============================================================================
# 15. SAVE
# ==============================================================================

write.csv(
  threshold_results,
  file.path(out, "update_PVA_threshold_probabilities.csv"),
  row.names = FALSE
)

write.csv(
  threshold_differences,
  file.path(out, "update_PVA_threshold_differences.csv"),
  row.names = FALSE
)

write.csv(
  trajectory_summary,
  file.path(out, "update_PVA_trajectory_summary.csv"),
  row.names = FALSE
)

write.csv(
  future_take_QA,
  file.path(out, "update_PVA_future_take_QA.csv"),
  row.names = FALSE
)

write.csv(
  loggerhead_RI_ANE_sensitivity,
  file.path(out, "update_PVA_loggerhead_RI_sensitivity_ANE.csv"),
  row.names = FALSE
)

write.csv(
  loggerhead_RI_threshold_results,
  file.path(out, "update_PVA_loggerhead_RI_sensitivity_thresholds.csv"),
  row.names = FALSE
)

write.csv(
  loggerhead_RI_effect_comparison,
  file.path(out, "update_PVA_loggerhead_RI_sensitivity_effect.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    RI_lower_bound_years = 1,
    future_take = CC_take_RImin1,
    PVA = PVA_CC_RImin1,
    ANE_summary = loggerhead_RI_ANE_sensitivity,
    threshold_results = loggerhead_RI_threshold_results,
    effect_comparison = loggerhead_RI_effect_comparison
  ),
  file.path(out, "update_PVA_loggerhead_RI_sensitivity.rds")
)

saveRDS(
  list(
    posterior = post$CC,
    future_take = CC_take,
    PVA = PVA_CC
  ),
  file.path(out, "update_PVA_loggerhead.rds")
)

saveRDS(
  list(
    MEDIAN = list(
      posterior = post$DC_MEDIAN,
      future_take = DC_take_med,
      PVA = PVA_DC_MEDIAN
    ),
    LOW = list(
      posterior = post$DC_LOW,
      future_take = DC_take_low,
      PVA = PVA_DC_LOW
    ),
    HIGH = list(
      posterior = post$DC_HIGH,
      future_take = DC_take_hi,
      PVA = PVA_DC_HIGH
    )
  ),
  file.path(out, "update_PVA_leatherback.rds")
)


# ==============================================================================
# 16. PRINT RESULTS
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("UPDATED SSLL TAKE / NO-TAKE PVA\n")
cat("============================================================\n")
cat("Simulations:", n_sim, "\n")
cat("Projection years:", n_y, "\n")
cat("Thresholds: 50%, 25%, 12.5% of starting N0\n")
cat("Horizons: 5, 10, 25, 50, 100 years\n\n")

cat("FUTURE TAKE QA\n")
print(future_take_QA, row.names = FALSE)

cat("\n100-YEAR THRESHOLD RESULTS\n")
print(
  threshold_results[
    threshold_results$horizon == 100,
    c(
      "species", "branch", "threshold", "scenario",
      "probability_below", "mean_year", "median_year",
      "L95_year", "U95_year"
    )
  ],
  row.names = FALSE
)

cat("\nTAKE - NO_TAKE PROBABILITY DIFFERENCES\n")
print(threshold_differences, row.names = FALSE)

cat("\nLOGGERHEAD RI SENSITIVITY - ANNUAL REALIZED ANE\n")
print(loggerhead_RI_ANE_sensitivity, row.names = FALSE)

cat("\nLOGGERHEAD RI SENSITIVITY - TAKE - NO_TAKE EFFECT\n")
print(loggerhead_RI_effect_comparison, row.names = FALSE)

cat("\nSaved to:\n", out, "\n", sep = "")
cat("============================================================\n")