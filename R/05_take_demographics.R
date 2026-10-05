# ==============================================================================
# 05_take_demographics.R
# Updated SSLL future-take demographics
#
# Martin et al. (2020) bivariate model:
#   log(SCL) and logit(post-release mortality)
#
# Mean log(SCL) varies with annual realized take.
# Mean logit(mortality) is constant.
#
# Updated PIRO SSLL interactions:
#   Loggerhead: complete years 2005-2025; 2022/2023 excluded as pending
#   Leatherback: complete years 2005-2025
#
# These fitted demographic distributions feed the future-take PVA.
# ==============================================================================

source("R/00_config.R")

if (!requireNamespace("rstan", quietly = TRUE)) {
  stop("Package 'rstan' is required.")
}
if (!requireNamespace("boot", quietly = TRUE)) {
  stop("Package 'boot' is required.")
}
if (!requireNamespace("truncnorm", quietly = TRUE)) {
  stop("Package 'truncnorm' is required.")
}

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

set.seed(132)

int_file <- file.path(
  cfg$output_dirs$audit,
  "ssll_interactions_prepared.csv"
)

annual_file <- file.path(
  cfg$output_dirs$audit,
  "ssll_annual_counts.csv"
)

out_dir <- cfg$output_dirs$update

if (!file.exists(int_file)) stop("Missing: ", int_file)
if (!file.exists(annual_file)) stop("Missing: ", annual_file)

td <- read.csv(
  int_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

annual <- read.csv(
  annual_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

needed_td <- c(
  "year", "species", "scl_cm",
  "m_low", "m_high"
)

needed_annual <- c(
  "year", "species", "n", "status"
)

if (!all(needed_td %in% names(td))) {
  stop(
    "Interaction file missing: ",
    paste(setdiff(needed_td, names(td)), collapse = ", ")
  )
}

if (!all(needed_annual %in% names(annual))) {
  stop(
    "Annual-count file missing: ",
    paste(setdiff(needed_annual, names(annual)), collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# Martin take-demographics Stan model
# ------------------------------------------------------------------------------

mvnorm <- "
functions {
  matrix cov_matrix_2d(vector sigma, real rho) {
    matrix[2,2] Sigma;

    Sigma[1,1] = square(sigma[1]);
    Sigma[2,2] = square(sigma[2]);
    Sigma[1,2] = sigma[1] * sigma[2] * rho;
    Sigma[2,1] = Sigma[1,2];

    return Sigma;
  }
}

data {
  int<lower=1> N;
  array[N] vector[2] x;

  int<lower=1> nyear;
  array[N] int<lower=1, upper=nyear> year;

  vector[nyear] rtl;
}

parameters {
  real<lower=-1, upper=1> rho;
  vector<lower=0>[2] sigma;

  real beta0;
  real beta1;
  real mu0;
}

transformed parameters {
  array[nyear] vector[2] mu;

  for (y in 1:nyear) {
    mu[y,1] = beta0 + beta1 * rtl[y];
    mu[y,2] = mu0;
  }
}

model {
  mu0  ~ normal(0, 2);
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 2);
  sigma ~ normal(0, 2);

  (rho + 1) / 2 ~ beta(2, 2);

  for (n in 1:N) {
    x[n] ~ multi_normal(
      mu[year[n]],
      cov_matrix_2d(sigma, rho)
    );
  }
}
"

stan_mod <- rstan::stan_model(model_code = mvnorm)

# ------------------------------------------------------------------------------
# Prepare one species
# ------------------------------------------------------------------------------

prepare_species <- function(species_name) {
  
  a <- annual[
    annual$species == species_name,
    ,
    drop = FALSE
  ]
  
  pending_years <- a$year[tolower(a$status) == "pending"]
  
  a <- a[
    tolower(a$status) == "complete",
    ,
    drop = FALSE
  ]
  
  a <- a[order(a$year), ]
  
  years <- a$year
  
  z <- td[
    td$species == species_name &
      td$year %in% years,
    ,
    drop = FALSE
  ]
  
  z <- z[order(z$year), ]
  
  if (!nrow(z)) {
    stop(species_name, ": no interaction records.")
  }
  
  # --------------------------------------------------------------------------
  # Length — Martin species-median imputation
  # --------------------------------------------------------------------------
  
  n_missing_scl <- sum(is.na(z$scl_cm))
  
  scl_median <- median(
    z$scl_cm,
    na.rm = TRUE
  )
  
  z$scl_cm[is.na(z$scl_cm)] <- scl_median
  
  # --------------------------------------------------------------------------
  # Mortality — Martin midpoint and species-mean imputation
  # --------------------------------------------------------------------------
  
  z$M_mu <- rowMeans(
    cbind(z$m_low, z$m_high),
    na.rm = TRUE
  )
  
  z$M_mu[is.nan(z$M_mu)] <- NA_real_
  
  n_missing_mort <- sum(is.na(z$M_mu))
  
  mortality_mean <- mean(
    z$M_mu,
    na.rm = TRUE
  )
  
  z$M_mu[is.na(z$M_mu)] <- mortality_mean
  
  # Martin explicitly constrained loggerhead mortality = 1 to 0.999
  n_one_adjusted <- 0L
  
  if (species_name == "Loggerhead") {
    n_one_adjusted <- sum(z$M_mu == 1)
    
    z$M_mu[z$M_mu == 1] <- 0.999
  }
  
  if (any(!is.finite(z$scl_cm)) || any(z$scl_cm <= 0)) {
    stop(species_name, ": invalid SCL after preparation.")
  }
  
  if (any(!is.finite(z$M_mu)) ||
      any(z$M_mu <= 0 | z$M_mu >= 1)) {
    stop(
      species_name,
      ": mortality outside open interval (0,1) after preparation."
    )
  }
  
  # --------------------------------------------------------------------------
  # Explicit calendar-year indexing
  #
  # This preserves complete zero-interaction years in rtl and avoids the
  # factor-index shift that can occur when a year has zero individual rows.
  # --------------------------------------------------------------------------
  
  rtl <- a$n
  
  year_pointer <- match(
    z$year,
    years
  )
  
  if (anyNA(year_pointer)) {
    stop(species_name, ": failed calendar-year matching.")
  }
  
  stan_data <- list(
    N = nrow(z),
    
    x = cbind(
      log(z$scl_cm),
      boot::logit(z$M_mu)
    ),
    
    nyear = length(years),
    
    year = as.integer(
      year_pointer
    ),
    
    rtl = as.numeric(
      rtl
    )
  )
  
  list(
    data = z,
    annual = a,
    years = years,
    pending_years = pending_years,
    stan_data = stan_data,
    n_missing_scl = n_missing_scl,
    scl_median = scl_median,
    n_missing_mort = n_missing_mort,
    mortality_mean = mortality_mean,
    n_one_adjusted = n_one_adjusted
  )
}

cc <- prepare_species("Loggerhead")
dc <- prepare_species("Leatherback")

# ------------------------------------------------------------------------------
# Martin-style initial values
# ------------------------------------------------------------------------------

make_init <- function(prep) {
  
  function(chain_id) {
    
    d <- prep$data
    x <- prep$stan_data$x
    
    lm0 <- lm(
      log(scl_cm) ~ year,
      data = d
    )
    
    b0 <- coef(lm0)[1]
    b1 <- coef(lm0)[2]
    
    sx <- apply(x, 2, sd)
    rx <- cor(x)[1,2]
    
    if (!is.finite(rx)) rx <- 0
    
    list(
      beta0 = rnorm(
        1,
        b0,
        max(abs(b0) * 0.1, 0.01)
      ),
      
      beta1 = rnorm(
        1,
        b1,
        max(abs(b1) * 0.1, 0.001)
      ),
      
      mu0 = rnorm(
        1,
        0,
        1
      ),
      
      sigma = truncnorm::rtruncnorm(
        2,
        a = 0,
        mean = sx,
        sd = pmax(sx * 0.1, 0.001)
      ),
      
      rho = truncnorm::rtruncnorm(
        1,
        a = -1,
        b = 1,
        mean = rx,
        sd = max(abs(rx) * 0.1, 0.01)
      )
    )
  }
}

# ------------------------------------------------------------------------------
# Fit
# ------------------------------------------------------------------------------

cat("\nFitting LOGGERHEAD take-demographics model...\n")

CC.td.fit <- rstan::sampling(
  stan_mod,
  data = cc$stan_data,
  init = make_init(cc),
  chains = 4,
  iter = 7500,
  warmup = 5000,
  seed = 132,
  refresh = 500
)

cat("\nFitting LEATHERBACK take-demographics model...\n")

DC.td.fit <- rstan::sampling(
  stan_mod,
  data = dc$stan_data,
  init = make_init(dc),
  chains = 4,
  iter = 7500,
  warmup = 5000,
  seed = 132,
  refresh = 500
)

# ------------------------------------------------------------------------------
# Extract Martin parameters
# ------------------------------------------------------------------------------

CC.sims <- rstan::extract(CC.td.fit)
DC.sims <- rstan::extract(DC.td.fit)

CC.mu0 <- median(CC.sims$mu0)
CC.beta0 <- median(CC.sims$beta0)
CC.beta1 <- median(CC.sims$beta1)
CC.sigma <- apply(CC.sims$sigma, 2, median)
CC.rho <- median(CC.sims$rho)

CC.cov <- matrix(
  c(
    CC.sigma[1]^2,
    CC.sigma[1] * CC.sigma[2] * CC.rho,
    CC.sigma[1] * CC.sigma[2] * CC.rho,
    CC.sigma[2]^2
  ),
  2, 2
)

DC.mu0 <- median(DC.sims$mu0)
DC.beta0 <- median(DC.sims$beta0)
DC.beta1 <- median(DC.sims$beta1)
DC.sigma <- apply(DC.sims$sigma, 2, median)
DC.rho <- median(DC.sims$rho)

DC.cov <- matrix(
  c(
    DC.sigma[1]^2,
    DC.sigma[1] * DC.sigma[2] * DC.rho,
    DC.sigma[1] * DC.sigma[2] * DC.rho,
    DC.sigma[2]^2
  ),
  2, 2
)

TD_CC_MVN <- list(
  mu0 = CC.mu0,
  beta0 = CC.beta0,
  beta1 = CC.beta1,
  cov = CC.cov
)

TD_DC_MVN <- list(
  mu0 = DC.mu0,
  beta0 = DC.beta0,
  beta1 = DC.beta1,
  cov = DC.cov
)

# ------------------------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------------------------

max_rhat <- function(fit) {

  draws <- rstan::extract(
    fit,
    permuted = FALSE, 
    inc_warmup = FALSE
  )
  
  mon <- rstan::monitor(
    draws,
    print = FALSE
  )
  
  rh <- mon[, "Rhat"]
  
  max(rh[is.finite(rh)], na.rm = TRUE)
}

CC.rhat <- max_rhat(CC.td.fit)
DC.rhat <- max_rhat(DC.td.fit)

summary_table <- data.frame(
  species = c(
    "Loggerhead",
    "Leatherback"
  ),
  
  n_interactions = c(
    nrow(cc$data),
    nrow(dc$data)
  ),
  
  n_years = c(
    length(cc$years),
    length(dc$years)
  ),
  
  scl_median_imputation = c(
    cc$scl_median,
    dc$scl_median
  ),
  
  n_scl_imputed = c(
    cc$n_missing_scl,
    dc$n_missing_scl
  ),
  
  mortality_mean_imputation = c(
    cc$mortality_mean,
    dc$mortality_mean
  ),
  
  n_mortality_imputed = c(
    cc$n_missing_mort,
    dc$n_missing_mort
  ),
  
  mu0 = c(
    CC.mu0,
    DC.mu0
  ),
  
  beta0 = c(
    CC.beta0,
    DC.beta0
  ),
  
  beta1 = c(
    CC.beta1,
    DC.beta1
  ),
  
  sigma_length = c(
    CC.sigma[1],
    DC.sigma[1]
  ),
  
  sigma_mortality = c(
    CC.sigma[2],
    DC.sigma[2]
  ),
  
  rho = c(
    CC.rho,
    DC.rho
  ),
  
  max_Rhat = c(
    CC.rhat,
    DC.rhat
  )
)

write.csv(
  summary_table,
  file.path(
    out_dir,
    "update_take_demographics_summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    CC = TD_CC_MVN,
    DC = TD_DC_MVN,
    
    metadata = list(
      CC_years = cc$years,
      DC_years = dc$years,
      CC_pending_years = cc$pending_years,
      DC_pending_years = dc$pending_years
    )
  ),
  file.path(
    out_dir,
    "update_take_demographics_parameters.rds"
  )
)

saveRDS(
  CC.td.fit,
  file.path(
    out_dir,
    "update_take_demographics_loggerhead_fit.rds"
  )
)

saveRDS(
  DC.td.fit,
  file.path(
    out_dir,
    "update_take_demographics_leatherback_fit.rds"
  )
)

# ------------------------------------------------------------------------------
# Console report
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("UPDATED SSLL TAKE DEMOGRAPHICS — MARTIN 2020 METHOD\n")
cat("============================================================\n\n")

cat("LOGGERHEAD\n")
cat("  Complete years: ",
    min(cc$years), "-", max(cc$years),
    " (n = ", length(cc$years), ")\n", sep = "")

cat("  Pending years excluded: ",
    if (length(cc$pending_years))
      paste(cc$pending_years, collapse = ",")
    else "none",
    "\n", sep = "")

cat("  Interactions used: ", nrow(cc$data), "\n", sep = "")
cat("  Missing SCL imputed: ", cc$n_missing_scl, "\n", sep = "")
cat("  Median SCL: ", round(cc$scl_median, 3), " cm\n", sep = "")
cat("  Missing mortality imputed: ", cc$n_missing_mort, "\n", sep = "")
cat("  Mean mortality: ", round(cc$mortality_mean, 4), "\n", sep = "")
cat("  Mortality = 1 adjusted to 0.999: ",
    cc$n_one_adjusted, "\n", sep = "")

cat("  mu0: ", round(CC.mu0, 6), "\n", sep = "")
cat("  beta0: ", round(CC.beta0, 6), "\n", sep = "")
cat("  beta1: ", round(CC.beta1, 6), "\n", sep = "")
cat("  sigma length: ", round(CC.sigma[1], 6), "\n", sep = "")
cat("  sigma mortality: ", round(CC.sigma[2], 6), "\n", sep = "")
cat("  rho: ", round(CC.rho, 6), "\n", sep = "")
cat("  max Rhat: ", round(CC.rhat, 6), "\n\n", sep = "")

cat("LEATHERBACK\n")
cat("  Complete years: ",
    min(dc$years), "-", max(dc$years),
    " (n = ", length(dc$years), ")\n", sep = "")

cat("  Pending years excluded: none\n")
cat("  Interactions used: ", nrow(dc$data), "\n", sep = "")
cat("  Missing SCL imputed: ", dc$n_missing_scl, "\n", sep = "")
cat("  Median SCL: ", round(dc$scl_median, 3), " cm\n", sep = "")
cat("  Missing mortality imputed: ", dc$n_missing_mort, "\n", sep = "")
cat("  Mean mortality: ", round(dc$mortality_mean, 4), "\n", sep = "")

cat("  mu0: ", round(DC.mu0, 6), "\n", sep = "")
cat("  beta0: ", round(DC.beta0, 6), "\n", sep = "")
cat("  beta1: ", round(DC.beta1, 6), "\n", sep = "")
cat("  sigma length: ", round(DC.sigma[1], 6), "\n", sep = "")
cat("  sigma mortality: ", round(DC.sigma[2], 6), "\n", sep = "")
cat("  rho: ", round(DC.rho, 6), "\n", sep = "")
cat("  max Rhat: ", round(DC.rhat, 6), "\n\n", sep = "")

cat("Saved:\n")
cat("  update_take_demographics_summary.csv\n")
cat("  update_take_demographics_parameters.rds\n")
cat("  update_take_demographics_loggerhead_fit.rds\n")
cat("  update_take_demographics_leatherback_fit.rds\n")

cat("\nThese parameters are the future-take demographic inputs for the PVA.\n")
cat("============================================================\n")