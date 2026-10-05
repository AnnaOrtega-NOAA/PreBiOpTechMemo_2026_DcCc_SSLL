# ==============================================================================
# 04_atl.R
# Updated SSLL Anticipated Take Level (ATL)
#
# Updated method:
#   Fit one annual Conway-Maxwell-Poisson (CMP) distribution to complete
#   observed annual SSLL interaction totals, following the later
#   Siders/take_integrated ATL update approach.
#
# CMP likelihood is evaluated directly in log space to avoid the numerical
# approximation switch at mu = 10 in the historical CMP R package.
#
# Martin 2020 segmented SSLL ATL is reconstructed for comparison only.
#
# Loggerhead 2022/2023 are PENDING and excluded, not treated as zero.
# ==============================================================================

source("R/00_config.R")

in_file <- file.path(cfg$output_dirs$audit, "ssll_annual_counts.csv")
out_dir <- cfg$output_dirs$update

if (!file.exists(in_file)) stop("Missing: ", in_file)

dat <- read.csv(in_file, check.names = FALSE, stringsAsFactors = FALSE)

needed <- c("year", "species", "n", "status")
if (!all(needed %in% names(dat))) {
  stop("Missing columns: ",
       paste(setdiff(needed, names(dat)), collapse = ", "))
}

# ------------------------------------------------------------------------------
# Stable CMP functions
#
# Same mu parameterization used by the historical CMP package:
#
#   P(X=x) proportional to (mu^x / x!)^nu
#
# Normalizing constant is evaluated directly rather than switching to the
# asymptotic approximation used by the historical package when mu > 10.
# ------------------------------------------------------------------------------

logspace_sum <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

cmp_log_terms <- function(mu, nu, max_x = 5000L) {
  x <- 0:max_x
  nu * (x * log(mu) - lgamma(x + 1))
}

cmp_logZ <- function(mu, nu, max_x = 5000L) {
  logspace_sum(cmp_log_terms(mu, nu, max_x))
}

cmp_nll <- function(theta, x, max_x = 5000L) {
  mu <- exp(theta[1])
  nu <- exp(theta[2])
  
  if (!is.finite(mu) || !is.finite(nu) ||
      mu <= 0 || nu <= 0) return(1e100)
  
  logZ <- cmp_logZ(mu, nu, max_x)
  
  ll <- nu * (x * log(mu) - lgamma(x + 1)) - logZ
  
  if (any(!is.finite(ll))) return(1e100)
  
  -sum(ll)
}

fit_cmp <- function(x) {
  starts <- list(
    c(log(5), log(0.05)),
    c(log(10), log(0.10)),
    c(log(15), log(0.10)),
    c(log(20), log(0.20))
  )
  
  fits <- lapply(starts, function(start) {
    optim(
      par = start,
      fn = cmp_nll,
      x = x,
      method = "Nelder-Mead",
      control = list(
        maxit = 10000,
        reltol = 1e-12
      )
    )
  })
  
  values <- vapply(fits, function(z) z$value, numeric(1))
  fit <- fits[[which.min(values)]]
  
  if (fit$convergence != 0) {
    stop("CMP optimization failed; convergence code = ",
         fit$convergence)
  }
  
  list(
    fit = fit,
    mu = exp(fit$par[1]),
    nu = exp(fit$par[2]),
    neg_loglik = fit$value,
    start_nll_range = range(values)
  )
}

cmp_pmf <- function(mu, nu, max_x = 5000L) {
  x <- 0:max_x
  lp <- cmp_log_terms(mu, nu, max_x)
  lp <- lp - logspace_sum(lp)
  p <- exp(lp)
  
  keep <- cumsum(p) < (1 - 1e-12)
  last <- min(max(which(keep)) + 1L, length(x))
  
  data.frame(
    ATL = x[1:last],
    probability = p[1:last]
  )
}

convolve_pmf <- function(a, b) {
  out <- convolve(a, rev(b), type = "open")
  out[out < 0] <- 0
  out / sum(out)
}

summarise_pmf <- function(pmf) {
  x <- pmf$ATL
  p <- pmf$probability / sum(pmf$probability)
  cdf <- cumsum(p)
  
  qfun <- function(prob) x[which(cdf >= prob)[1]]
  
  mn <- sum(x * p)
  vr <- sum((x - mn)^2 * p)
  
  c(
    mean = mn,
    sd = sqrt(vr),
    median = qfun(0.5),
    L95 = qfun(0.025),
    U95 = qfun(0.975),
    P0 = p[x == 0]
  )
}

# ------------------------------------------------------------------------------
# Updated annual CMP fits
# ------------------------------------------------------------------------------

fit_species <- function(species_name) {
  z <- dat[dat$species == species_name, , drop = FALSE]
  
  pending <- z$year[toupper(z$status) == "PENDING"]
  
  complete <- z[
    tolower(z$status) == "complete",
    ,
    drop = FALSE
  ]
  
  if (!nrow(complete)) {
    stop(species_name, ": no complete annual counts.")
  }
  
  fit <- fit_cmp(complete$n)
  pmf <- cmp_pmf(fit$mu, fit$nu)
  sm <- summarise_pmf(pmf)
  
  list(
    species = species_name,
    years = complete$year,
    counts = complete$n,
    pending = pending,
    fit = fit,
    pmf = pmf,
    summary = sm
  )
}

cc_new <- fit_species("Loggerhead")
dc_new <- fit_species("Leatherback")

# ------------------------------------------------------------------------------
# Martin 2020 segmented SSLL ATL
#
# Loggerhead:
#   time block 1 + combined time blocks 2/3
#
# Leatherback:
#   time block 1 + time block 2 + time block 3
# ------------------------------------------------------------------------------

cc_old_1 <- cmp_pmf(
  mu = 3.444104,
  nu = 0.06451325
)

cc_old_2 <- cmp_pmf(
  mu = 7.506372e-05,
  nu = 0.01686626
)

cc_old_prob <- convolve_pmf(
  cc_old_1$probability,
  cc_old_2$probability
)

cc_old <- data.frame(
  ATL = 0:(length(cc_old_prob) - 1L),
  probability = cc_old_prob
)

dc_old_1 <- cmp_pmf(
  mu = 2.124568,
  nu = 0.4805365
)

dc_old_2 <- cmp_pmf(
  mu = 2.344938,
  nu = 0.141262
)

dc_old_3 <- cmp_pmf(
  mu = 0.03930914,
  nu = 0.1149811
)

dc_old_prob <- convolve_pmf(
  convolve_pmf(
    dc_old_1$probability,
    dc_old_2$probability
  ),
  dc_old_3$probability
)

dc_old <- data.frame(
  ATL = 0:(length(dc_old_prob) - 1L),
  probability = dc_old_prob
)

cc_old_sum <- summarise_pmf(cc_old)
dc_old_sum <- summarise_pmf(dc_old)

# ------------------------------------------------------------------------------
# QA
# ------------------------------------------------------------------------------

qa <- data.frame(
  species = c("Loggerhead", "Leatherback"),
  sample_mean = c(
    mean(cc_new$counts),
    mean(dc_new$counts)
  ),
  fitted_mean = c(
    cc_new$summary["mean"],
    dc_new$summary["mean"]
  ),
  mean_difference = c(
    cc_new$summary["mean"] - mean(cc_new$counts),
    dc_new$summary["mean"] - mean(dc_new$counts)
  ),
  nll_start_range = c(
    diff(cc_new$fit$start_nll_range),
    diff(dc_new$fit$start_nll_range)
  )
)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

write.csv(
  cc_new$pmf,
  file.path(out_dir, "atl_loggerhead_update_CMP.csv"),
  row.names = FALSE
)

write.csv(
  dc_new$pmf,
  file.path(out_dir, "atl_leatherback_update_CMP.csv"),
  row.names = FALSE
)

write.csv(
  cc_old,
  file.path(out_dir, "atl_loggerhead_Martin2020_CMP.csv"),
  row.names = FALSE
)

write.csv(
  dc_old,
  file.path(out_dir, "atl_leatherback_Martin2020_CMP.csv"),
  row.names = FALSE
)

make_row <- function(species, method, mu, nu, sm, years) {
  data.frame(
    species = species,
    method = method,
    years = years,
    mu = mu,
    nu = nu,
    mean_ATL = sm["mean"],
    sd_ATL = sm["sd"],
    median_ATL = sm["median"],
    L95 = sm["L95"],
    U95 = sm["U95"],
    P_zero = sm["P0"],
    row.names = NULL
  )
}

summary_table <- rbind(
  make_row(
    "Loggerhead",
    "Martin 2020 segmented CMP",
    NA, NA,
    cc_old_sum,
    "published"
  ),
  make_row(
    "Loggerhead",
    "Updated annual CMP",
    cc_new$fit$mu,
    cc_new$fit$nu,
    cc_new$summary,
    paste(range(cc_new$years), collapse = "-")
  ),
  make_row(
    "Leatherback",
    "Martin 2020 segmented CMP",
    NA, NA,
    dc_old_sum,
    "published"
  ),
  make_row(
    "Leatherback",
    "Updated annual CMP",
    dc_new$fit$mu,
    dc_new$fit$nu,
    dc_new$summary,
    paste(range(dc_new$years), collapse = "-")
  )
)

write.csv(
  summary_table,
  file.path(out_dir, "atl_update_summary.csv"),
  row.names = FALSE
)

write.csv(
  qa,
  file.path(out_dir, "atl_update_QA.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# Console report
# ------------------------------------------------------------------------------

print_atl <- function(x, old_summary) {
  cat(toupper(x$species), "\n")
  
  cat("  Complete years used: ",
      min(x$years), "-", max(x$years),
      " (n = ", length(x$years), ")\n", sep = "")
  
  if (length(x$pending)) {
    cat("  Pending years excluded: ",
        paste(x$pending, collapse = ","), "\n", sep = "")
  } else {
    cat("  Pending years excluded: none\n")
  }
  
  cat("  Observed annual counts: mean ",
      round(mean(x$counts), 3),
      "; median ", median(x$counts),
      "; range ", min(x$counts), "-", max(x$counts), "\n", sep = "")
  
  cat("  Updated CMP mu: ",
      round(x$fit$mu, 6), "\n", sep = "")
  
  cat("  Updated CMP nu: ",
      round(x$fit$nu, 6), "\n", sep = "")
  
  cat("  Updated negative log-likelihood: ",
      round(x$fit$neg_loglik, 6), "\n", sep = "")
  
  cat("  NLL range across starting values: ",
      format(diff(x$fit$start_nll_range), scientific = TRUE),
      "\n", sep = "")
  
  cat("  Updated ATL mean: ",
      round(x$summary["mean"], 3), "\n", sep = "")
  
  cat("  Updated ATL SD: ",
      round(x$summary["sd"], 3), "\n", sep = "")
  
  cat("  Updated ATL median: ",
      x$summary["median"], "\n", sep = "")
  
  cat("  Updated ATL 95% interval: ",
      x$summary["L95"], "-", x$summary["U95"],
      "\n", sep = "")
  
  cat("  Updated P(ATL=0): ",
      round(x$summary["P0"], 4), "\n", sep = "")
  
  cat("  Martin ATL mean: ",
      round(old_summary["mean"], 3), "\n", sep = "")
  
  cat("  Martin ATL median: ",
      old_summary["median"], "\n", sep = "")
  
  cat("  Martin ATL 95% interval: ",
      old_summary["L95"], "-", old_summary["U95"],
      "\n\n", sep = "")
}

cat("\n============================================================\n")
cat("UPDATED SSLL ANTICIPATED TAKE LEVEL — STABLE CMP FIT\n")
cat("============================================================\n\n")

print_atl(cc_new, cc_old_sum)
print_atl(dc_new, dc_old_sum)

cat("QA\n")
print(qa, row.names = FALSE, digits = 8)

cat("\nMETHOD:\n")
cat("  Updated ATL = annual CMP fit to complete observed SSLL totals.\n")
cat("  CMP likelihood evaluated directly in log space.\n")
cat("  No asymptotic switch at mu = 10 is used.\n")
cat("  Martin 2020 segmented CMP retained as the comparison ATL.\n")
cat("  Loggerhead 2022/2023 pending years excluded, not zero.\n")
cat("  Loggerhead 2011 retained in the annual updated fit.\n")
cat("  No trend, abundance, historical ANE or PVA method changed.\n")
cat("============================================================\n")