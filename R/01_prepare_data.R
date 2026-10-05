# ==============================================================================
# 01_prepare_data.R
#
# Team-facing data preparation for the 2026 SSLL Technical Memorandum update.
#
# Combines the frozen input-audit step with the nesting preparation required by
# the updated population models. Source CSVs are never modified. No imputation
# or statistical model fitting occurs here.
# ==============================================================================

source("R/00_config.R")

out_dir <- cfg$output_dirs$audit
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# HELPERS
# ==============================================================================

read_input <- function(file) {
  
  path <- file.path(cfg$input_dir, file)
  
  if (!file.exists(path)) {
    stop("Missing input: ", path)
  }
  
  d <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "N/A")
  )
  
  d$.source_row <- seq_len(nrow(d)) + 1L
  d
}


need <- function(d, cols, label) {
  
  missing <- setdiff(cols, names(d))
  
  if (length(missing)) {
    stop(
      label,
      " missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
}


num <- function(x, label) {
  
  out <- suppressWarnings(as.numeric(x))
  
  bad <- is.na(out) & !is.na(x)
  
  if (any(bad)) {
    stop(
      label,
      " contains nonnumeric values: ",
      paste(unique(x[bad]), collapse = ", ")
    )
  }
  
  out
}


write_out <- function(x, file) {
  
  write.csv(
    x,
    file.path(out_dir, file),
    row.names = FALSE,
    na = "NA"
  )
}


# ==============================================================================
# 1. LOGGERHEAD NESTING
# ==============================================================================

cc <- read_input(cfg$files$cc_nesting)

need(
  cc,
  c(
    "Year",
    "Inakahama",
    "Maehama",
    "Yotsusehama"
  ),
  "Loggerhead nesting"
)


for (x in c(
  "Year",
  "Inakahama",
  "Maehama",
  "Yotsusehama"
)) {
  
  cc[[x]] <- num(
    cc[[x]],
    paste("Loggerhead", x)
  )
}


if (anyDuplicated(cc$Year)) {
  stop("Duplicate loggerhead years.")
}


if (any(
  as.matrix(
    cc[c(
      "Inakahama",
      "Maehama",
      "Yotsusehama"
    )]
  ) < 0,
  na.rm = TRUE
)) {
  stop("Negative loggerhead nest count.")
}


cc_clean <- data.frame(
  year = as.integer(cc$Year),
  inakahama = cc$Inakahama,
  maehama = cc$Maehama,
  yotsusehama = cc$Yotsusehama
)


cc_clean <- cc_clean[
  order(cc_clean$year),
]


write_out(
  cc_clean,
  "cc_nesting_prepared.csv"
)


# ==============================================================================
# 2. LEATHERBACK NESTING
# ==============================================================================

dc <- read_input(cfg$files$dc_nesting)

need(
  dc,
  c(
    "Year_begin",
    "Month_begin",
    "JM_Nests",
    "W_Nests"
  ),
  "Leatherback nesting"
)


for (x in c(
  "Year_begin",
  "Month_begin",
  "JM_Nests",
  "W_Nests"
)) {
  
  dc[[x]] <- num(
    dc[[x]],
    paste("Leatherback", x)
  )
}


if (any(!dc$Month_begin %in% 1:12)) {
  stop("Leatherback Month_begin must be 1-12.")
}


if (any(
  as.matrix(
    dc[c(
      "JM_Nests",
      "W_Nests"
    )]
  ) < 0,
  na.rm = TRUE
)) {
  stop("Negative leatherback nest count.")
}


# ------------------------------------------------------------------------------
# Merge complementary duplicate year-month rows.
# Never sum duplicate observations.
# ------------------------------------------------------------------------------

groups <- split(
  seq_len(nrow(dc)),
  paste(
    dc$Year_begin,
    dc$Month_begin,
    sep = "-"
  )
)


dc_clean <- do.call(
  rbind,
  lapply(groups, function(i) {
    
    b <- dc[i, , drop = FALSE]
    
    get_value <- function(x) {
      
      z <- unique(x[!is.na(x)])
      
      if (length(z) > 1) {
        stop(
          "Conflicting leatherback values for ",
          b$Year_begin[1],
          "-",
          b$Month_begin[1]
        )
      }
      
      if (length(z)) z else NA_real_
    }
    
    
    data.frame(
      year = as.integer(b$Year_begin[1]),
      month = as.integer(b$Month_begin[1]),
      jm_nests = get_value(b$JM_Nests),
      w_nests = get_value(b$W_Nests)
    )
  })
)


rownames(dc_clean) <- NULL

dc_clean <- dc_clean[
  order(
    dc_clean$year,
    dc_clean$month
  ),
]


write_out(
  dc_clean,
  "dc_nesting_monthly_prepared.csv"
)


# ==============================================================================
# 3. PIRO SSLL INTERACTIONS
#
# Exact PIRO headers are used below, including:
#
#   "Length  (SCL cm)"              <- two spaces after Length
#   "Mortality Coefficent (low)"    <- source spelling
#   "Mortality Coefficent (high)"   <- source spelling
#
# The source CSV is never modified.
# ==============================================================================

ssll <- read_input(
  cfg$files$interactions
)


need(
  ssll,
  c(
    "Year",
    "Species",
    "Length  (SCL cm)",
    "Mortality Coefficent (low)",
    "Mortality Coefficent (high)"
  ),
  "PIRO SSLL"
)


# ------------------------------------------------------------------------------
# Convert analytical fields to numeric
# ------------------------------------------------------------------------------

ssll$Year <- num(
  ssll$Year,
  "SSLL Year"
)


ssll[["Length  (SCL cm)"]] <- num(
  ssll[["Length  (SCL cm)"]],
  "SSLL SCL cm"
)


ssll[["Mortality Coefficent (low)"]] <- num(
  ssll[["Mortality Coefficent (low)"]],
  "SSLL mortality low"
)


ssll[["Mortality Coefficent (high)"]] <- num(
  ssll[["Mortality Coefficent (high)"]],
  "SSLL mortality high"
)


# ------------------------------------------------------------------------------
# Validate interaction years
# ------------------------------------------------------------------------------

if (any(
  ssll$Year < cfg$update$interaction_start |
  ssll$Year > cfg$update$interaction_end,
  na.rm = TRUE
)) {
  
  stop(
    "SSLL interaction year outside configured ",
    cfg$update$interaction_start,
    "-",
    cfg$update$interaction_end,
    " range."
  )
}


# ------------------------------------------------------------------------------
# Standardize species
# ------------------------------------------------------------------------------

s <- tolower(
  trimws(ssll$Species)
)


species <- ifelse(
  
  s %in% c(
    "loggerhead",
    "loggerhead turtle",
    "caretta caretta",
    "cc"
  ),
  
  "Loggerhead",
  
  ifelse(
    
    s %in% c(
      "leatherback",
      "leatherback turtle",
      "dermochelys coriacea",
      "dc"
    ),
    
    "Leatherback",
    
    NA_character_
  )
)


if (anyNA(species)) {
  
  stop(
    "Unrecognized species: ",
    paste(
      unique(ssll$Species[is.na(species)]),
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# SCL
#
# Use Length (SCL cm) only.
# No ft-to-cm conversion is performed.
# ------------------------------------------------------------------------------

scl_cm <- ssll[[
  "Length  (SCL cm)"
]]


if (any(
  scl_cm <= 0,
  na.rm = TRUE
)) {
  stop("Nonpositive SCL value.")
}


# ==============================================================================
# 4. MORTALITY
#
# PIRO's supplied "mean" column is NOT used.
#
# Analytical mortality:
#
#   both bounds present -> mean(low, high)
#   one bound present   -> available value
#   both missing        -> NA for later Martin-method imputation
#
# One source m_high value is known to exceed the logical probability bound.
# The original value is retained, while the analytical value is capped at 1.
# ==============================================================================

m_low <- ssll[[
  "Mortality Coefficent (low)"
]]


m_high <- ssll[[
  "Mortality Coefficent (high)"
]]


# Preserve original source values before any analytical correction.
m_high_original <- m_high


# Low mortality values must already be valid probabilities.
if (any(
  m_low < 0 | m_low > 1,
  na.rm = TRUE
)) {
  stop("Mortality low outside 0-1.")
}


# Negative high values are invalid.
if (any(
  m_high < 0,
  na.rm = TRUE
)) {
  stop("Negative mortality high value.")
}


# Flag and cap high mortality values above 1.
mortality_capped <- !is.na(m_high) & m_high > 1


if (any(mortality_capped)) {
  
  warning(
    sum(mortality_capped),
    " mortality high value(s) exceeded 1 and were capped at 1. ",
    "Original values retained in m_high_original."
  )
  
  m_high[
    mortality_capped
  ] <- 1
}


# Check low/high relationship after correction.
if (any(
  m_low > m_high,
  na.rm = TRUE
)) {
  stop("Mortality low exceeds mortality high.")
}


# Calculate mean mortality following the Martin approach.
m_mean <- rowMeans(
  cbind(
    m_low,
    m_high
  ),
  na.rm = TRUE
)


# Both missing produces NaN; retain as NA for later imputation.
m_mean[
  is.nan(m_mean)
] <- NA_real_


# ==============================================================================
# 5. CLEAN INTERACTION TABLE
# ==============================================================================

ssll_clean <- data.frame(
  
  year = as.integer(ssll$Year),
  
  species = species,
  
  scl_cm = scl_cm,
  
  m_low = m_low,
  
  m_high = m_high,
  
  m_high_original = m_high_original,
  
  mortality_capped = mortality_capped,
  
  m_mean = m_mean,
  
  source_row = ssll$.source_row
)


ssll_clean <- ssll_clean[
  order(
    ssll_clean$year,
    ssll_clean$species,
    ssll_clean$source_row
  ),
]


rownames(ssll_clean) <- NULL


write_out(
  ssll_clean,
  "ssll_interactions_prepared.csv"
)


# ==============================================================================
# 6. MARTIN-METHOD BASELINE FISHERY SUBSET
#
# Current PIRO ground-truth dataset restricted to 2005-2018.
#
# This is NOT being presented as Martin's original fishery dataset.
# ==============================================================================

baseline <- ssll_clean[
  
  ssll_clean$year >=
    cfg$baseline$interaction_start &
    
    ssll_clean$year <=
    cfg$baseline$interaction_end,
  
]


if (!nrow(baseline)) {
  stop(
    "No SSLL records found in baseline 2005-2018 window."
  )
}


write_out(
  baseline,
  "ssll_interactions_baseline_2005_2018.csv"
)


# ==============================================================================
# 7. SIMPLE INTERACTION AUDIT
# ==============================================================================

audit <- do.call(
  rbind,
  lapply(
    c(
      "Loggerhead",
      "Leatherback"
    ),
    function(spp) {
      
      x <- ssll_clean[
        ssll_clean$species == spp,
      ]
      
      b <- baseline[
        baseline$species == spp,
      ]
      
      
      data.frame(
        
        species = spp,
        
        rows_all =
          nrow(x),
        
        rows_2005_2018 =
          nrow(b),
        
        first_year =
          if (nrow(x)) min(x$year) else NA_integer_,
        
        last_year =
          if (nrow(x)) max(x$year) else NA_integer_,
        
        missing_scl =
          sum(is.na(x$scl_cm)),
        
        missing_mortality =
          sum(is.na(x$m_mean)),
        
        mortality_capped =
          sum(x$mortality_capped)
      )
    }
  )
)


write_out(
  audit,
  "interaction_audit_summary.csv"
)


# ==============================================================================
# 8. ANNUAL INTERACTION COUNTS
#
# Inspection only.
#
# Loggerhead 2022 and 2023 are known to be pending and must not later be
# interpreted as zero interactions.
# ==============================================================================

count_years <- seq(
  cfg$update$interaction_start,
  cfg$update$interaction_end
)


counts <- expand.grid(
  year = count_years,
  species = c(
    "Loggerhead",
    "Leatherback"
  ),
  stringsAsFactors = FALSE
)


counts$n <- mapply(
  function(y, spp) {
    
    sum(
      ssll_clean$year == y &
        ssll_clean$species == spp
    )
    
  },
  counts$year,
  counts$species
)


counts$status <- "complete"


counts$status[
  counts$species == "Loggerhead" &
    counts$year %in%
    cfg$data_status$cc_interaction_years_pending
] <- "PENDING"


write_out(
  counts,
  "ssll_annual_counts.csv"
)


# ==============================================================================
# 9. REPORT
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("INPUT AUDIT COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Loggerhead nesting: ",
  min(cc_clean$year),
  "-",
  max(cc_clean$year),
  "\n",
  sep = ""
)


cat(
  "Leatherback nesting: ",
  min(dc_clean$year),
  "-",
  max(dc_clean$year),
  "\n\n",
  sep = ""
)


cat("SSLL interaction audit:\n")

print(
  audit,
  row.names = FALSE
)


cat("\nAnnual SSLL interaction counts:\n")

print(
  counts,
  row.names = FALSE
)


cat(
  "\nBaseline fishery subset: ",
  cfg$baseline$interaction_start,
  "-",
  cfg$baseline$interaction_end,
  "\n",
  sep = ""
)


cat(
  "  Loggerhead: ",
  sum(
    baseline$species == "Loggerhead"
  ),
  " records\n",
  sep = ""
)


cat(
  "  Leatherback: ",
  sum(
    baseline$species == "Leatherback"
  ),
  " records\n",
  sep = ""
)


if (any(mortality_capped)) {
  
  cat("\nMortality correction:\n")
  
  print(
    ssll_clean[
      ssll_clean$mortality_capped,
      c(
        "source_row",
        "year",
        "species",
        "m_low",
        "m_high_original",
        "m_high",
        "m_mean"
      )
    ],
    row.names = FALSE
  )
}


cat("\nRules applied:\n")
cat("  - Used Length  (SCL cm) only.\n")
cat("  - No ft-to-cm conversion performed.\n")
cat("  - Calculated m_mean from low/high mortality.\n")
cat("  - PIRO mean column was not used.\n")
cat("  - Mortality high >1 capped at 1; original retained.\n")
cat("  - Missing SCL and mortality were NOT imputed.\n")
cat("  - Loggerhead 2022 and 2023 remain PENDING.\n")
cat("  - No statistical models were run.\n")

cat("\n============================================================\n")


# ==============================================================================
# 10. PREPARE PRODUCTION NESTING INPUTS
#
# The audited tables above use standardized lower-case field names. This section
# creates the exact nesting files consumed by the frozen updated models.
# ==============================================================================

nest_dir <- cfg$output_dirs$nesting
dir.create(nest_dir, recursive = TRUE, showWarnings = FALSE)

write_nesting <- function(x, file) {
  write.csv(x, file.path(nest_dir, file), row.names = FALSE, na = "NA")
}

# ------------------------------------------------------------------------------
# LOGGERHEAD: annual nesters, 1986-2025
# ------------------------------------------------------------------------------

cc_prep <- read.csv(
  file.path(cfg$output_dirs$audit, "cc_nesting_prepared.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cc_required <- c("year", "inakahama", "maehama", "yotsusehama")
if (!all(cc_required %in% names(cc_prep))) {
  stop("Audited loggerhead file is missing required columns.")
}

cc_grid <- data.frame(
  Year = seq.int(cfg$baseline$cc_nesting_start, cfg$update$cc_nesting_end)
)

cc_for_merge <- data.frame(
  Year = as.integer(cc_prep$year),
  Inakahama = cc_prep$inakahama,
  Maehama = cc_prep$maehama,
  Yotsusehama = cc_prep$yotsusehama
)

cc_grid <- merge(
  cc_grid,
  cc_for_merge,
  by = "Year",
  all.x = TRUE,
  sort = TRUE
)

cc_beaches <- c("Inakahama", "Maehama", "Yotsusehama")

if (!identical(cc_grid$Year, 1986:2025)) {
  stop("Loggerhead production grid must contain exactly 1986-2025.")
}

if (any(as.matrix(cc_grid[cc_beaches]) < 0, na.rm = TRUE)) {
  stop("Negative loggerhead nest count in production grid.")
}

write_nesting(cc_grid, "cc_annual_nests.csv")

cc_nesters <- cc_grid
cc_nesters[cc_beaches] <- cc_nesters[cc_beaches] / cfg$biology$cc$CF
write_nesting(cc_nesters, "cc_annual_nesters.csv")

saveRDS(
  list(
    years = cc_nesters$Year,
    annual_nesters = t(as.matrix(cc_nesters[cc_beaches])),
    reference_beach = "Inakahama",
    CF = cfg$biology$cc$CF
  ),
  file.path(nest_dir, "cc_trend_input.rds")
)

# ------------------------------------------------------------------------------
# LEATHERBACK: complete April-March monthly grid through last observed March
# ------------------------------------------------------------------------------

dc_prep <- read.csv(
  file.path(cfg$output_dirs$audit, "dc_nesting_monthly_prepared.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dc_required <- c("year", "month", "jm_nests", "w_nests")
if (!all(dc_required %in% names(dc_prep))) {
  stop("Audited leatherback file is missing required columns.")
}

dc_prep$date <- as.Date(sprintf("%04d-%02d-01", dc_prep$year, dc_prep$month))
if (anyNA(dc_prep$date) || anyDuplicated(dc_prep$date)) {
  stop("Invalid or duplicate audited leatherback month.")
}

march_dates <- dc_prep$date[dc_prep$month == 3L]
if (!length(march_dates)) stop("No March observations available to close a leatherback season.")

last_march <- max(march_dates)
last_dc_season <- as.integer(format(last_march, "%Y")) - 1L

first_date <- as.Date(sprintf("%04d-04-01", cfg$baseline$dc_nesting_start))
grid <- data.frame(date = seq(first_date, last_march, by = "month"))

dc_merge <- data.frame(
  date = dc_prep$date,
  JM_Nests = dc_prep$jm_nests,
  W_Nests = dc_prep$w_nests
)

grid <- merge(grid, dc_merge, by = "date", all.x = TRUE, sort = TRUE)
grid$Year <- as.integer(format(grid$date, "%Y"))
grid$Month <- as.integer(format(grid$date, "%m"))
grid$Season <- ifelse(grid$Month < 4L, grid$Year - 1L, grid$Year)
grid$Seq_month <- ifelse(grid$Month < 4L, grid$Month + 9L, grid$Month - 3L)
grid$row_present_in_source <- grid$date %in% dc_prep$date
grid$beyond_source_calendar <- grid$date > max(dc_prep$date)

# Martin site window: Wermon begins in season 2006.
grid$W_Nests[grid$Season < 2006L] <- NA_real_

season_counts <- table(grid$Season)
if (any(season_counts != 12L)) {
  stop("Leatherback production grid must contain 12 months per season.")
}

write_nesting(grid, "dc_monthly_season_grid.csv")

season_audit <- do.call(
  rbind,
  lapply(sort(unique(grid$Season)), function(s) {
    b <- grid[grid$Season == s, ]
    data.frame(
      Season = s,
      calendar_start = as.character(min(b$date)),
      calendar_end = as.character(max(b$date)),
      JM_observed_months = sum(!is.na(b$JM_Nests)),
      W_observed_months = sum(!is.na(b$W_Nests)),
      JM_zero_months = sum(b$JM_Nests == 0, na.rm = TRUE),
      W_zero_months = sum(b$W_Nests == 0, na.rm = TRUE),
      months_beyond_source_calendar = sum(b$beyond_source_calendar)
    )
  })
)

write_nesting(season_audit, "dc_season_coverage.csv")

saveRDS(
  list(
    monthly = grid,
    seasons = sort(unique(grid$Season)),
    period = c(JM = 12L, W = 6L),
    notes = "Calendar dates converted to April-March seasons. No imputation yet."
  ),
  file.path(nest_dir, "dc_imputation_input.rds")
)

cat("\nPRODUCTION NESTING INPUTS PREPARED\n")
cat("  Loggerhead years: 1986-2025; CF =", cfg$biology$cc$CF, "\n")
cat("  Leatherback seasons: 2001-", last_dc_season, "\n", sep = "")
cat("  No nesting values were imputed in this script.\n")
