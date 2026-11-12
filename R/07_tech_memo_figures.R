# ==============================================================================
# 07_tech_memo_figures.R
#
# 2026 SSLL TECHNICAL MEMORANDUM
#
# NO MODELS ARE FITTED HERE.
# EXISTING ACCEPTED OUTPUTS ARE READ ONLY.
# FIGURES ARE PNG ONLY.
#
# FIGURE 1 — NESTING DATA USED IN THE 2026 UPDATE
#   A. North Pacific loggerhead annual nest counts
#   B. West Pacific leatherback annual imputed nest counts
#
# FIGURE 2 — UPDATED POPULATION TREND MODEL FITS
#   A. North Pacific loggerhead
#   B. West Pacific leatherback
#
# FIGURE 3 — 2020 ASSESSMENT VS 2026 POSTERIORS
#   A. Loggerhead trend
#   B. Loggerhead abundance
#   C. Leatherback trend
#   D. Leatherback abundance
#
# FIGURE 4 — SSLL INTERACTIONS AND ANTICIPATED TAKE
#   A. Loggerhead annual interactions
#   B. Loggerhead ATL
#   C. Leatherback annual interactions
#   D. Leatherback ATL
#
# FIGURE 5 — 10-YEAR FUTURE POPULATION PROJECTIONS (95% CrI)
#   A. North Pacific loggerhead (10-year projection)
#   B. West Pacific leatherback (10-year projection)
#
# TABLE 1  — PVA THRESHOLD PROBABILITIES (STATUS QUO)
# TABLE 1B — MULTI-SCENARIO RISK COMPARISON MATRIX
#
# FIGURE 6 — SCL AND POST-INTERACTION MORTALITY DISTRIBUTIONS
#   A. Loggerhead SCL
#   B. Loggerhead mortality coefficient
#   C. Leatherback SCL
#   D. Leatherback mortality coefficient
#
# FIGURE 7 — COUNTERFACTUAL TAKE AND ANE DISTRIBUTIONS
#   A. Loggerhead annual interactions
#   B. Loggerhead Adult Nester Equivalents (ANE)
#   C. Leatherback annual interactions
#   D. Leatherback Adult Nester Equivalents (ANE)
#
# RUN FROM PROJECT ROOT:
#   source("R/07_tech_memo_figures.R")
# ==============================================================================


# ==============================================================================
# 0. SETUP
# ==============================================================================

source("R/00_config.R")

options(stringsAsFactors = FALSE)

required_packages <- c(
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

library(ggplot2)
library(patchwork)


# ==============================================================================
# 1. OUTPUT DIRECTORY
# ==============================================================================

fig_dir <- file.path(
  "outputs",
  "tech_memo_figures"
)

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("\n")
cat("============================================================\n")
cat("2026 TECH MEMO FIGURES\n")
cat("============================================================\n")
cat("PNG only.\n")
cat("No models will be fitted.\n")
cat("Output directory: ", fig_dir, "\n", sep = "")
cat("============================================================\n\n")


# ==============================================================================
# 2. FIND EXISTING OUTPUT
# ==============================================================================

find_output <- function(
    filename,
    required = TRUE
) {
  
  hits <- list.files(
    path = "outputs",
    recursive = TRUE,
    full.names = TRUE
  )
  
  hits <- hits[
    basename(hits) == filename
  ]
  
  if (length(hits) == 0L) {
    if (required) {
      stop(
        "\nCould not find existing analysis output:\n  ",
        filename,
        "\n\nSearched recursively under outputs/.\n",
        "No model was rerun."
      )
    }
    return(NA_character_)
  }
  
  if (length(hits) > 1L) {
    cat(
      "NOTE: multiple copies found for ",
      filename,
      ".\nUsing:\n  ",
      hits[1],
      "\n\n",
      sep = ""
    )
  }
  
  normalizePath(
    hits[1],
    winslash = "/",
    mustWork = TRUE
  )
}


# ==============================================================================
# 3. COLORS
# ==============================================================================

col_cc <- "#C56A2D"
col_dc <- "#2F6B8A"

col_take    <- "#B2473E"   # Status quo future take (rust red)
col_no_take <- "#3F6485"   # No future take (slate blue)
col_no_trip <- "#E69F00"   # No trip limits (amber / gold)
col_no_cap  <- "#7B4173"   # No hard caps (plum / purple)
col_both    <- "#009E73"   # Both removed (teal / green)

col_dark     <- "#262626"
col_mid      <- "#707070"
col_previous <- "#777777"


# ==============================================================================
# 4. COMMON THEME
# ==============================================================================

theme_memo <- function(base_size = 10.5) {
  
  theme_classic(
    base_size = base_size,
    base_family = "sans"
  ) +
    
    theme(
      axis.title = element_text(
        size = base_size + 0.5,
        colour = col_dark
      ),
      axis.text = element_text(
        size = base_size - 0.5,
        colour = col_dark
      ),
      axis.line = element_line(
        linewidth = 0.4,
        colour = col_dark
      ),
      axis.ticks = element_line(
        linewidth = 0.35,
        colour = col_dark
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(
        size = base_size - 0.5
      ),
      plot.title = element_text(
        face = "bold",
        size = base_size + 0.5
      ),
      plot.tag = element_text(
        face = "bold",
        size = base_size + 1.5
      ),
      strip.background = element_blank(),
      strip.text = element_text(
        face = "bold",
        size = base_size + 0.5
      ),
      plot.margin = margin(
        8,
        10,
        8,
        8
      )
    )
}


# ==============================================================================
# 5. SAVE FIGURE — PNG ONLY
# ==============================================================================

save_memo_figure <- function(
    plot,
    filename,
    width,
    height
) {
  
  png_file <- file.path(
    fig_dir,
    paste0(filename, ".png")
  )
  
  ggsave(
    filename = png_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 400,
    bg = "white"
  )
  
  cat(
    "Saved: ",
    png_file,
    "\n",
    sep = ""
  )
}


# ==============================================================================
# 6. GENERAL HELPER
# ==============================================================================

q3 <- function(x) {
  unname(
    quantile(
      x,
      probs = c(
        0.025,
        0.50,
        0.975
      ),
      na.rm = TRUE
    )
  )
}


# ##############################################################################
# FIGURE 1: NESTING DATA USED IN THE 2026 UPDATE
# ##############################################################################

# ==============================================================================
# 7. FIGURE 1A — LOGGERHEAD ANNUAL NEST COUNTS
# ==============================================================================

cc_nests_file <- find_output("cc_annual_nests.csv")
cc_nests      <- read.csv(cc_nests_file, check.names = FALSE)
cc_beaches    <- c("Inakahama", "Maehama", "Yotsusehama")
needed        <- c("Year", cc_beaches)

if (!all(needed %in% names(cc_nests))) {
  stop(
    "Loggerhead annual nest file missing column(s): ",
    paste(setdiff(needed, names(cc_nests)), collapse = ", ")
  )
}

cc_nests_long <- do.call(
  rbind,
  lapply(
    cc_beaches,
    function(site) {
      data.frame(
        Year = cc_nests$Year,
        Site = site,
        Nests = cc_nests[[site]]
      )
    }
  )
)

cc_nests_long$Site <- factor(cc_nests_long$Site, levels = cc_beaches)
cc_nests_long$Used_in_trend <- cc_nests_long$Year != 2021

cc_nests_long$Line_segment <- NA_character_
cc_nests_long$Line_segment[cc_nests_long$Year <= 2020] <- "pre_gap"
cc_nests_long$Line_segment[cc_nests_long$Year >= 2022] <- "post_gap"

cc_line_data <- cc_nests_long[
  !is.na(cc_nests_long$Line_segment) & !is.na(cc_nests_long$Nests),
]

# p1a <- ggplot() +
#   geom_line(
#     data = cc_line_data,
#     aes(
#       x = Year,
#       y = Nests,
#       group = interaction(Site, Line_segment),
#       linetype = Site
#     ),
#     colour = "grey35",
#     linewidth = 0.65,
#     na.rm = TRUE
#   ) +
p1a <- ggplot() +
  # Add subtle grey rectangle for extended years behind the data:
  annotate(
    "rect",
    xmin = 2015.5, xmax = 2025.5,
    ymin = -Inf, ymax = Inf,
    fill = "grey92", colour = NA
  ) +
  geom_line(
    data = cc_line_data,
    aes(x = Year, y = Nests, group = interaction(Site, Line_segment), linetype = Site),
    colour = "grey35", linewidth = 0.65, na.rm = TRUE
  ) +
  geom_point(
    data = cc_nests_long[cc_nests_long$Used_in_trend, ],
    aes(x = Year, y = Nests, shape = Site),
    colour = "black",
    size = 1.5,
    na.rm = TRUE
  ) +
  geom_point(
    data = cc_nests_long[cc_nests_long$Year == 2021, ],
    aes(x = Year, y = Nests),
    shape = 21,
    fill = "white",
    colour = "black",
    size = 2.1,
    stroke = 0.75,
    na.rm = TRUE
  ) +
  scale_linetype_manual(
    values = c(
      "Inakahama" = "solid",
      "Maehama" = "dashed",
      "Yotsusehama" = "dotted"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Inakahama" = 16,
      "Maehama" = 17,
      "Yotsusehama" = 15
    )
  ) +
  scale_x_continuous(breaks = seq(1990, 2025, 5), minor_breaks = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "North Pacific loggerhead",
    x = "Year",
    y = "Annual nest count",
    linetype = NULL,
    shape = NULL
  ) +
  theme_memo() +
  theme(legend.position = "bottom")

# ==============================================================================
# 8. FIGURE 1B — LEATHERBACK ANNUAL IMPUTED NEST COUNTS
# ==============================================================================

dc_nests_file <- find_output("update_leatherback_annual_imputed_nests_2001_2025.csv")
dc_nests      <- read.csv(dc_nests_file)
needed        <- c("Season", "JM_median", "W_median")

if (!all(needed %in% names(dc_nests))) {
  stop(
    "Leatherback annual imputed nest file missing column(s): ",
    paste(setdiff(needed, names(dc_nests)), collapse = ", ")
  )
}

dc_nests_long <- rbind(
  data.frame(Season = dc_nests$Season, Site = "Jeen Yessa", Nests = dc_nests$JM_median),
  data.frame(Season = dc_nests$Season, Site = "Jeen Syuab", Nests = dc_nests$W_median)
)

dc_nests_long$Nests[dc_nests_long$Site == "Jeen Syuab" & dc_nests_long$Season < 2006] <- NA_real_
dc_nests_long$Site <- factor(dc_nests_long$Site, levels = c("Jeen Yessa", "Jeen Syuab"))

p1b <- ggplot(
  data = dc_nests_long,
  aes(x = Season, y = Nests, group = Site)
) +
  # 1. Background shading for extended years added since 2020 assessment (seasons 2018-2025)
  annotate(
    "rect",
    xmin = 2017.5, xmax = 2025.5,
    ymin = -Inf, ymax = Inf,
    fill = "grey92", colour = NA
  ) +
  # 2. Line layer
  geom_line(
    aes(linetype = Site),
    colour = "grey35", linewidth = 0.7, na.rm = TRUE
  ) +
  # 3. Point layer (with shape mapped to Site)
  geom_point(
    aes(shape = Site),
    colour = "black", size = 1.6, na.rm = TRUE
  ) +
  scale_linetype_manual(values = c("Jeen Yessa" = "solid", "Jeen Syuab" = "dashed")) +
  scale_shape_manual(values = c("Jeen Yessa" = 16, "Jeen Syuab" = 17)) +
  scale_x_continuous(breaks = seq(2001, 2025, 3), minor_breaks = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "West Pacific leatherback",
    x = "Nesting season",
    y = "Annual nest count",
    linetype = NULL,
    shape = NULL
  ) +
  theme_memo() +
  theme(legend.position = "bottom")

# ==============================================================================
# 9. COMBINE FIGURE 1
# ==============================================================================

figure1 <- (p1a / p1b) + plot_annotation(tag_levels = "A")

save_memo_figure(
  figure1,
  "Figure_01_nesting_data",
  width = 7.5,
  height = 7.2
)


# ##############################################################################
# FIGURE 2: UPDATED POPULATION TREND MODEL FITS
# ##############################################################################

# ==============================================================================
# 10. TREND FIT DATA HELPER
# ==============================================================================

make_trend_data <- function(fit, observed, observed_years, species) {
  
  X  <- as.matrix(fit$sims.list$X)
  X0 <- as.numeric(fit$sims.list$X0)
  A  <- as.matrix(fit$sims.list$A)
  
  if (ncol(X) != length(observed_years)) {
    stop(species, ": number of modeled X years does not match observed years.")
  }
  if (nrow(A) != nrow(X)) {
    stop(species, ": A and X posterior draws do not align.")
  }
  
  n_draw <- nrow(X)
  n_site <- ncol(A)
  
  X0_total <- rep(0, n_draw)
  for (j in seq_len(n_site)) {
    X0_total <- X0_total + exp(X0 + A[, j])
  }
  
  X_total <- matrix(0, nrow = n_draw, ncol = ncol(X))
  for (j in seq_len(n_site)) {
    X_total <- X_total + exp(X + matrix(A[, j], nrow = n_draw, ncol = ncol(X)))
  }
  
  all_states <- cbind(X0_total, X_total)
  all_years  <- c(min(observed_years) - 1L, observed_years)
  log_states <- log(all_states)
  
  qq <- apply(log_states, 2, quantile, probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  
  modeled <- data.frame(
    Year = all_years,
    L95 = qq[1, ],
    Median = qq[2, ],
    U95 = qq[3, ]
  )
  
  observed_df <- data.frame(
    Year = observed_years,
    LogAnnualNesters = ifelse(observed > 0, log(observed), NA_real_)
  )
  
  list(modeled = modeled, observed = observed_df)
}


# ==============================================================================
# 11. CLEAN TREND PLOT
# ==============================================================================

plot_trend_fit <- function(trend_data, species, x_breaks) {
  
  model <- trend_data$modeled
  obs   <- trend_data$observed
  
  ggplot() +
    geom_ribbon(
      data = model,
      aes(x = Year, ymin = L95, ymax = U95),
      fill = "grey85",
      colour = NA
    ) +
    geom_line(
      data = model,
      aes(x = Year, y = Median),
      colour = "black",
      linewidth = 0.9
    ) +
    geom_point(
      data = obs,
      aes(x = Year, y = LogAnnualNesters),
      colour = "black",
      shape = 16,
      size = 1.7,
      na.rm = TRUE
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      minor_breaks = NULL,
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    labs(
      title = species,
      x = "Year",
      y = "ln(Annual Nesters)"
    ) +
    theme_memo() +
    theme(legend.position = "none")
}


# ==============================================================================
# 12. LOGGERHEAD OBSERVED ANNUAL NESTERS
# ==============================================================================

cc_nesters_file <- find_output("cc_annual_nesters.csv")
cc_nesters      <- read.csv(cc_nesters_file, check.names = FALSE)
cc_nesters      <- cc_nesters[cc_nesters$Year >= 1986 & cc_nesters$Year <= 2025, ]

cc_observed <- rowSums(cc_nesters[, cc_beaches], na.rm = TRUE)
cc_all_missing <- rowSums(!is.na(cc_nesters[, cc_beaches])) == 0
cc_observed[cc_all_missing] <- NA_real_


# ==============================================================================
# 13. LOGGERHEAD TREND FIT
# ==============================================================================

cc_fit_file <- find_output("update_trend_loggerhead_1986_2025.rds")
cc_fit      <- readRDS(cc_fit_file)

cc_trend <- make_trend_data(
  fit = cc_fit,
  observed = cc_observed,
  observed_years = cc_nesters$Year,
  species = "North Pacific loggerhead"
)

p2a <- plot_trend_fit(
  trend_data = cc_trend,
  species = "North Pacific loggerhead",
  x_breaks = seq(1985, 2025, 5)
)


# ==============================================================================
# 14. LEATHERBACK OBSERVED ANNUAL NESTERS
# ==============================================================================

CF_dc <- cfg$biology$dc$CF

dc_trend_input <- data.frame(
  Season = dc_nests$Season,
  JM = dc_nests$JM_median / CF_dc,
  W  = dc_nests$W_median / CF_dc
)
dc_trend_input$W[dc_trend_input$Season < 2006] <- NA_real_

dc_observed <- rowSums(dc_trend_input[, c("JM", "W")], na.rm = TRUE)


# ==============================================================================
# 15. LEATHERBACK TREND FIT
# ==============================================================================

dc_fit_file <- find_output("update_trend_leatherback_MEDIAN_2001_2025.rds")
dc_fit      <- readRDS(dc_fit_file)

dc_trend <- make_trend_data(
  fit = dc_fit,
  observed = dc_observed,
  observed_years = dc_trend_input$Season,
  species = "West Pacific leatherback"
)

p2b <- plot_trend_fit(
  trend_data = dc_trend,
  species = "West Pacific leatherback",
  x_breaks = seq(2000, 2025, 4)
)


# ==============================================================================
# 16. COMBINE FIGURE 2
# ==============================================================================

figure2 <- (p2a / p2b) + plot_annotation(tag_levels = "A")

save_memo_figure(
  figure2,
  "Figure_02_updated_population_trend_fits",
  width = 7.5,
  height = 7.5
)


# ##############################################################################
# FIGURE 3: 2020 ASSESSMENT VS 2026 POSTERIORS
# ##############################################################################

# ==============================================================================
# 17. RECONSTRUCT VALIDATION POSTERIOR IF CSV ABSENT
# ==============================================================================

posterior_from_validation_fit <- function(rds_file, years, n_timeseries, label) {
  
  cat("Reconstructing posterior draws from existing saved fit:\n  ", rds_file, "\n", sep = "")
  fit <- readRDS(rds_file)
  
  if (is.null(fit$sims.list$U) || is.null(fit$sims.list$X) || is.null(fit$sims.list$A)) {
    stop(label, ": saved validation fit does not contain expected posterior objects.")
  }
  
  nsim  <- 10000L
  U     <- as.numeric(fit$sims.list$U)
  X     <- as.matrix(fit$sims.list$X)
  A     <- as.matrix(fit$sims.list$A)
  X_len <- nrow(X)
  
  if (X_len < nsim) stop(label, ": saved validation fit has fewer than 10,000 posterior draws.")
  
  X_thin <- floor(seq(from = 1, to = X_len, length.out = nsim))
  X_thin <- unique(pmax(1L, pmin(X_len, X_thin)))
  
  if (length(X_thin) != nsim) {
    step   <- X_len / nsim
    X_thin <- as.integer(seq(from = 1, by = step, length.out = nsim))
  }
  
  fy <- length(years)
  regional_nesters <- function(year_index) {
    x <- X[X_thin, year_index]
    total <- exp(x)
    if (n_timeseries > 1L) {
      for (i in 2:n_timeseries) {
        total <- total + exp(x + A[X_thin, i])
      }
    }
    total
  }
  
  N_fym0 <- regional_nesters(fy)
  N_fym1 <- regional_nesters(fy - 1L)
  N_fym2 <- regional_nesters(fy - 2L)
  N_fym3 <- regional_nesters(fy - 3L)
  
  Q <- if (!is.null(fit$sims.list$Q)) as.numeric(fit$sims.list$Q)[X_thin] else rep(NA_real_, nsim)
  
  data.frame(
    U = U[X_thin],
    Q = Q,
    N_fym0 = N_fym0,
    N_fym1 = N_fym1,
    N_fym2 = N_fym2,
    N_fym3 = N_fym3
  )
}


# ==============================================================================
# 18. 2020 ASSESSMENT LOGGERHEAD POSTERIOR
# ==============================================================================

cc_previous_csv <- find_output("validation_posteriors_loggerhead.csv", required = FALSE)

if (!is.na(cc_previous_csv)) {
  cc_previous <- read.csv(cc_previous_csv)
} else {
  cc_previous_rds <- find_output("validation_trend_loggerhead_1986_2015.rds", required = FALSE)
  if (is.na(cc_previous_rds)) stop("Cannot make Figure 3: loggerhead validation posterior not found.")
  cc_previous <- posterior_from_validation_fit(
    rds_file = cc_previous_rds, years = 1986:2015, n_timeseries = 3L,
    label = "Loggerhead 2020 assessment"
  )
}


# ==============================================================================
# 19. 2020 ASSESSMENT LEATHERBACK POSTERIOR
# ==============================================================================

dc_previous_csv <- find_output("validation_posteriors_leatherback_MEDIAN.csv", required = FALSE)

if (!is.na(dc_previous_csv)) {
  dc_previous <- read.csv(dc_previous_csv)
} else {
  dc_previous_rds <- find_output("validation_trend_leatherback_MEDIAN_2001_2017.rds", required = FALSE)
  if (is.na(dc_previous_rds)) stop("Cannot make Figure 3: leatherback validation posterior not found.")
  dc_previous <- posterior_from_validation_fit(
    rds_file = dc_previous_rds, years = 2001:2017, n_timeseries = 2L,
    label = "Leatherback 2020 assessment"
  )
}


# ==============================================================================
# 20. UPDATED POSTERIORS
# ==============================================================================

cc_updated <- read.csv(find_output("update_posteriors_loggerhead_1986_2025.csv"))
dc_updated <- read.csv(find_output("update_posteriors_leatherback_MEDIAN_2001_2025.csv"))


# ==============================================================================
# 21. CHECK POSTERIOR STRUCTURE
# ==============================================================================

required_posterior_columns <- c("U", "N_fym0", "N_fym1", "N_fym2", "N_fym3")

check_posterior <- function(x, label) {
  missing_columns <- setdiff(required_posterior_columns, names(x))
  if (length(missing_columns) > 0L) {
    stop(label, " is missing posterior column(s): ", paste(missing_columns, collapse = ", "))
  }
  invisible(TRUE)
}

check_posterior(cc_previous, "Loggerhead 2020 assessment")
check_posterior(cc_updated,  "Loggerhead 2026 update")
check_posterior(dc_previous, "Leatherback 2020 assessment")
check_posterior(dc_updated,  "Leatherback 2026 update")


# ==============================================================================
# 22. CONVERT POSTERIORS
# ==============================================================================

make_posterior_data <- function(x, species, assessment, period, RI) {
  annual_percent    <- (exp(x$U) - 1) * 100
  current_abundance <- rowSums(x[, c("N_fym0", "N_fym1", "N_fym2", "N_fym3")]) * RI / 4
  
  rbind(
    data.frame(
      Species = species, Assessment = assessment, Period = period,
      Quantity = "Annual population trend", Value = annual_percent
    ),
    data.frame(
      Species = species, Assessment = assessment, Period = period,
      Quantity = "Current abundance", Value = current_abundance
    )
  )
}

posterior_compare <- rbind(
  make_posterior_data(cc_previous, "North Pacific loggerhead", "2020 assessment", "1986-2015", cfg$biology$cc$RI),
  make_posterior_data(cc_updated,  "North Pacific loggerhead", "2026 update",      "1986-2025", cfg$biology$cc$RI),
  make_posterior_data(dc_previous, "West Pacific leatherback", "2020 assessment", "2001-2017", cfg$biology$dc$RI),
  make_posterior_data(dc_updated,  "West Pacific leatherback", "2026 update",      "2001-2025", cfg$biology$dc$RI)
)

posterior_compare$Assessment <- factor(
  posterior_compare$Assessment,
  levels = c("2020 assessment", "2026 update")
)


# ==============================================================================
# 23. POSTERIOR QA SUMMARY
# ==============================================================================

posterior_groups <- split(
  posterior_compare,
  interaction(posterior_compare$Species, posterior_compare$Quantity, posterior_compare$Assessment, drop = TRUE)
)

posterior_summary <- do.call(
  rbind,
  lapply(posterior_groups, function(x) {
    qq <- q3(x$Value)
    data.frame(
      Species    = x$Species[1],
      Assessment = as.character(x$Assessment[1]),
      Period     = x$Period[1],
      Quantity   = x$Quantity[1],
      L95        = qq[1],
      Median     = qq[2],
      U95        = qq[3],
      n_draws    = nrow(x)
    )
  })
)
rownames(posterior_summary) <- NULL

write.csv(
  posterior_summary,
  file.path(fig_dir, "posterior_comparison_QA.csv"),
  row.names = FALSE
)


# ==============================================================================
# 24. POSTERIOR DENSITY FUNCTION
# ==============================================================================

posterior_density_plot <- function(data, species, quantity, update_colour, x_label, panel_title, zero_line = FALSE) {
  
  dat     <- data[data$Species == species & data$Quantity == quantity, ]
  medians <- aggregate(Value ~ Assessment, data = dat, FUN = median)
  
  p <- ggplot(dat, aes(x = Value, fill = Assessment, colour = Assessment, linetype = Assessment)) +
    geom_density(alpha = 0.22, linewidth = 0.9, adjust = 1, trim = TRUE) +
    geom_vline(
      data = medians,
      aes(xintercept = Value, colour = Assessment, linetype = Assessment),
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c("2020 assessment" = "#BDBDBD", "2026 update" = update_colour)) +
    scale_colour_manual(values = c("2020 assessment" = col_previous, "2026 update" = update_colour)) +
    scale_linetype_manual(values = c("2020 assessment" = "dashed", "2026 update" = "solid")) +
    labs(title = panel_title, x = x_label, y = "Posterior density") +
    theme_memo()
  
  if (zero_line) {
    p <- p + geom_vline(xintercept = 0, linewidth = 0.4, colour = "#333333")
  }
  p
}


# ==============================================================================
# 25. FIGURE 3 PANELS
# ==============================================================================

p3a <- posterior_density_plot(
  posterior_compare, "North Pacific loggerhead", "Annual population trend",
  col_cc, "Annual population change (%)", "North Pacific loggerhead — trend", zero_line = TRUE
)

p3b <- posterior_density_plot(
  posterior_compare, "North Pacific loggerhead", "Current abundance",
  col_cc, "Current total nester abundance", "North Pacific loggerhead — abundance"
)

p3c <- posterior_density_plot(
  posterior_compare, "West Pacific leatherback", "Annual population trend",
  col_dc, "Annual population change (%)", "West Pacific leatherback — trend", zero_line = TRUE
)

p3d <- posterior_density_plot(
  posterior_compare, "West Pacific leatherback", "Current abundance",
  col_dc, "Current total nester abundance", "West Pacific leatherback — abundance"
)


# ==============================================================================
# 26. COMBINE FIGURE 3
# ==============================================================================

figure3 <- ((p3a + p3b) / (p3c + p3d)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_memo_figure(
  figure3,
  "Figure_03_previous_vs_2026_posteriors",
  width = 8.4,
  height = 6.8
)


# ##############################################################################
# FIGURE 4: SSLL INTERACTIONS AND ANTICIPATED TAKE
# ##############################################################################

# ==============================================================================
# 27. READ ANNUAL SSLL INTERACTIONS
# ==============================================================================

fish_file <- find_output("ssll_annual_counts.csv")
fish      <- read.csv(fish_file, check.names = FALSE)
needed    <- c("year", "species", "n", "status")

if (!all(needed %in% names(fish))) {
  stop("Annual interaction file missing column(s): ", paste(setdiff(needed, names(fish)), collapse = ", "))
}

fish$species_std <- ifelse(
  grepl("logger|cc", fish$species, ignore.case = TRUE), "Loggerhead",
  ifelse(grepl("leather|dc", fish$species, ignore.case = TRUE), "Leatherback", NA_character_)
)
fish$complete <- toupper(fish$status) == "COMPLETE"


# ==============================================================================
# 28. INTERACTION PANEL FUNCTION
# ==============================================================================

make_interaction_panel <- function(dat, species, fill_colour) {
  x <- dat[dat$species_std == species, ]
  
  ggplot(x, aes(x = year, y = n)) +
    geom_col(data = x[x$complete, ], fill = fill_colour, width = 0.72) +
    geom_col(data = x[!x$complete, ], fill = "white", colour = fill_colour, linewidth = 0.7, width = 0.72) +
    scale_x_continuous(breaks = seq(2005, 2025, 4), minor_breaks = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
    labs(x = NULL, y = "Observed interactions") +
    theme_memo() +
    theme(legend.position = "none")
}

p4a <- make_interaction_panel(fish, "Loggerhead", col_cc)
p4c <- make_interaction_panel(fish, "Leatherback", col_dc)


# ==============================================================================
# 29. READ ATL DISTRIBUTIONS
# ==============================================================================

read_atl <- function(filename, method, species) {
  f <- find_output(filename)
  x <- read.csv(f)
  
  if (!all(c("ATL", "probability") %in% names(x))) {
    stop(filename, " must contain ATL and probability.")
  }
  
  x <- x[is.finite(x$ATL) & is.finite(x$probability), ]
  x$probability <- x$probability / sum(x$probability)
  x$Method      <- method
  x$Species     <- species
  x
}

atl <- rbind(
  read_atl("atl_loggerhead_Martin2020_CMP.csv", "2020 assessment", "Loggerhead"),
  read_atl("atl_loggerhead_update_CMP.csv",     "2026 update",     "Loggerhead"),
  read_atl("atl_leatherback_Martin2020_CMP.csv", "2020 assessment", "Leatherback"),
  read_atl("atl_leatherback_update_CMP.csv",     "2026 update",     "Leatherback")
)

atl$Method <- factor(atl$Method, levels = c("2020 assessment", "2026 update"))


# ==============================================================================
# 30. TRIM ATL PMFs
# ==============================================================================

trim_pmf <- function(x) {
  x  <- x[order(x$ATL), ]
  cp <- cumsum(x$probability)
  cutoff <- min(x$ATL[cp >= 0.999])
  x[x$ATL <= cutoff, ]
}

atl_plot <- do.call(
  rbind,
  lapply(split(atl, interaction(atl$Species, atl$Method, drop = TRUE)), trim_pmf)
)


# ==============================================================================
# 31. ATL PANEL FUNCTION (Comment 31 / 45: Axis labelled 'ATL')
# ==============================================================================

make_atl_panel <- function(dat, species, updated_colour) {
  x <- dat[dat$Species == species, ]
  
  ggplot(x, aes(x = ATL, y = probability, colour = Method, linetype = Method)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = c("2020 assessment" = col_mid, "2026 update" = updated_colour)) +
    scale_linetype_manual(values = c("2020 assessment" = "dashed", "2026 update" = "solid")) +
    labs(x = "ATL", y = "Probability") +
    theme_memo()
}

p4b <- make_atl_panel(atl_plot, "Loggerhead",  col_cc)
p4d <- make_atl_panel(atl_plot, "Leatherback", col_dc)


# ==============================================================================
# 32. COMBINE FIGURE 4
# ==============================================================================

figure4 <- ((p4a + p4b) / (p4c + p4d)) + plot_annotation(tag_levels = "A")

save_memo_figure(
  figure4,
  "Figure_04_SSLL_interactions_and_ATL",
  width = 8.2,
  height = 6.8
)


# ##############################################################################
# FIGURE 5: 10-YEAR FUTURE POPULATION PROJECTIONS (95% CrI)
# ##############################################################################

# ==============================================================================
# 33. LOAD PVA OBJECTS
# ==============================================================================

cc_pva_file <- find_output("update_PVA_loggerhead.rds")
dc_pva_file <- find_output("update_PVA_leatherback.rds")
notrip_file <- find_output("update_PVA_no_trip_limits_scenario.rds", required = FALSE)
nocap_file  <- find_output("update_PVA_no_hard_caps_scenario.rds",   required = FALSE)
both_file   <- find_output("update_PVA_both_removed_scenario.rds",   required = FALSE)

cc_pva_obj <- readRDS(cc_pva_file)
dc_pva_obj <- readRDS(dc_pva_file)

PVA_CC <- cc_pva_obj$PVA
PVA_DC <- dc_pva_obj$MEDIAN$PVA

notrip_obj <- NULL
nocap_obj  <- NULL
both_obj   <- NULL

has_notrip <- !is.na(notrip_file)
if (has_notrip) notrip_obj <- readRDS(notrip_file)

has_nocap <- !is.na(nocap_file)
if (has_nocap) nocap_obj <- readRDS(nocap_file)

has_both <- !is.na(both_file)
if (has_both) both_obj <- readRDS(both_file)


# ==============================================================================
# 34. SUMMARISE PVA TRAJECTORIES (RESTORED ENGINE)
# ==============================================================================

safe_log <- function(x) {
  out <- rep(NA_real_, length(x))
  pos <- which(is.finite(x) & x > 0)
  out[pos] <- log(x[pos])
  out
}

summarise_pva <- function(N, scenario, species) {
  N  <- as.matrix(N)
  qq <- apply(N, 1, quantile, probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  
  L95    <- qq[1, ]
  Median <- qq[2, ]
  U95    <- qq[3, ]
  
  data.frame(
    ProjectionYear = seq_len(nrow(N)),
    Species        = species,
    Scenario       = scenario,
    L95            = L95,
    Median         = Median,
    U95            = U95,
    LogL95         = safe_log(L95),
    LogMedian      = safe_log(Median),
    LogU95         = safe_log(U95)
  )
}


# ==============================================================================
# 35. COMBINE ALL PVA TRAJECTORIES ACROSS ACTIVE SCENARIOS
# ==============================================================================

build_scenario_plot_data <- function(pva_base, notrip, nocap, both, species) {
  res <- list(
    summarise_pva(pva_base$N_no_take, "No future SSLL take", species),
    summarise_pva(pva_base$N_take,    "Status quo SSLL take", species)
  )
  if (has_notrip && !is.null(notrip)) res[[length(res) + 1L]] <- summarise_pva(notrip$PVA$N_take, "No trip limits", species)
  if (has_nocap  && !is.null(nocap))  res[[length(res) + 1L]] <- summarise_pva(nocap$PVA$N_take,  "No hard caps",   species)
  if (has_both   && !is.null(both))   res[[length(res) + 1L]] <- summarise_pva(both$PVA$N_take,   "Both removed",   species)
  do.call(rbind, res)
}

cc_combined <- build_scenario_plot_data(PVA_CC, notrip_obj$CC, nocap_obj$CC, both_obj$CC, "North Pacific loggerhead")
dc_combined <- build_scenario_plot_data(PVA_DC, notrip_obj$DC, nocap_obj$DC, both_obj$DC, "West Pacific leatherback")

pva_plot_data <- rbind(cc_combined, dc_combined)

scenario_levels <- c("No future SSLL take", "Status quo SSLL take")
if (has_notrip) scenario_levels <- c(scenario_levels, "No trip limits")
if (has_nocap)  scenario_levels <- c(scenario_levels, "No hard caps")
if (has_both)   scenario_levels <- c(scenario_levels, "Both removed")

pva_plot_data$Scenario <- factor(pva_plot_data$Scenario, levels = scenario_levels)


# ==============================================================================
# 36. SPECIES-SPECIFIC PLOTTING FLOORS
# ==============================================================================

pva_plot_data$LogL95_plot <- pmax(pva_plot_data$LogL95, 0)
pva_plot_data$LogL95_plot[!is.finite(pva_plot_data$LogL95_plot)] <- 0


# ==============================================================================
# 37. SPECIES-SPECIFIC Y LIMITS
# ==============================================================================

get_pva_y_limits <- function(data, species) {
  x <- data[data$Species == species, ]
  upper_values <- c(x$LogU95, x$LogMedian)
  upper_values <- upper_values[is.finite(upper_values)]
  
  if (length(upper_values) == 0L) {
    stop("No finite upper PVA values found for ", species, ".")
  }
  c(0, ceiling(max(upper_values)))
}

cc_pva_ylim <- get_pva_y_limits(pva_plot_data, "North Pacific loggerhead")
dc_pva_ylim <- get_pva_y_limits(pva_plot_data, "West Pacific leatherback")


# ==============================================================================
# 38. PVA PANEL — MEDIAN + 95% CrI (10-YEAR PROJECTION)
# ==============================================================================

make_pva_ci <- function(data, species, panel_title, y_limits) {
  
  x <- data[data$Species == species, ]
  
  ggplot(x, aes(x = ProjectionYear)) +
    geom_ribbon(
      aes(ymin = LogL95_plot, ymax = LogU95, fill = Scenario, group = Scenario),
      alpha = 0.16,
      colour = NA,
      na.rm = TRUE
    ) +
    geom_line(
      aes(y = LogMedian, colour = Scenario, linetype = Scenario),
      linewidth = 0.95,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      values = c(
        "No future SSLL take"  = col_no_take,
        "Status quo SSLL take" = col_take,
        "No trip limits"       = col_no_trip,
        "No hard caps"         = col_no_cap,
        "Both removed"         = col_both
      )
    ) +
    scale_colour_manual(
      values = c(
        "No future SSLL take"  = col_no_take,
        "Status quo SSLL take" = col_take,
        "No trip limits"       = col_no_trip,
        "No hard caps"         = col_no_cap,
        "Both removed"         = col_both
      )
    ) +
    scale_linetype_manual(
      values = c(
        "No future SSLL take"  = "dashed",
        "Status quo SSLL take" = "solid",
        "No trip limits"       = "dotdash",
        "No hard caps"         = "dotted",
        "Both removed"         = "longdash"
      )
    ) +
    scale_x_continuous(breaks = seq(1, 10, 1), limits = c(1, 10), minor_breaks = NULL) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.01, 0.03))) +
    labs(
      title = panel_title,
      x = "Projection year",
      y = "ln(Annual Nesters)",
      fill = NULL,
      colour = NULL,
      linetype = NULL
    ) +
    theme_memo()
}


# ==============================================================================
# 39. BUILD TWO FIGURE 5 PANELS
# ==============================================================================

p5a <- make_pva_ci(
  pva_plot_data,
  species = "North Pacific loggerhead",
  panel_title = "North Pacific loggerhead",
  y_limits = cc_pva_ylim
)

p5b <- make_pva_ci(
  pva_plot_data,
  species = "West Pacific leatherback",
  panel_title = "West Pacific leatherback",
  y_limits = dc_pva_ylim
)


# ==============================================================================
# 40. COMBINE FIGURE 5 (2-PANEL VERTICAL)
# ==============================================================================

figure5 <- (p5a / p5b) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")


# ==============================================================================
# 41. SAVE FIGURE 5
# ==============================================================================

save_memo_figure(
  figure5,
  "Figure_05_future_PVA_take_vs_no_take",
  width = 7.5,
  height = 7.5
)


# ##############################################################################
# TABLE 1: PVA THRESHOLD PROBABILITIES
# ##############################################################################

# ==============================================================================
# 42. READ PVA THRESHOLD RESULTS
# ==============================================================================

pva_diff <- read.csv(find_output("update_PVA_threshold_differences.csv"))
needed   <- c(
  "species", "branch", "threshold", "horizon",
  "probability_below_take", "probability_below_no_take", "difference"
)

if (!all(needed %in% names(pva_diff))) {
  stop("PVA threshold output missing column(s): ", paste(setdiff(needed, names(pva_diff)), collapse = ", "))
}


# ==============================================================================
# 43. PRIMARY ANALYSES ONLY
# ==============================================================================

pva_table <- pva_diff[
  (pva_diff$species == "Loggerhead"  & pva_diff$branch == "PRIMARY") |
    (pva_diff$species == "Leatherback" & pva_diff$branch == "MEDIAN"),
]


# ==============================================================================
# 44. FORMAT TABLE
# ==============================================================================

pva_table$Species <- ifelse(
  pva_table$species == "Loggerhead",
  "North Pacific loggerhead",
  "West Pacific leatherback"
)

pva_table$Threshold <- ifelse(
  abs(pva_table$threshold - 0.50) < 1e-8, "50%",
  ifelse(abs(pva_table$threshold - 0.25) < 1e-8, "25%", "12.5%")
)

pva_table$Horizon_years                 <- pva_table$horizon
pva_table$No_future_SSLL_take_percent   <- round(100 * pva_table$probability_below_no_take, 1)
pva_table$Future_SSLL_take_percent      <- round(100 * pva_table$probability_below_take, 1)
pva_table$Difference_percentage_points  <- round(100 * pva_table$difference, 1)

pva_table <- pva_table[
  ,
  c(
    "Species",
    "Threshold",
    "Horizon_years",
    "No_future_SSLL_take_percent",
    "Future_SSLL_take_percent",
    "Difference_percentage_points"
  )
]

pva_table$Species   <- factor(pva_table$Species, levels = c("North Pacific loggerhead", "West Pacific leatherback"))
pva_table$Threshold <- factor(pva_table$Threshold, levels = c("50%", "25%", "12.5%"))

pva_table <- pva_table[order(pva_table$Species, pva_table$Horizon_years, pva_table$Threshold), ]
pva_table$Species   <- as.character(pva_table$Species)
pva_table$Threshold <- as.character(pva_table$Threshold)


# ==============================================================================
# 45. SAVE TABLE 1
# ==============================================================================

pva_table_file <- file.path(fig_dir, "Table_01_PVA_threshold_results.csv")
write.csv(pva_table, pva_table_file, row.names = FALSE)
cat("Saved: ", pva_table_file, "\n", sep = "")


# ==============================================================================
# 46. PRINT TABLE 1
# ==============================================================================

cat("\n============================================================\n")
cat("TABLE 1 — PVA THRESHOLD RESULTS (STATUS QUO)\n")
cat("============================================================\n")
print(pva_table, row.names = FALSE)
cat("============================================================\n\n")


# ==============================================================================
# 47. TABLE 1B — MULTI-SCENARIO RISK COMPARISON MATRIX
# ==============================================================================

scenario_table <- pva_table

add_scenario_col <- function(current_tab, csv_name, col_name) {
  f <- find_output(csv_name, required = FALSE)
  if (is.na(f)) return(current_tab)
  
  df      <- read.csv(f)
  df_take <- subset(df, scenario == "TAKE")
  
  df_take$Species <- ifelse(
    df_take$species == "Loggerhead",
    "North Pacific loggerhead",
    "West Pacific leatherback"
  )
  df_take$Threshold <- ifelse(
    abs(df_take$threshold - 0.50) < 1e-8, "50%",
    ifelse(abs(df_take$threshold - 0.25) < 1e-8, "25%", "12.5%")
  )
  df_take$Horizon_years <- df_take$horizon
  df_take[[col_name]]   <- round(100 * df_take$probability_below, 1)
  
  merge(
    current_tab,
    df_take[, c("Species", "Threshold", "Horizon_years", col_name)],
    by = c("Species", "Threshold", "Horizon_years"),
    all.x = TRUE
  )
}

scenario_table <- add_scenario_col(scenario_table, "update_PVA_threshold_probabilities_NO_TRIP_LIMITS.csv", "No_trip_limits_take_percent")
scenario_table <- add_scenario_col(scenario_table, "update_PVA_threshold_probabilities_NO_HARD_CAPS.csv",   "No_hard_caps_take_percent")
scenario_table <- add_scenario_col(scenario_table, "update_PVA_threshold_probabilities_BOTH_REMOVED.csv",   "Both_removed_take_percent")

if ("No_trip_limits_take_percent" %in% names(scenario_table)) {
  scenario_table$Delta_trip_limits <- round(scenario_table$No_trip_limits_take_percent - scenario_table$Future_SSLL_take_percent, 2)
}
if ("No_hard_caps_take_percent" %in% names(scenario_table)) {
  scenario_table$Delta_hard_caps <- round(scenario_table$No_hard_caps_take_percent - scenario_table$Future_SSLL_take_percent, 2)
}
if ("Both_removed_take_percent" %in% names(scenario_table)) {
  scenario_table$Delta_both_removed <- round(scenario_table$Both_removed_take_percent - scenario_table$Future_SSLL_take_percent, 2)
}

scenario_table <- scenario_table[
  order(scenario_table$Species, scenario_table$Horizon_years, scenario_table$Threshold),
]

scenario_table_file <- file.path(fig_dir, "Table_01_scenario_comparison.csv")
write.csv(scenario_table, scenario_table_file, row.names = FALSE)
cat("Saved: ", scenario_table_file, "\n", sep = "")

cat("\n============================================================\n")
cat("TABLE 1B — MULTI-SCENARIO RISK COMPARISON\n")
cat("============================================================\n")
print(scenario_table, row.names = FALSE)
cat("============================================================\n\n")


# ##############################################################################
# FIGURE 6: SCL AND POST-INTERACTION MORTALITY DISTRIBUTIONS
# ##############################################################################

# ==============================================================================
# 48. LOAD INTERACTION & MORTALITY DATA
# ==============================================================================

prep_file <- find_output("ssll_interactions_prepared.csv", required = FALSE)
raw_file  <- file.path(cfg$input_dir, cfg$files$interactions)

if (!is.na(prep_file) && file.exists(prep_file)) {
  df_raw <- read.csv(prep_file, check.names = FALSE, stringsAsFactors = FALSE)
  
  df_clean <- data.frame(
    Year    = as.integer(df_raw$year),
    Species = ifelse(
      grepl("logger|cc", df_raw$species, ignore.case = TRUE),
      "North Pacific loggerhead",
      "West Pacific leatherback"
    ),
    SCL     = as.numeric(df_raw$scl_cm),
    M_mean  = as.numeric(df_raw$m_mean)
  )
} else if (file.exists(raw_file)) {
  df_raw     <- read.csv(raw_file, check.names = FALSE, stringsAsFactors = FALSE)
  col_scl    <- grep("length.*scl.*cm", names(df_raw), ignore.case = TRUE, value = TRUE)[1]
  col_m_low  <- grep("mortality.*low",  names(df_raw), ignore.case = TRUE, value = TRUE)[1]
  col_m_high <- grep("mortality.*high", names(df_raw), ignore.case = TRUE, value = TRUE)[1]
  
  df_clean <- data.frame(
    Year    = as.integer(df_raw$Year),
    Species = ifelse(
      grepl("logger|cc", df_raw$Species, ignore.case = TRUE),
      "North Pacific loggerhead",
      "West Pacific leatherback"
    ),
    SCL     = as.numeric(df_raw[[col_scl]]),
    M_mean  = (as.numeric(df_raw[[col_m_low]]) + as.numeric(df_raw[[col_m_high]])) / 2
  )
} else {
  stop("Could not locate interaction dataset in outputs/input_audit/ or inputs/.")
}

df_clean <- df_clean[
  !is.na(df_clean$Species) &
    !is.na(df_clean$SCL)     &
    !is.na(df_clean$M_mean)  &
    df_clean$Year >= 2005    &
    df_clean$Year <= 2025,
]

df_18 <- df_clean[df_clean$Year <= 2018, ]
df_18$Timeframe <- "2005–2018"

df_25 <- df_clean[df_clean$Year <= 2025, ]
df_25$Timeframe <- "2005–2025"

df_dist <- rbind(df_18, df_25)
df_dist$Timeframe <- factor(df_dist$Timeframe, levels = c("2005–2018", "2005–2025"))


# ==============================================================================
# 49. DISTRIBUTION PANEL FUNCTION
# ==============================================================================

make_distribution_panel <- function(
    dat, species_name, variable_name, x_label,
    panel_title, update_color, x_breaks, x_limits = NULL
) {
  
  sub_d     <- dat[dat$Species == species_name, ]
  sub_d$Val <- sub_d[[variable_name]]
  means     <- aggregate(Val ~ Timeframe, data = sub_d, FUN = mean)
  
  ggplot(sub_d, aes(x = Val, fill = Timeframe, colour = Timeframe, linetype = Timeframe)) +
    # Scale density peak to 1.0 for each group:
    geom_density(
      aes(y = after_stat(density / max(density))),
      alpha = 0.22, linewidth = 0.9, adjust = 1.1, trim = TRUE
    ) +
    geom_vline(
      data = means,
      aes(xintercept = Val, colour = Timeframe, linetype = Timeframe),
      linewidth = 0.75,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c("2005–2018" = "#BDBDBD", "2005–2025" = update_color)) +
    scale_colour_manual(values = c("2005–2018" = col_previous, "2005–2025" = update_color)) +
    scale_linetype_manual(values = c("2005–2018" = "dashed", "2005–2025" = "solid")) +
    scale_x_continuous(breaks = x_breaks, limits = x_limits, expand = expansion(mult = c(0.01, 0.02))) +
    # Standardize y-axis across all 4 panels:
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.02))) +
    labs(title = panel_title, x = x_label, y = "Relative density") +
    theme_memo() +
    theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.35))
}


# ==============================================================================
# 50. ASSEMBLE & SAVE FIGURE 6 (Comment 34 / 48: Axis labelled 'SCL (cm)')
# ==============================================================================

p6a <- make_distribution_panel(
  dat = df_dist, species_name = "North Pacific loggerhead", variable_name = "SCL",
  x_label = "SCL (cm)", panel_title = "North Pacific loggerhead — length",
  update_color = col_cc, x_breaks = seq(40, 100, 10), x_limits = c(38, 102)
)

p6b <- make_distribution_panel(
  dat = df_dist, species_name = "North Pacific loggerhead", variable_name = "M_mean",
  x_label = "Mortality coefficient", panel_title = "North Pacific loggerhead — mortality",
  update_color = col_cc, x_breaks = seq(0, 1, 0.2), x_limits = c(0, 1.0)
)

p6c <- make_distribution_panel(
  dat = df_dist, species_name = "West Pacific leatherback", variable_name = "SCL",
  x_label = "SCL (cm)", panel_title = "West Pacific leatherback — length",
  update_color = col_dc, x_breaks = seq(100, 200, 20), x_limits = c(95, 205)
)

p6d <- make_distribution_panel(
  dat = df_dist, species_name = "West Pacific leatherback", variable_name = "M_mean",
  x_label = "Mortality coefficient", panel_title = "West Pacific leatherback — mortality",
  update_color = col_dc, x_breaks = seq(0, 1, 0.2), x_limits = c(0, 1.0)
)

figure6 <- ((p6a + p6b) / (p6c + p6d)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_memo_figure(
  figure6,
  "Figure_06_SCL_and_mortality_distributions",
  width = 8.4,
  height = 6.8
)


# ##############################################################################
# FIGURE 7: COUNTERFACTUAL TAKE AND ANE DISTRIBUTIONS
# ##############################################################################

# ==============================================================================
# 51. ASSEMBLE & SAVE FIGURE 7
# ==============================================================================

if (has_notrip || has_nocap || has_both) {
  
  cat("\nGenerating Figure 7 (Counterfactual take & ANE comparison)...\n")
  
  make_take_df <- function(obj_base, obj_notrip, obj_nocap, obj_both, species) {
    df_list <- list(
      data.frame(
        Species  = species,
        Scenario = "Status quo SSLL take",
        Takes    = as.vector(obj_base$future_take$ATL),
        ANE      = as.vector(obj_base$future_take$ANE)
      )
    )
    if (has_notrip && !is.null(obj_notrip)) {
      df_list[[length(df_list) + 1L]] <- data.frame(
        Species  = species,
        Scenario = "No trip limits",
        Takes    = as.vector(obj_notrip$future_take$ATL),
        ANE      = as.vector(obj_notrip$future_take$ANE)
      )
    }
    if (has_nocap && !is.null(obj_nocap)) {
      df_list[[length(df_list) + 1L]] <- data.frame(
        Species  = species,
        Scenario = "No hard caps",
        Takes    = as.vector(obj_nocap$future_take$ATL),
        ANE      = as.vector(obj_nocap$future_take$ANE)
      )
    }
    if (has_both && !is.null(obj_both)) {
      df_list[[length(df_list) + 1L]] <- data.frame(
        Species  = species,
        Scenario = "Both removed",
        Takes    = as.vector(obj_both$future_take$ATL),
        ANE      = as.vector(obj_both$future_take$ANE)
      )
    }
    do.call(rbind, df_list)
  }
  
  df_take_cc <- make_take_df(cc_pva_obj, notrip_obj$CC, nocap_obj$CC, both_obj$CC, "North Pacific loggerhead")
  df_take_dc <- make_take_df(dc_pva_obj$MEDIAN, notrip_obj$DC, nocap_obj$DC, both_obj$DC, "West Pacific leatherback")
  
  df_takes_all <- rbind(df_take_cc, df_take_dc)
  df_takes_all$Scenario <- factor(df_takes_all$Scenario, levels = scenario_levels[-1]) # omit No Take
  
  plot_take_comp <- function(dat, spp, var, x_lab, title_txt, max_x) {
    sub_d <- dat[dat$Species == spp & dat[[var]] <= max_x, ]
    
    ggplot(sub_d, aes(x = .data[[var]], fill = Scenario, colour = Scenario, linetype = Scenario)) +
      geom_density(alpha = 0.18, linewidth = 0.85, adjust = 1.2, trim = TRUE) +
      scale_fill_manual(
        values = c(
          "Status quo SSLL take" = col_take,
          "No trip limits"       = col_no_trip,
          "No hard caps"         = col_no_cap,
          "Both removed"         = col_both
        )
      ) +
      scale_colour_manual(
        values = c(
          "Status quo SSLL take" = col_take,
          "No trip limits"       = col_no_trip,
          "No hard caps"         = col_no_cap,
          "Both removed"         = col_both
        )
      ) +
      scale_linetype_manual(
        values = c(
          "Status quo SSLL take" = "solid",
          "No trip limits"       = "dotdash",
          "No hard caps"         = "dotted",
          "Both removed"         = "longdash"
        )
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      labs(
        title    = title_txt,
        x        = x_lab,
        y        = "Density",
        fill     = NULL,
        colour   = NULL,
        linetype = NULL
      ) +
      theme_memo() +
      theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.35))
  }
  
  p7a <- plot_take_comp(
    df_takes_all, "North Pacific loggerhead", "Takes",
    "Simulated annual interactions", "North Pacific loggerhead — annual takes", 55
  )
  p7b <- plot_take_comp(
    df_takes_all, "North Pacific loggerhead", "ANE",
    "Adult Nester Equivalents (ANE)", "North Pacific loggerhead — annual ANE", 8
  )
  p7c <- plot_take_comp(
    df_takes_all, "West Pacific leatherback", "Takes",
    "Simulated annual interactions", "West Pacific leatherback — annual takes", 30
  )
  p7d <- plot_take_comp(
    df_takes_all, "West Pacific leatherback", "ANE",
    "Adult Nester Equivalents (ANE)", "West Pacific leatherback — annual ANE", 4
  )
  
  figure7 <- ((p7a + p7b) / (p7c + p7d)) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")
  
  save_memo_figure(
    figure7,
    "Figure_07_take_and_ANE_no_trip_limits_comparison",
    width = 8.6,
    height = 7.0
  )
}


# ==============================================================================
# 52. FINISHED
# ==============================================================================

cat("\n============================================================\n")
cat("TECH MEMO FIGURES COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Created:\n",
  "  Figure 1  — Nesting data used in 2026 update (Panels A, B)\n",
  "  Figure 2  — Updated population trend model fits (Panels A, B)\n",
  "  Figure 3  — 2020 assessment vs 2026 posteriors (Panels A-D)\n",
  "  Figure 4  — SSLL interactions and ATL (Panels A-D)\n",
  "  Figure 5  — 10-year future PVA projections with 95% CrI (Panels A, B)\n",
  "  Table 1   — PVA threshold probabilities (Status Quo)\n",
  "  Table 1B  — Multi-scenario risk comparison matrix\n",
  "  Figure 6  — SCL and mortality distributions (Panels A-D)\n",
  "  Figure 7  — Counterfactual take and ANE distributions (Panels A-D)\n\n",
  sep = ""
)

cat("PNG only. No PDFs created.\n")
cat("No models were fitted or rerun.\n")
cat("No accepted analysis outputs were modified.\n\n")
cat("Files saved in:\n  ", fig_dir, "\n", sep = "")
cat("============================================================\n")