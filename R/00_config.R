# ==============================================================================
# 00_config.R
#
# Central configuration for the 2026 SSLL update.
#
# IMPORTANT:
#   - Source CSVs are never modified.
#   - "baseline" means Martin-method analysis using the current PIRO dataset,
#     restricted to the endpoints available to the previous assessment.
#   - It is NOT intended to recreate Martin's original fishery input dataset.
#   - "update" extends the same validated methods to the new data.
# ==============================================================================

cfg <- list(
  
  # --------------------------------------------------------------------------
  # FILES
  # --------------------------------------------------------------------------
  
  input_dir = "inputs",
  
  files = list(
    cc_nesting = "Yakushima_Cc_1985-2025.csv",
    dc_nesting = "WestPapua_Dc_1980-2025.csv",
    interactions = "ssll_CcDc_2005-2025.csv"
  ),
  
  output_dirs = list(
    audit      = "outputs/input_audit",
    nesting    = "outputs/nesting",
    validation = "outputs/validation_martin_baseline",
    update     = "outputs/update_2026"
  ),
  
  # --------------------------------------------------------------------------
  # BASELINE ANALYSIS WINDOWS
  #
  # These reproduce the ENDPOINTS of the previous analysis using the PIRO
  # interaction dataset currently supplied to us.
  #
  # The PIRO interaction file begins in 2005, so 2004 is not manufactured,
  # imputed, or assumed to be zero.
  # --------------------------------------------------------------------------
  
  baseline = list(
    cc_nesting_start = 1986L,
    cc_nesting_end   = 2015L,
    
    dc_nesting_start = 2001L,
    dc_nesting_end   = 2017L,
    
    interaction_start = 2005L,
    interaction_end   = 2018L
  ),
  
  # --------------------------------------------------------------------------
  # UPDATED ANALYSIS WINDOWS
  # --------------------------------------------------------------------------
  
  update = list(
    cc_nesting_end = 2025L,
    
    # Leatherback endpoint will be determined from the last complete
    # April-March season supported by the source file.
    dc_nesting_end = NA_integer_,
    
    interaction_start = 2005L,
    interaction_end   = 2025L
  ),
  
  # --------------------------------------------------------------------------
  # KNOWN DATA STATUS
  # --------------------------------------------------------------------------
  
  data_status = list(
    
    # Loggerhead SSLL interactions for these years have not yet been supplied.
    # These MUST NOT be interpreted as zero interactions in the updated
    # analysis.
    cc_interaction_years_pending = c(2022L, 2023L),
    
    # No currently known pending leatherback years.
    dc_interaction_years_pending = integer(0)
  ),
  
  # --------------------------------------------------------------------------
  # MARTIN 2020 BIOLOGICAL PARAMETERS
  #
  # Values below follow the implementation used in the Martin Appendix code.
  # --------------------------------------------------------------------------
  
  biology = list(
    
    cc = list(
      species = "Loggerhead",
      
      # Nesting
      CF = 3.0,
      RI = 3.3,
      
      # Survival / sex
      Pj = 0.80,
      Pa = 0.895,
      PF = 0.65,
      
      # von Bertalanffy growth model
      Linf  = 80.4473850,
      K     = 0.1396317,
      Lknot = 4.7363329,
      Amat  = 26.4950786
    ),
    
    dc = list(
      species = "Leatherback",
      
      # Nesting
      CF = 5.5,
      RI = 3.06,
      
      # Survival / sex
      Pj = 0.81,
      Pa = 0.893,
      PF = 0.73,
      
      # von Bertalanffy growth model
      Linf  = 142.7,
      K     = 0.2262,
      tknot = -0.17,
      Amat  = 16.1
    )
  ),
  
  # --------------------------------------------------------------------------
  # MARTIN TREND MODEL SETTINGS
  #
  # These are the settings from the implementation that successfully
  # reproduced the Martin 2020 baseline.
  # --------------------------------------------------------------------------
  
  trend = list(
    a_mean = 0,
    a_sd   = 4,
    
    u_mean = 0,
    u_sd   = 0.5,
    
    q_alpha = 0.01,
    q_beta  = 0.01,
    
    r_alpha = 0.01,
    r_beta  = 0.01,
    
    x0_sd = 10,
    
    chains  = 2L,
    samples = 10000L,
    burnin  = 5000L,
    thin    = 50L,
    
    seed = 132L
  ),
  
  # --------------------------------------------------------------------------
  # MARTIN LEATHERBACK IMPUTATION SETTINGS
  # --------------------------------------------------------------------------
  
  imputation = list(
    jm_period = 12L,
    w_period  = 6L,
    
    chains  = 5L,
    samples = 100000L,
    burnin  = 50000L,
    thin    = 5L,
    
    seed = 132L
  ),
  
  # --------------------------------------------------------------------------
  # PVA SETTINGS
  # --------------------------------------------------------------------------
  
  pva = list(
    projection_years = 100L,
    
    abundance_thresholds = c(
      0.50,
      0.25,
      0.125
    ),
    
    reporting_years = c(
      5L,
      10L,
      25L,
      50L,
      100L
    )
  )
)

# ==============================================================================
# BASIC CONFIGURATION CHECKS
# ==============================================================================

stopifnot(
  cfg$baseline$interaction_start == 2005L,
  cfg$baseline$interaction_end   == 2018L,
  
  cfg$baseline$cc_nesting_end == 2015L,
  cfg$baseline$dc_nesting_end == 2017L,
  
  cfg$update$interaction_end == 2025L,
  
  cfg$biology$cc$CF == 3.0,
  cfg$biology$dc$CF == 5.5,
  
  cfg$trend$u_sd == 0.5,
  cfg$trend$q_alpha == 0.01,
  cfg$trend$r_alpha == 0.01
)

message("00_config.R loaded successfully.")
message(
  "Baseline: CC nesting through ", cfg$baseline$cc_nesting_end,
  "; DC nesting through ", cfg$baseline$dc_nesting_end,
  "; PIRO SSLL interactions ",
  cfg$baseline$interaction_start, "-", cfg$baseline$interaction_end, "."
)
message(
  "Update: nesting through current endpoints; PIRO SSLL interactions ",
  cfg$update$interaction_start, "-", cfg$update$interaction_end, "."
)
message(
  "Pending loggerhead interaction years: ",
  paste(cfg$data_status$cc_interaction_years_pending, collapse = ", "),
  " — these will NOT be treated as zero."
)