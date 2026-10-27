#!/usr/bin/env Rscript

# ==============================================================================
# COLNE P-FASE 2 km INPUT PREPARATION
# UKCEH LCM2024 25 m + CEH-GEAR DAILY 2019 NETCDF
# ==============================================================================
#
# ----------------------
# CEH_GEAR_daily_GB_2019.nc is ONE annual NetCDF file containing 365 daily
# rainfall grids. The current NanoFASE `nfdata` workflow is most reliable with
# one raster per timestep for spatiotemporal raster variables. This script reads
# the annual CEH-GEAR NetCDF one day at a time and writes:
#
#   precip/precip_0.tif
#   ...
#   precip/precip_364.tif
#
# Each CEH-GEAR 1 km daily rainfall field is AREA-AVERAGED to the exact Colne
# 2 km model grid, then converted from mm/day to m/day.
#
# UKCEH LCM2024 is handled differently:
# - band 1 = land-cover class identifier (1..21)
# - bands 2/3 = probability/purity QA layers
# - the class raster is clipped to Colne but KEPT at 25 m.
#   nfdata then maps those classes and computes fractional NanoFASE land use
#   within each 2 km model cell.
#
# JULES hydrology
# ----------------
# The completed JULES production run is read directly from the twelve monthly
# pfase_hydrology NetCDF files. Daily JULES fluxes are spatially averaged from
# the 1 km GB grid to the same exact Colne 2 km BNG grid used by all other
# model inputs and converted from kg m-2 s-1 to m/day.
#
# Mapping into the current P-FASE/NanoFASE input names:
#   JULES runoff         -> runoff
#   JULES surf_roff      -> quickflow (overland runoff used by erosion)
#   JULES fqw_gb         -> evap (negative condensation/dew clipped to zero)
#
# JULES sub_surf_roff and drain are also written as separate 2 km daily
# products for QA and later PFAS leaching work, but the current nfdata/P-FASE
# input interface does not expose them as separate time-series variables.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(ncdf4)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------------------------

ROOT <- "/users/pollution/canulu/Desktop/nanofase/data.example/PFAS_data"

SCENARIO_DIR <- file.path(ROOT, "Colne_pfase_2km_2019")
MODEL_DIR <- file.path(SCENARIO_DIR, "processed_2km")
QA_DIR <- file.path(SCENARIO_DIR, "qa")
CALIB_DIR <- file.path(SCENARIO_DIR, "candidate_calibration_2km")
PRECIP_DIR <- file.path(MODEL_DIR, "precip")
RUNOFF_DIR_OUT <- file.path(MODEL_DIR, "runoff")
QUICKFLOW_DIR_OUT <- file.path(MODEL_DIR, "quickflow")
EVAP_DIR_OUT <- file.path(MODEL_DIR, "evap")
JULES_SUBSURF_DIR_OUT <- file.path(MODEL_DIR, "jules_subsurface_runoff")
JULES_DRAIN_DIR_OUT <- file.path(MODEL_DIR, "jules_drainage")
OUTPUT_DIR <- file.path(SCENARIO_DIR, "output")
MODEL_OUTPUT_DIR <- file.path(SCENARIO_DIR, "model_output")

for (d in c(SCENARIO_DIR, MODEL_DIR, QA_DIR, CALIB_DIR, PRECIP_DIR,
            RUNOFF_DIR_OUT, QUICKFLOW_DIR_OUT, EVAP_DIR_OUT,
            JULES_SUBSURF_DIR_OUT, JULES_DRAIN_DIR_OUT,
            OUTPUT_DIR, MODEL_OUTPUT_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

TARGET_CRS <- "EPSG:27700"
TARGET_RES <- 2000

SIM_YEAR <- 2019L
N_TIMESTEPS <- 365L
START_DATE <- "2019-01-01"

# Catchment polygon
CATCHMENT_SHP <- file.path(
  ROOT, "Colne", "WFD_Surface_Water_Management_Catchments_Cycle_3.shp"
)
CATCHMENT_PATTERN <- "Colne"
ALLOW_UNION_ALL_IF_NO_NAME_MATCH <- FALSE

# Original HydroSHEDS 15 arc-second DEM used by the earlier NanoFASE
# preprocessing notebook. This is the correct source for the new 2 km Colne
# grid; do NOT resample the old 10 km flowdir.tif or dem_conditioned.tif.
#
# The code checks both likely locations so you can keep the file either directly
# under Eu_data or under Eu_data/input_data.
HYDROSHEDS_DEM_CANDIDATES <- c(
  file.path(ROOT, "hydrosheds_eu_dem_15s.tif"),
  file.path(ROOT, "input_data", "hydrosheds_eu_dem_15s.tif")
)

HYDROSHEDS_DEM <- HYDROSHEDS_DEM_CANDIDATES[
  file.exists(HYDROSHEDS_DEM_CANDIDATES)
][1]

# Buffer around the Colne polygon before cropping the native HydroSHEDS DEM.
DEM_BUFFER_M <- 10000

# Existing 10 km DEM is retained only for QA/comparison and is never used to
# create the scientific 2 km routing grid.
OLD_DEM_10KM <- file.path(ROOT, "dem_conditioned.tif")

# New rainfall source
GEAR_NC <- file.path(ROOT, "CEH_GEAR_daily_GB_2019.nc")
GEAR_MISSING_NC <- file.path(ROOT, "CEH_GEAR_daily_GB_2019_missing_records.nc")
GEAR_RAIN_VAR_CANDIDATES <- c(
  "rainfall_amount", "rainfall", "precipitation", "precip"
)
GEAR_MISSING_VAR_CANDIDATES <- c(
  "missing_records", "missing_estimates", "n_missing"
)

# New land-cover source
LCM2024 <- file.path(ROOT, "gblcm2024_25m.tif")

# JULES 2019 hydrology production output.
# This run was warm-started and successfully completed 365 daily records.
JULES_OUTPUT_DIR <- "/users/pollution/canulu/Desktop/nanofase/data.example/PFAS_data/jules_hydrology_2019"
JULES_FILE_STEM <- "gb2019_hydrology_2019.pfase_hydrology"
JULES_MONTHLY_FILES <- file.path(
  JULES_OUTPUT_DIR,
  sprintf("%s.%04d%02d.nc", JULES_FILE_STEM, SIM_YEAR, 1:12)
)

# JULES output variables used here.
JULES_RUNOFF_VAR <- "runoff"
JULES_QUICKFLOW_VAR <- "surf_roff"
JULES_SUBSURF_VAR <- "sub_surf_roff"
JULES_DRAIN_VAR <- "drain"
JULES_EVAP_VAR <- "fqw_gb"

# JULES rates are kg m-2 s-1. For liquid water, 1 kg m-2 = 1 mm, so:
# kg m-2 s-1 * 86400 / 1000 = m/day.
JULES_FLUX_TO_M_DAY <- 86400 / 1000

# Expected native JULES/CHESS grid spacing, used only as a geometry sanity check.
JULES_EXPECTED_RES_M <- 1000
JULES_GRID_RES_TOL_M <- 25

# First run options
USE_ESTUARY <- FALSE

# Supplied P-factor raster was found to be a copy of soil-loss data.
P_FACTOR_VALUE <- 1.0

# Existing cfvo represents VOLUME fraction while NanoFASE documentation labels
# coarse fragments kg/kg. Preserve real cfvo in QA and use zero in model until
# that semantic mismatch is resolved.
COARSE_FRAG_MODEL_VALUE <- 0.0

# PFAS software-test areal loading to water.
PFAS_TEST_LOAD_KG_M2_DAY <- 1.0e-12

# ------------------------------------------------------------------------------
# STATIC SOURCE DATA
# ------------------------------------------------------------------------------

static_src <- list(
  field_capacity = file.path(ROOT, "field_capacity_processed.tif"),
  saturation = file.path(ROOT, "saturation_processed.tif"),
  hydraulic_conductivity = file.path(ROOT, "hydraulic_conductivity_processed.tif"),
  
  c_factor = file.path(ROOT, "filled_cfactor_processed.tif"),
  ls_factor = file.path(ROOT, "filled_lsfactor_processed.tif"),
  
  # Candidate calibration layers (processed for QA, not enabled in first YAML)
  bank_erosion_alpha = file.path(ROOT, "bank_erosion_alpha.tif"),
  bank_erosion_beta = file.path(ROOT, "bank_erosion_beta.tif"),
  deposition_alpha = file.path(ROOT, "deposition_alpha.tif"),
  deposition_beta = file.path(ROOT, "deposition_beta_spatial.tif"),
  resuspension_alpha = file.path(ROOT, "resuspension_alpha.tif"),
  resuspension_beta = file.path(ROOT, "resuspension_beta.tif"),
  sediment_transport_a = file.path(ROOT, "sediment_transport_a.tif"),
  sediment_transport_b = file.path(ROOT, "sediment_transport_b.tif"),
  sediment_transport_c = file.path(ROOT, "sediment_transport_c.tif"),
  
  bad_p_factor = file.path(ROOT, "filled_pfactor_processed.tif"),
  soil_loss = file.path(ROOT, "soilloss_processed.tif")
)

SOIL_RAW_DIR <- file.path(ROOT, "input_data")

soil_depth_files <- list(
  sand = c(
    "sand_0-5cm_mean_1000.tif",
    "sand_5-15cm_mean_1000.tif",
    "sand_15-30cm_mean_1000.tif"
  ),
  clay = c(
    "clay_0-5cm_mean_1000.tif",
    "clay_5-15cm_mean_1000.tif",
    "clay_15-30cm_mean_1000.tif"
  ),
  silt = c(
    "silt_0-5cm_mean_1000.tif",
    "silt_5-15cm_mean_1000.tif",
    "silt_15-30cm_mean_1000.tif"
  ),
  cfvo = c(
    "cfvo_0-5cm_mean_1000.tif",
    "cfvo_5-15cm_mean_1000.tif",
    "cfvo_15-30cm_mean_1000.tif"
  ),
  bdod = c(
    "bdod_0-5cm_mean_1000.tif",
    "bdod_5-15cm_mean_1000.tif",
    "bdod_15-30cm_mean_1000.tif"
  )
)

depth_weights <- c(5, 10, 15)

soil_fallback <- list(
  sand = file.path(ROOT, "sand_0-30cm_mean_1000.tif"),
  clay = file.path(ROOT, "clay_0-30cm_mean_1000.tif"),
  silt = file.path(ROOT, "silt_0-30cm_mean_1000.tif"),
  cfvo = file.path(ROOT, "cfvo_0-30cm_mean_1000.tif"),
  bdod = file.path(ROOT, "bdod_0-30cm_mean_1000.tif")
)

# ------------------------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------------------------

stop_missing <- function(paths, label = "required file") {
  paths <- as.character(paths)
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      paste0(label, "(s) missing:\n", paste(" -", missing, collapse = "\n")),
      call. = FALSE
    )
  }
}

stats4 <- function(r) {
  c(
    min = global(r, "min", na.rm = TRUE)[1, 1],
    median = global(r, "median", na.rm = TRUE)[1, 1],
    mean = global(r, "mean", na.rm = TRUE)[1, 1],
    max = global(r, "max", na.rm = TRUE)[1, 1]
  )
}

show_stats <- function(label, r) {
  s <- stats4(r)
  cat(sprintf(
    "%-48s min=%g median=%g mean=%g max=%g\n",
    label, s["min"], s["median"], s["mean"], s["max"]
  ))
}

write_float <- function(r, path) {
  writeRaster(
    r, path,
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2")
  )
}

write_int <- function(r, path, datatype = "INT4S") {
  naflag <- if (toupper(datatype) == "INT1U") 255 else -9999
  writeRaster(
    r, path,
    overwrite = TRUE,
    datatype = datatype,
    NAflag = naflag,
    gdal = c("COMPRESS=DEFLATE")
  )
}

make_2km_template <- function(catchment_sf) {
  b <- st_bbox(catchment_sf)
  
  xmin <- floor(as.numeric(b["xmin"]) / TARGET_RES) * TARGET_RES
  xmax <- ceiling(as.numeric(b["xmax"]) / TARGET_RES) * TARGET_RES
  ymin <- floor(as.numeric(b["ymin"]) / TARGET_RES) * TARGET_RES
  ymax <- ceiling(as.numeric(b["ymax"]) / TARGET_RES) * TARGET_RES
  
  rast(
    xmin = xmin, xmax = xmax,
    ymin = ymin, ymax = ymax,
    resolution = TARGET_RES,
    crs = TARGET_CRS
  )
}

to_2km <- function(path, template, catch_mask, method = "bilinear") {
  if (!file.exists(path)) stop("Missing raster: ", path, call. = FALSE)
  r <- rast(path)
  if (nlyr(r) > 1) r <- r[[1]]
  
  if (is.na(crs(r)) || crs(r) == "") {
    stop("Raster has no CRS: ", path, call. = FALSE)
  }
  
  z <- project(r, template, method = method, threads = TRUE)
  mask(z, catch_mask)
}

depth_weighted_mean <- function(paths, weights) {
  stop_missing(paths, "soil depth layer")
  rr <- lapply(paths, rast)
  
  weighted <- NULL
  weight_total <- NULL
  
  for (i in seq_along(rr)) {
    valid <- ifel(!is.na(rr[[i]]), 1, 0)
    contribution <- rr[[i]] * weights[i] * valid
    w <- weights[i] * valid
    
    if (is.null(weighted)) {
      weighted <- contribution
      weight_total <- w
    } else {
      weighted <- weighted + contribution
      weight_total <- weight_total + w
    }
  }
  
  ifel(weight_total > 0, weighted / weight_total, NA)
}

convert_water_depth_to_m_day <- function(r, units) {
  u <- tolower(gsub(" ", "", units))
  
  if (u %in% c("m/day", "m/d", "mday-1")) return(r)
  if (u %in% c("mm/day", "mm/d", "mmday-1", "mm")) return(r / 1000)
  if (u %in% c("cm/day", "cm/d", "cmday-1")) return(r / 100)
  
  stop("Unsupported hydrological units: ", units, call. = FALSE)
}

# Hydraulic-conductivity raster audit showed values are encoded as
# (cm/day) * 100.
ks_encoded_to_m_s <- function(r) {
  (r / 100) * 0.01 / 86400
}

clear_daily_outputs <- function(out_dir, prefix) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  old <- list.files(
    out_dir,
    pattern = paste0("^", prefix, "_[0-9]+\\.tif$"),
    full.names = TRUE
  )
  if (length(old) > 0) unlink(old)
}

normalise_units <- function(x) {
  tolower(gsub("[[:space:]_]+", "", x))
}

check_jules_flux_units <- function(nc, var_name) {
  att <- ncatt_get(nc, var_name, "units")
  if (!isTRUE(att$hasatt)) {
    stop("JULES variable '", var_name, "' has no units attribute.", call. = FALSE)
  }
  
  u <- normalise_units(att$value)
  acceptable <- c(
    "kgm-2s-1", "kgm^-2s^-1", "kg/m2/s", "kgm**-2s**-1"
  )
  
  if (!(u %in% acceptable)) {
    stop(
      "Expected JULES ", var_name, " in kg m-2 s-1, found: ", att$value,
      call. = FALSE
    )
  }
  
  invisible(att$value)
}

nc_var_dim_names <- function(nc, var_name) {
  vapply(nc$var[[var_name]]$dim, function(d) d$name, character(1))
}

# Return any 2-D NetCDF variable as a matrix ordered [y, x].
nc_read_static_yx <- function(nc, var_name) {
  v <- nc$var[[var_name]]
  dn <- nc_var_dim_names(nc, var_name)
  z <- ncvar_get(nc, var_name, collapse_degen = TRUE)
  
  if (length(dim(z)) != 2) {
    stop("Expected 2-D variable '", var_name, "'.", call. = FALSE)
  }
  
  iy <- match("y", dn)
  ix <- match("x", dn)
  
  if (any(is.na(c(iy, ix)))) {
    stop(
      "Variable '", var_name, "' does not use x/y dimensions: ",
      paste(dn, collapse = ", "), call. = FALSE
    )
  }
  
  if (identical(dn, c("y", "x"))) return(z)
  if (identical(dn, c("x", "y"))) return(t(z))
  
  stop(
    "Unexpected dimension order for '", var_name, "': ",
    paste(dn, collapse = ", "), call. = FALSE
  )
}

# Read one daily JULES field and return a matrix ordered [y, x].
nc_read_daily_yx <- function(nc, var_name, time_index) {
  v <- nc$var[[var_name]]
  dn <- nc_var_dim_names(nc, var_name)
  
  it <- match("time", dn)
  iy <- match("y", dn)
  ix <- match("x", dn)
  
  if (any(is.na(c(it, iy, ix)))) {
    stop(
      "JULES variable '", var_name, "' dimensions are ",
      paste(dn, collapse = ", "), "; expected time,y,x.",
      call. = FALSE
    )
  }
  
  start <- rep(1L, length(dn))
  count <- vapply(v$dim, function(d) d$len, integer(1))
  start[it] <- time_index
  count[it] <- 1L
  
  z <- ncvar_get(
    nc,
    var_name,
    start = start,
    count = count,
    collapse_degen = TRUE
  )
  
  rem <- dn[-it]
  
  if (length(dim(z)) != 2) {
    stop(
      "Unexpected daily slice shape for JULES variable '", var_name, "'.",
      call. = FALSE
    )
  }
  
  if (identical(rem, c("y", "x"))) return(z)
  if (identical(rem, c("x", "y"))) return(t(z))
  
  stop(
    "Unexpected remaining dimension order for '", var_name, "': ",
    paste(rem, collapse = ", "), call. = FALSE
  )
}

# Reconstruct the regular 1 km British National Grid geometry from the 2-D
# latitude/longitude fields written by JULES. The JULES output file uses x/y
# dimensions but does not need to expose BNG x/y coordinate variables for this
# conversion to work.
infer_jules_bng_geometry <- function(nc) {
  if (!("latitude" %in% names(nc$var)) || !("longitude" %in% names(nc$var))) {
    stop("JULES output lacks latitude/longitude coordinate variables.", call. = FALSE)
  }
  
  lat <- nc_read_static_yx(nc, "latitude")
  lon <- nc_read_static_yx(nc, "longitude")
  
  if (!identical(dim(lat), dim(lon))) {
    stop("JULES latitude/longitude shapes differ.", call. = FALSE)
  }
  
  ny <- nrow(lat)
  nx <- ncol(lat)
  mid_y <- max(1L, as.integer(round(ny / 2)))
  mid_x <- max(1L, as.integer(round(nx / 2)))
  
  row_df <- data.frame(lon = lon[mid_y, ], lat = lat[mid_y, ])
  col_df <- data.frame(lon = lon[, mid_x], lat = lat[, mid_x])
  
  row_sf <- st_as_sf(row_df, coords = c("lon", "lat"), crs = "EPSG:4326")
  col_sf <- st_as_sf(col_df, coords = c("lon", "lat"), crs = "EPSG:4326")
  
  row_xy <- st_coordinates(st_transform(row_sf, TARGET_CRS))
  col_xy <- st_coordinates(st_transform(col_sf, TARGET_CRS))
  
  x_centres <- row_xy[, 1]
  y_centres <- col_xy[, 2]
  
  dx <- median(abs(diff(x_centres)), na.rm = TRUE)
  dy <- median(abs(diff(y_centres)), na.rm = TRUE)
  
  if (!is.finite(dx) || !is.finite(dy)) {
    stop("Could not infer JULES grid spacing from latitude/longitude.", call. = FALSE)
  }
  
  if (abs(dx - JULES_EXPECTED_RES_M) > JULES_GRID_RES_TOL_M ||
      abs(dy - JULES_EXPECTED_RES_M) > JULES_GRID_RES_TOL_M) {
    stop(
      "Unexpected JULES grid spacing inferred from lat/lon: ",
      sprintf("dx=%.3f m, dy=%.3f m", dx, dy),
      call. = FALSE
    )
  }
  
  src <- rast(
    nrows = ny,
    ncols = nx,
    xmin = min(x_centres, na.rm = TRUE) - dx / 2,
    xmax = max(x_centres, na.rm = TRUE) + dx / 2,
    ymin = min(y_centres, na.rm = TRUE) - dy / 2,
    ymax = max(y_centres, na.rm = TRUE) + dy / 2,
    crs = TARGET_CRS
  )
  
  list(
    template = src,
    flip_rows = isTRUE(y_centres[1] < y_centres[length(y_centres)]),
    flip_cols = isTRUE(x_centres[1] > x_centres[length(x_centres)]),
    dx = dx,
    dy = dy,
    nx = nx,
    ny = ny
  )
}

jules_matrix_to_raster <- function(z_yx, geom) {
  z <- z_yx
  
  if (geom$flip_rows) {
    z <- z[nrow(z):1, , drop = FALSE]
  }
  
  if (geom$flip_cols) {
    z <- z[, ncol(z):1, drop = FALSE]
  }
  
  if (nrow(z) != geom$ny || ncol(z) != geom$nx) {
    stop("JULES data slice does not match inferred grid dimensions.", call. = FALSE)
  }
  
  rast(z, extent = ext(geom$template), crs = TARGET_CRS)
}

jules_to_colne_2km <- function(r, model_template, catch_mask) {
  # These variables are water-depth fluxes. Spatial aggregation from 1 km to
  # 2 km therefore uses an area mean, not a sum.
  z <- project(r, model_template, method = "average", threads = TRUE)
  mask(z, catch_mask)
}

process_jules_hydrology <- function(
    monthly_files,
    model_template,
    catch_mask,
    year,
    expected_days
) {
  stop_missing(monthly_files, "JULES monthly hydrology NetCDF")
  
  vars <- c(
    JULES_RUNOFF_VAR,
    JULES_QUICKFLOW_VAR,
    JULES_SUBSURF_VAR,
    JULES_DRAIN_VAR,
    JULES_EVAP_VAR
  )
  
  clear_daily_outputs(RUNOFF_DIR_OUT, "runoff")
  clear_daily_outputs(QUICKFLOW_DIR_OUT, "quickflow")
  clear_daily_outputs(EVAP_DIR_OUT, "evap")
  clear_daily_outputs(JULES_SUBSURF_DIR_OUT, "sub_surf_roff")
  clear_daily_outputs(JULES_DRAIN_DIR_OUT, "drain")
  
  nc0 <- nc_open(monthly_files[1])
  on.exit(try(nc_close(nc0), silent = TRUE), add = TRUE)
  
  missing_vars <- vars[!vars %in% names(nc0$var)]
  if (length(missing_vars) > 0) {
    stop(
      "JULES output is missing required variable(s): ",
      paste(missing_vars, collapse = ", "), call. = FALSE
    )
  }
  
  for (v in vars) check_jules_flux_units(nc0, v)
  
  geom <- infer_jules_bng_geometry(nc0)
  nc_close(nc0)
  
  cat("\nJULES HYDROLOGY GRID\n")
  cat("====================\n")
  cat("Native rows/cols: ", geom$ny, " x ", geom$nx, "\n", sep = "")
  cat(sprintf("Inferred BNG resolution: %.3f x %.3f m\n", geom$dx, geom$dy))
  cat("Runoff mapping:    runoff -> runoff\n")
  cat("Quickflow mapping: surf_roff -> quickflow\n")
  cat("Evap mapping:      positive fqw_gb -> evap\n")
  cat("Extra QA outputs:  sub_surf_roff, drain\n\n")
  
  expected_dates <- seq.Date(
    as.Date(sprintf("%04d-01-01", year)),
    as.Date(sprintf("%04d-12-31", year)),
    by = "day"
  )
  
  if (length(expected_dates) != expected_days) {
    stop("Internal expected-day count mismatch.", call. = FALSE)
  }
  
  month_expected <- as.integer(table(format(expected_dates, "%m")))
  day_counter <- 0L
  summary_rows <- vector("list", expected_days)
  
  for (m in seq_along(monthly_files)) {
    f <- monthly_files[m]
    nc <- nc_open(f)
    
    ntime <- if ("time" %in% names(nc$dim)) nc$dim$time$len else NA_integer_
    if (is.na(ntime)) {
      nc_close(nc)
      stop("No time dimension in JULES file: ", f, call. = FALSE)
    }
    
    if (ntime != month_expected[m]) {
      nc_close(nc)
      stop(
        "JULES file ", basename(f), " has ", ntime,
        " records; expected ", month_expected[m], ".",
        call. = FALSE
      )
    }
    
    for (d in seq_len(ntime)) {
      day_counter <- day_counter + 1L
      
      r_runoff <- jules_matrix_to_raster(
        nc_read_daily_yx(nc, JULES_RUNOFF_VAR, d), geom
      )
      r_surf <- jules_matrix_to_raster(
        nc_read_daily_yx(nc, JULES_QUICKFLOW_VAR, d), geom
      )
      r_sub <- jules_matrix_to_raster(
        nc_read_daily_yx(nc, JULES_SUBSURF_VAR, d), geom
      )
      r_drain <- jules_matrix_to_raster(
        nc_read_daily_yx(nc, JULES_DRAIN_VAR, d), geom
      )
      r_fqw <- jules_matrix_to_raster(
        nc_read_daily_yx(nc, JULES_EVAP_VAR, d), geom
      )
      
      daily_stack <- c(r_runoff, r_surf, r_sub, r_drain, r_fqw)
      names(daily_stack) <- c("runoff", "surf", "sub", "drain", "fqw")
      
      daily_2km <- jules_to_colne_2km(
        daily_stack,
        model_template,
        catch_mask
      )
      
      # Convert daily-mean rates from kg m-2 s-1 to m/day.
      daily_2km <- daily_2km * JULES_FLUX_TO_M_DAY
      
      runoff <- ifel(daily_2km[["runoff"]] < 0, 0, daily_2km[["runoff"]])
      quickflow <- ifel(daily_2km[["surf"]] < 0, 0, daily_2km[["surf"]])
      subsurf <- ifel(daily_2km[["sub"]] < 0, 0, daily_2km[["sub"]])
      drain <- ifel(daily_2km[["drain"]] < 0, 0, daily_2km[["drain"]])
      
      # JULES fqw_gb is the moisture flux FROM the surface. Positive values
      # represent evaporation/sublimation loss. Negative values represent
      # condensation/deposition and are clipped to zero because the current
      # P-FASE 'evap' input is a non-negative water loss term.
      fqw_signed <- daily_2km[["fqw"]]
      evap <- ifel(fqw_signed < 0, 0, fqw_signed)
      
      ti <- day_counter - 1L
      
      write_float(runoff, file.path(RUNOFF_DIR_OUT, sprintf("runoff_%d.tif", ti)))
      write_float(quickflow, file.path(QUICKFLOW_DIR_OUT, sprintf("quickflow_%d.tif", ti)))
      write_float(evap, file.path(EVAP_DIR_OUT, sprintf("evap_%d.tif", ti)))
      write_float(subsurf, file.path(JULES_SUBSURF_DIR_OUT, sprintf("sub_surf_roff_%d.tif", ti)))
      write_float(drain, file.path(JULES_DRAIN_DIR_OUT, sprintf("drain_%d.tif", ti)))
      
      residual <- runoff - (quickflow + subsurf)
      residual_abs <- abs(residual)
      
      sr <- stats4(runoff)
      sq <- stats4(quickflow)
      ss <- stats4(subsurf)
      sd <- stats4(drain)
      se <- stats4(evap)
      sres <- stats4(residual_abs)
      
      neg_fqw_n <- global(fqw_signed < 0, "sum", na.rm = TRUE)[1, 1]
      
      summary_rows[[day_counter]] <- data.frame(
        timestep = ti,
        date = as.character(expected_dates[day_counter]),
        runoff_mean_m_day = sr["mean"],
        runoff_max_m_day = sr["max"],
        quickflow_mean_m_day = sq["mean"],
        subsurface_mean_m_day = ss["mean"],
        drain_mean_m_day = sd["mean"],
        evap_mean_m_day = se["mean"],
        runoff_identity_max_abs_m_day = sres["max"],
        negative_fqw_cells_clipped = neg_fqw_n
      )
      
      if (day_counter == 1L || day_counter %% 50L == 0L || day_counter == expected_days) {
        cat(sprintf("  JULES hydrology %d/%d\n", day_counter, expected_days))
      }
    }
    
    nc_close(nc)
  }
  
  if (day_counter != expected_days) {
    stop(
      "Processed ", day_counter, " JULES days; expected ", expected_days, ".",
      call. = FALSE
    )
  }
  
  summary_df <- do.call(rbind, summary_rows)
  write.csv(
    summary_df,
    file.path(QA_DIR, sprintf("JULES_%d_hydrology_daily_summary.csv", year)),
    row.names = FALSE
  )
  
  max_identity_error <- max(
    summary_df$runoff_identity_max_abs_m_day,
    na.rm = TRUE
  )
  
  cat(sprintf(
    "JULES runoff identity QA: max |runoff-(surface+subsurface)| = %.6g m/day\n",
    max_identity_error
  ))
  
  invisible(TRUE)
}

nc_find_var <- function(nc, candidates, required = TRUE) {
  vars <- names(nc$var)
  hit <- candidates[candidates %in% vars]
  
  if (length(hit) == 0) {
    if (required) {
      stop(
        "None of these NetCDF variables were found: ",
        paste(candidates, collapse = ", "),
        "\nAvailable variables: ",
        paste(vars, collapse = ", "),
        call. = FALSE
      )
    }
    return(NA_character_)
  }
  
  hit[1]
}

# Construct a correctly oriented SpatRaster from one CEH-GEAR x,y matrix.
gear_matrix_to_raster <- function(z_xy, x, y) {
  dx <- median(abs(diff(x)))
  dy <- median(abs(diff(y)))
  
  # ncvar_get returns [x, y]. Convert to [row(y), col(x)].
  m <- t(z_xy)
  
  # terra row 1 is north/top. Ensure rows run north -> south.
  if (length(y) > 1 && y[1] < y[length(y)]) {
    m <- m[nrow(m):1, , drop = FALSE]
  }
  
  # Ensure columns run west -> east.
  if (length(x) > 1 && x[1] > x[length(x)]) {
    m <- m[, ncol(m):1, drop = FALSE]
  }
  
  r <- rast(
    m,
    extent = ext(
      min(x) - dx / 2,
      max(x) + dx / 2,
      min(y) - dy / 2,
      max(y) + dy / 2
    ),
    crs = TARGET_CRS
  )
  
  r
}

process_gear_daily <- function(
    gear_nc,
    out_dir,
    model_template,
    catch_mask,
    year,
    expected_days
) {
  stop_missing(gear_nc, "CEH-GEAR daily NetCDF")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  nc <- nc_open(gear_nc)
  on.exit(nc_close(nc), add = TRUE)
  
  rain_var <- nc_find_var(nc, GEAR_RAIN_VAR_CANDIDATES)
  
  cat("\nCEH-GEAR rainfall variable: ", rain_var, "\n", sep = "")
  cat("NetCDF dimensions: ", paste(names(nc$dim), collapse = ", "), "\n", sep = "")
  
  units_att <- ncatt_get(nc, rain_var, "units")
  rain_units <- if (isTRUE(units_att$hasatt)) units_att$value else NA_character_
  
  cat("CEH-GEAR rainfall units attribute: ", rain_units, "\n", sep = "")
  
  if (!is.na(rain_units) && !grepl("mm", rain_units, ignore.case = TRUE)) {
    stop(
      "Expected CEH-GEAR daily rainfall units containing 'mm', but found: ",
      rain_units,
      call. = FALSE
    )
  }
  
  # CEH-GEAR daily files use x, y and time dimensions.
  x_name <- if ("x" %in% names(nc$dim)) "x" else
    if ("Easting" %in% names(nc$dim)) "Easting" else NA_character_
  y_name <- if ("y" %in% names(nc$dim)) "y" else
    if ("Northing" %in% names(nc$dim)) "Northing" else NA_character_
  time_name <- if ("time" %in% names(nc$dim)) "time" else
    if ("Time" %in% names(nc$dim)) "Time" else NA_character_
  
  if (any(is.na(c(x_name, y_name, time_name)))) {
    stop(
      "Could not identify CEH-GEAR x/y/time dimensions. Found: ",
      paste(names(nc$dim), collapse = ", "),
      call. = FALSE
    )
  }
  
  x <- ncvar_get(nc, x_name)
  y <- ncvar_get(nc, y_name)
  ntime <- length(nc$dim[[time_name]]$vals)
  
  if (ntime != expected_days) {
    stop(
      "CEH-GEAR file has ", ntime, " daily layers, expected ",
      expected_days, " for ", year, ".",
      call. = FALSE
    )
  }
  
  # Establish dimension order of rainfall_amount dynamically.
  v <- nc$var[[rain_var]]
  dim_names <- vapply(v$dim, function(d) d$name, character(1))
  
  ix <- match(x_name, dim_names)
  iy <- match(y_name, dim_names)
  it <- match(time_name, dim_names)
  
  if (any(is.na(c(ix, iy, it)))) {
    stop(
      "rainfall variable dimensions are ",
      paste(dim_names, collapse = ", "),
      "; expected x, y, time.",
      call. = FALSE
    )
  }
  
  # We read one day at a time, so the large annual file is never loaded fully.
  start <- rep(1L, length(dim_names))
  count <- vapply(v$dim, function(d) d$len, integer(1))
  count[it] <- 1L
  
  annual_total <- model_template
  values(annual_total) <- 0
  
  daily_summary <- vector("list", expected_days)
  
  for (i in seq_len(expected_days)) {
    start[it] <- i
    
    z <- ncvar_get(
      nc,
      rain_var,
      start = start,
      count = count,
      collapse_degen = TRUE
    )
    
    # After removing the singleton time dimension the remaining array can be
    # [x,y] or [y,x]. Re-orient to [x,y] if necessary.
    rem_names <- dim_names[-it]
    
    if (length(dim(z)) != 2) {
      stop("Unexpected rainfall slice dimensionality on day ", i, call. = FALSE)
    }
    
    if (identical(rem_names, c(x_name, y_name))) {
      z_xy <- z
    } else if (identical(rem_names, c(y_name, x_name))) {
      z_xy <- t(z)
    } else {
      stop(
        "Unexpected spatial dimension order: ",
        paste(rem_names, collapse = ", "),
        call. = FALSE
      )
    }
    
    r1 <- gear_matrix_to_raster(z_xy, x, y)
    
    # Rainfall is a DEPTH, therefore when going 1 km -> 2 km use the spatial
    # MEAN, not the sum.
    r2_mm <- project(
      r1,
      model_template,
      method = "average",
      threads = TRUE
    )
    r2_mm <- mask(r2_mm, catch_mask)
    
    # Each layer is a 24-hour rainfall amount in mm.
    r2_m <- r2_mm / 1000
    r2_m <- ifel(r2_m < 0, 0, r2_m)
    
    write_float(
      r2_m,
      file.path(out_dir, sprintf("precip_%d.tif", i - 1L))
    )
    
    annual_total <- annual_total + r2_m
    
    s <- stats4(r2_m)
    daily_summary[[i]] <- data.frame(
      timestep = i - 1L,
      date = as.character(as.Date(sprintf("%d-01-01", year)) + (i - 1L)),
      min_m_day = s["min"],
      median_m_day = s["median"],
      mean_m_day = s["mean"],
      max_m_day = s["max"]
    )
    
    if (i == 1 || i %% 50 == 0 || i == expected_days) {
      cat(sprintf("  CEH-GEAR precipitation %d/%d\n", i, expected_days))
    }
  }
  
  write_float(
    annual_total,
    file.path(QA_DIR, sprintf("CEH_GEAR_%d_annual_total_m_2km.tif", year))
  )
  
  write.csv(
    do.call(rbind, daily_summary),
    file.path(QA_DIR, sprintf("CEH_GEAR_%d_daily_summary.csv", year)),
    row.names = FALSE
  )
  
  show_stats("CEH-GEAR annual rainfall [m/year]", annual_total)
  
  invisible(list(variable = rain_var, units = rain_units))
}

check_gear_missing_records <- function(
    missing_nc,
    model_template,
    catch_mask
) {
  if (!file.exists(missing_nc)) {
    warning(
      "CEH-GEAR missing-records QA file not found: ", missing_nc,
      ". Rainfall can still be processed, but annual missing-record QA is skipped."
    )
    return(invisible(NULL))
  }
  
  nc <- nc_open(missing_nc)
  on.exit(nc_close(nc), add = TRUE)
  
  miss_var <- nc_find_var(nc, GEAR_MISSING_VAR_CANDIDATES, required = FALSE)
  
  if (is.na(miss_var)) {
    warning(
      "Could not identify missing-record variable in ", missing_nc,
      ". Variables: ", paste(names(nc$var), collapse = ", ")
    )
    return(invisible(NULL))
  }
  
  x_name <- if ("x" %in% names(nc$dim)) "x" else
    if ("Easting" %in% names(nc$dim)) "Easting" else NA_character_
  y_name <- if ("y" %in% names(nc$dim)) "y" else
    if ("Northing" %in% names(nc$dim)) "Northing" else NA_character_
  
  if (is.na(x_name) || is.na(y_name)) {
    warning("Could not identify x/y in missing-record file.")
    return(invisible(NULL))
  }
  
  x <- ncvar_get(nc, x_name)
  y <- ncvar_get(nc, y_name)
  v <- nc$var[[miss_var]]
  z <- ncvar_get(nc, miss_var, collapse_degen = TRUE)
  
  # Remove singleton time dimension if present.
  while (length(dim(z)) > 2 && any(dim(z) == 1)) {
    singleton <- which(dim(z) == 1)[1]
    z <- drop(z)
  }
  
  if (length(dim(z)) != 2) {
    warning("Unexpected missing-record data shape; QA skipped.")
    return(invisible(NULL))
  }
  
  dim_names <- vapply(v$dim, function(d) d$name, character(1))
  spatial_names <- dim_names[dim_names %in% c(x_name, y_name)]
  
  if (length(spatial_names) == 2 && identical(spatial_names, c(y_name, x_name))) {
    z <- t(z)
  }
  
  r1 <- gear_matrix_to_raster(z, x, y)
  
  # Use MAX when aggregating missing counts: a 2 km cell is flagged if any
  # contributing 1 km cell has missing estimates.
  r2 <- project(r1, model_template, method = "max", threads = TRUE)
  r2 <- mask(r2, catch_mask)
  
  write_float(
    r2,
    file.path(QA_DIR, "CEH_GEAR_2019_missing_records_max_2km.tif")
  )
  
  mx <- global(r2, "max", na.rm = TRUE)[1, 1]
  
  if (!is.na(mx) && mx > 0) {
    warning(
      "CEH-GEAR has missing daily estimates within the Colne model domain. ",
      "Maximum missing-record count = ", mx,
      ". Inspect qa/CEH_GEAR_2019_missing_records_max_2km.tif."
    )
  } else {
    message("CEH-GEAR missing-record QA: no missing daily estimates in Colne.")
  }
}

# ------------------------------------------------------------------------------
# 1. SELECT COLNE AND BUILD EXACT 2 km BNG GRID
# ------------------------------------------------------------------------------

stop_missing(CATCHMENT_SHP, "Colne catchment shapefile")

catch_all <- st_read(CATCHMENT_SHP, quiet = TRUE)

if (is.na(st_crs(catch_all))) {
  stop("Catchment shapefile has no CRS.", call. = FALSE)
}

catch_all <- st_make_valid(catch_all)
attrs <- st_drop_geometry(catch_all)

char_cols <- names(attrs)[vapply(
  attrs,
  function(x) is.character(x) || is.factor(x),
  logical(1)
)]

hit <- rep(FALSE, nrow(catch_all))
matching_fields <- character()

for (nm in char_cols) {
  h <- grepl(
    CATCHMENT_PATTERN,
    as.character(attrs[[nm]]),
    ignore.case = TRUE
  )
  h[is.na(h)] <- FALSE
  
  if (any(h)) {
    hit <- hit | h
    matching_fields <- c(matching_fields, nm)
  }
}

if (any(hit)) {
  catch_sel <- catch_all[hit, ]
  message(
    "Selected ", sum(hit), " Colne feature(s); matching fields: ",
    paste(unique(matching_fields), collapse = ", ")
  )
} else if (nrow(catch_all) == 1) {
  catch_sel <- catch_all
  warning("No 'Colne' name match, but shapefile has one feature; using it.")
} else if (ALLOW_UNION_ALL_IF_NO_NAME_MATCH) {
  catch_sel <- catch_all
  warning("No 'Colne' name match; unioning all features by explicit override.")
} else {
  for (nm in char_cols) {
    cat("\nFIELD ", nm, "\n", sep = "")
    print(head(unique(as.character(attrs[[nm]])), 40))
  }
  
  stop(
    "Could not identify the Colne feature. Set CATCHMENT_PATTERN to the exact name.",
    call. = FALSE
  )
}

catch_sf <- st_sf(
  catchment = "Colne",
  geometry = st_union(st_geometry(catch_sel))
)

catch_sf <- st_make_valid(catch_sf)
catch_sf <- st_transform(catch_sf, TARGET_CRS)
catch_v <- vect(catch_sf)

st_write(
  catch_sf,
  file.path(QA_DIR, "colne_catchment_used.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

grid_2km <- make_2km_template(catch_sf)

catch_mask <- rasterize(
  catch_v,
  grid_2km,
  field = 1,
  background = NA,
  touches = TRUE
)

write_int(
  catch_mask,
  file.path(MODEL_DIR, "colne_mask_2km.tif"),
  "INT1U"
)

n_valid <- global(!is.na(catch_mask), "sum", na.rm = TRUE)[1, 1]

cat("\nCOLNE MODEL GRID\n")
cat("================\n")
cat("CRS:        EPSG:27700\n")
cat("Resolution: 2000 x 2000 m\n")
cat("Rows/cols:  ", nrow(grid_2km), " x ", ncol(grid_2km), "\n", sep = "")
cat("Valid cells:", n_valid, "\n")

# ------------------------------------------------------------------------------
# 2. DEM FOR NEW 2 km ROUTING
# ------------------------------------------------------------------------------

routing_dem_ready <- FALSE

if (length(HYDROSHEDS_DEM) == 0 || is.na(HYDROSHEDS_DEM)) {
  stop(
    paste0(
      "HydroSHEDS source DEM not found.\n",
      "Expected one of:\n  ",
      paste(HYDROSHEDS_DEM_CANDIDATES, collapse = "\n  "),
      "\n\nCopy hydrosheds_eu_dem_15s.tif to one of these locations."
    ),
    call. = FALSE
  )
}

message("Using HydroSHEDS source DEM: ", HYDROSHEDS_DEM)

dem_src <- rast(HYDROSHEDS_DEM)

if (is.na(crs(dem_src)) || crs(dem_src) == "") {
  stop(
    "HydroSHEDS DEM has no CRS metadata: ",
    HYDROSHEDS_DEM,
    call. = FALSE
  )
}

cat("\nHYDROSHEDS SOURCE DEM\n")
cat("=====================\n")
cat("File:       ", HYDROSHEDS_DEM, "\n")
cat("CRS:        ", crs(dem_src, proj = TRUE), "\n")
cat("Resolution: ", paste(res(dem_src), collapse = " x "), "\n")
cat("Rows/cols:  ", nrow(dem_src), " x ", ncol(dem_src), "\n")
cat("Extent:     ", paste(as.vector(ext(dem_src)), collapse = ", "), "\n\n")

# Crop in source CRS before reprojection.
catch_buffer_sf <- st_buffer(
  catch_sf,
  dist = DEM_BUFFER_M
)

catch_buffer_dem <- project(
  vect(catch_buffer_sf),
  crs(dem_src)
)

dem_crop <- crop(
  dem_src,
  catch_buffer_dem,
  snap = "out"
)

dem_crop <- mask(
  dem_crop,
  catch_buffer_dem,
  touches = TRUE
)

if (all(is.na(values(dem_crop)))) {
  stop(
    "HydroSHEDS crop contains no valid cells. Check catchment and DEM CRS.",
    call. = FALSE
  )
}

# Preserve native-resolution Colne + buffer DEM for QA.
writeRaster(
  dem_crop,
  file.path(QA_DIR, "hydrosheds_colne_buffer_native.tif"),
  overwrite = TRUE,
  datatype = "INT2S",
  NAflag = -32768,
  gdal = c("COMPRESS=DEFLATE")
)

# Directly aggregate/project from native HydroSHEDS to the exact 2 km BNG grid.
dem_2km <- project(
  dem_crop,
  grid_2km,
  method = "average",
  threads = TRUE
)

dem_2km <- mask(
  dem_2km,
  catch_mask
)

valid_dem <- global(
  !is.na(dem_2km),
  "sum",
  na.rm = TRUE
)[1, 1]

if (is.na(valid_dem) || valid_dem == 0) {
  stop(
    "No valid HydroSHEDS elevation cells overlap the Colne 2 km grid.",
    call. = FALSE
  )
}

write_float(
  dem_2km,
  file.path(MODEL_DIR, "dem_2km.tif")
)

show_stats(
  "HydroSHEDS DEM on Colne 2 km grid [m]",
  dem_2km
)

# Old 10 km DEM retained only for QA comparison.
if (file.exists(OLD_DEM_10KM)) {
  old_dem_qa <- to_2km(
    OLD_DEM_10KM,
    grid_2km,
    catch_mask,
    method = "bilinear"
  )
  
  write_float(
    old_dem_qa,
    file.path(QA_DIR, "old_10km_dem_interpolated_2km_QA_ONLY.tif")
  )
  
  dem_difference <- dem_2km - old_dem_qa
  
  write_float(
    dem_difference,
    file.path(QA_DIR, "hydrosheds2km_minus_old10km_dem_QA.tif")
  )
}

routing_dem_ready <- TRUE

# ------------------------------------------------------------------------------
# 3. SOIL: ORIGINAL DEPTH LAYERS -> 0-30 cm -> COLNE 2 km
# ------------------------------------------------------------------------------

all_raw_soil <- unlist(lapply(
  soil_depth_files,
  function(x) file.path(SOIL_RAW_DIR, x)
))

use_raw_soil <- all(file.exists(all_raw_soil))

if (use_raw_soil) {
  message("Using original depth-resolved soil layers from input_data/.")
  
  raw <- list()
  
  for (prop in names(soil_depth_files)) {
    raw[[prop]] <- depth_weighted_mean(
      file.path(SOIL_RAW_DIR, soil_depth_files[[prop]]),
      depth_weights
    )
  }
  
  # SoilGrids source encodings
  raw$sand <- raw$sand * 0.001
  raw$clay <- raw$clay * 0.001
  raw$silt <- raw$silt * 0.001
  raw$cfvo <- raw$cfvo * 0.001
  
  # SoilGrids bdod numerical representation:
  # e.g. raw 130 -> 1300 kg/m3.
  raw$bdod <- raw$bdod * 10
  
  sand <- project(raw$sand, grid_2km, method = "average", threads = TRUE)
  clay <- project(raw$clay, grid_2km, method = "average", threads = TRUE)
  silt <- project(raw$silt, grid_2km, method = "average", threads = TRUE)
  cfvol <- project(raw$cfvo, grid_2km, method = "average", threads = TRUE)
  bd <- project(raw$bdod, grid_2km, method = "average", threads = TRUE)
  
  sand <- mask(sand, catch_mask)
  clay <- mask(clay, catch_mask)
  silt <- mask(silt, catch_mask)
  cfvol <- mask(cfvol, catch_mask)
  bd <- mask(bd, catch_mask)
  
} else {
  warning(
    "Depth-resolved soil files are incomplete; using existing coarse 0-30 cm rasters."
  )
  
  stop_missing(unlist(soil_fallback), "fallback soil raster")
  
  sand <- to_2km(soil_fallback$sand, grid_2km, catch_mask, "bilinear")
  clay <- to_2km(soil_fallback$clay, grid_2km, catch_mask, "bilinear")
  silt <- to_2km(soil_fallback$silt, grid_2km, catch_mask, "bilinear")
  cfvol <- to_2km(soil_fallback$cfvo, grid_2km, catch_mask, "bilinear")
  
  bd_old <- to_2km(soil_fallback$bdod, grid_2km, catch_mask, "bilinear")
  bd_med <- stats4(bd_old)["median"]
  
  if (bd_med > 50000 && bd_med < 500000) {
    bd <- bd_old / 100
  } else if (bd_med >= 500 && bd_med <= 3000) {
    bd <- bd_old
  } else {
    stop(
      "Cannot safely infer existing bulk-density units; median = ",
      bd_med,
      call. = FALSE
    )
  }
}

# Fine-earth fractions must sum to 1.
texture_total <- sand + clay + silt

sand <- ifel(texture_total > 0, sand / texture_total, NA)
clay <- ifel(texture_total > 0, clay / texture_total, NA)
silt <- ifel(texture_total > 0, silt / texture_total, NA)

write_float(bd, file.path(MODEL_DIR, "soil_bulk_density.tif"))
write_float(sand, file.path(MODEL_DIR, "soil_texture_sand_content.tif"))
write_float(clay, file.path(MODEL_DIR, "soil_texture_clay_content.tif"))
write_float(silt, file.path(MODEL_DIR, "soil_texture_silt_content.tif"))

write_float(
  cfvol,
  file.path(QA_DIR, "soil_coarse_fragment_VOLUME_fraction_2km.tif")
)

coarse_model <- catch_mask * COARSE_FRAG_MODEL_VALUE

write_float(
  coarse_model,
  file.path(MODEL_DIR, "soil_texture_coarse_frag_content.tif")
)

show_stats("bulk density [kg/m3]", bd)
show_stats("sand [kg/kg]", sand)
show_stats("clay [kg/kg]", clay)
show_stats("silt [kg/kg]", silt)
show_stats("sand+clay+silt", sand + clay + silt)

# ------------------------------------------------------------------------------
# 4. HYDRAULIC SOIL PROPERTIES
# ------------------------------------------------------------------------------

stop_missing(
  c(
    static_src$field_capacity,
    static_src$saturation,
    static_src$hydraulic_conductivity
  ),
  "soil hydraulic raster"
)

fc <- to_2km(
  static_src$field_capacity, grid_2km, catch_mask, "bilinear"
) / 100

sat <- to_2km(
  static_src$saturation, grid_2km, catch_mask, "bilinear"
) / 100

ks <- to_2km(
  static_src$hydraulic_conductivity, grid_2km, catch_mask, "bilinear"
)

ks <- ks_encoded_to_m_s(ks)

fc <- clamp(fc, 0, 1, values = TRUE)
sat <- clamp(sat, 0, 1, values = TRUE)
fc <- ifel(fc > sat, sat, fc)
ks <- ifel(ks < 0, 0, ks)

write_float(
  fc,
  file.path(MODEL_DIR, "soil_water_content_field_capacity.tif")
)
write_float(
  sat,
  file.path(MODEL_DIR, "soil_water_content_saturation.tif")
)
write_float(
  ks,
  file.path(MODEL_DIR, "soil_hydraulic_conductivity.tif")
)

show_stats("field capacity [cm3/cm3]", fc)
show_stats("saturation [cm3/cm3]", sat)
show_stats("hydraulic conductivity [m/s]", ks)

# ------------------------------------------------------------------------------
# 5. USLE C, LS AND TEMPORARY P
# ------------------------------------------------------------------------------

stop_missing(
  c(static_src$c_factor, static_src$ls_factor),
  "USLE source raster"
)

c_factor <- to_2km(
  static_src$c_factor, grid_2km, catch_mask, "bilinear"
)
ls_factor <- to_2km(
  static_src$ls_factor, grid_2km, catch_mask, "bilinear"
)

if (stats4(c_factor)["min"] < 0 || stats4(c_factor)["max"] > 1.001) {
  stop("C factor is outside expected [0,1] range.", call. = FALSE)
}

if (stats4(ls_factor)["min"] < 0) {
  stop("LS factor contains negative values.", call. = FALSE)
}

p_factor <- catch_mask * P_FACTOR_VALUE

write_float(c_factor, file.path(MODEL_DIR, "soil_usle_c_factor.tif"))
write_float(ls_factor, file.path(MODEL_DIR, "soil_usle_ls_factor.tif"))
write_float(p_factor, file.path(MODEL_DIR, "soil_usle_p_factor.tif"))

# ------------------------------------------------------------------------------
# 6. UKCEH LCM2024 25 m LAND COVER
# ------------------------------------------------------------------------------

stop_missing(LCM2024, "UKCEH LCM2024 25 m raster")

lcm <- rast(LCM2024)

if (nlyr(lcm) < 1) {
  stop("LCM2024 contains no raster bands.", call. = FALSE)
}

if (is.na(crs(lcm)) || crs(lcm) == "") {
  stop("LCM2024 has no CRS.", call. = FALSE)
}

# Product band 1 = class identifier.
lcm_class <- lcm[[1]]

catch_lcm <- project(catch_v, crs(lcm_class))
lcm_colne <- crop(lcm_class, catch_lcm, snap = "out")
lcm_colne <- mask(lcm_colne, catch_lcm, touches = TRUE)
lcm_colne <- round(lcm_colne)

# Preserve 25 m resolution for nfdata fractional land-use aggregation.
write_int(
  lcm_colne,
  file.path(MODEL_DIR, "land_use_lcm2024_25m_colne.tif"),
  "INT1U"
)

# Optional LCM QA bands.
if (nlyr(lcm) >= 2) {
  probability <- crop(lcm[[2]], catch_lcm, snap = "out")
  probability <- mask(probability, catch_lcm, touches = TRUE)
  write_int(
    probability,
    file.path(QA_DIR, "LCM2024_probability_25m_colne.tif"),
    "INT1U"
  )
}

if (nlyr(lcm) >= 3) {
  purity <- crop(lcm[[3]], catch_lcm, snap = "out")
  purity <- mask(purity, catch_lcm, touches = TRUE)
  write_int(
    purity,
    file.path(QA_DIR, "LCM2024_purity_25m_colne.tif"),
    "INT1U"
  )
}

freq_lcm <- freq(lcm_colne, digits = 0)
write.csv(
  as.data.frame(freq_lcm),
  file.path(QA_DIR, "LCM2024_Colne_class_counts.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. FRESHWATER ESTUARY FLAG
# ------------------------------------------------------------------------------

if (USE_ESTUARY) {
  stop(
    "This first configuration is not parameterised for the Colne estuary. ",
    "Set USE_ESTUARY=FALSE until estuary geometry/constants are prepared.",
    call. = FALSE
  )
}

is_estuary <- catch_mask * 0

write_int(
  is_estuary,
  file.path(MODEL_DIR, "is_estuary_2km.tif"),
  "INT1U"
)

# ------------------------------------------------------------------------------
# 8. CEH-GEAR 2019 ONE-NETCDF -> 365 COLNE 2 km PRECIP TIFFS
# ------------------------------------------------------------------------------

process_gear_daily(
  GEAR_NC,
  PRECIP_DIR,
  grid_2km,
  catch_mask,
  SIM_YEAR,
  N_TIMESTEPS
)

check_gear_missing_records(
  GEAR_MISSING_NC,
  grid_2km,
  catch_mask
)

# ------------------------------------------------------------------------------
# 9. JULES 2019 HYDROLOGY -> COLNE 2 km DAILY TIFFS
# ------------------------------------------------------------------------------

matched_hydro_ready <- process_jules_hydrology(
  JULES_MONTHLY_FILES,
  grid_2km,
  catch_mask,
  SIM_YEAR,
  N_TIMESTEPS
)

# ------------------------------------------------------------------------------
# 10. CANDIDATE SEDIMENT CALIBRATION LAYERS
# ------------------------------------------------------------------------------

candidate_names <- c(
  "bank_erosion_alpha",
  "bank_erosion_beta",
  "deposition_alpha",
  "deposition_beta",
  "resuspension_alpha",
  "resuspension_beta",
  "sediment_transport_a",
  "sediment_transport_b",
  "sediment_transport_c"
)

for (nm in candidate_names) {
  p <- static_src[[nm]]
  
  if (!file.exists(p)) {
    warning("Candidate calibration source missing: ", p)
    next
  }
  
  r <- to_2km(p, grid_2km, catch_mask, "bilinear")
  
  write_float(
    r,
    file.path(CALIB_DIR, paste0(nm, "_2km.tif"))
  )
}

# ------------------------------------------------------------------------------
# 11. PFAS TEST EMISSION
# ------------------------------------------------------------------------------

pfas_test <- catch_mask * PFAS_TEST_LOAD_KG_M2_DAY

write_float(
  pfas_test,
  file.path(MODEL_DIR, "pfas_test_emission_water.tif")
)

# ------------------------------------------------------------------------------
# 12. FINAL QA
# ------------------------------------------------------------------------------

required_2km <- c(
  "soil_bulk_density.tif",
  "soil_texture_sand_content.tif",
  "soil_texture_clay_content.tif",
  "soil_texture_silt_content.tif",
  "soil_texture_coarse_frag_content.tif",
  "soil_water_content_field_capacity.tif",
  "soil_water_content_saturation.tif",
  "soil_hydraulic_conductivity.tif",
  "soil_usle_c_factor.tif",
  "soil_usle_ls_factor.tif",
  "soil_usle_p_factor.tif",
  "is_estuary_2km.tif",
  "pfas_test_emission_water.tif"
)

if (routing_dem_ready) {
  required_2km <- c("dem_2km.tif", required_2km)
}

qa_rows <- list()

for (fn in required_2km) {
  r <- rast(file.path(MODEL_DIR, fn))
  
  geom_ok <- compareGeom(
    r,
    grid_2km,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  )
  
  s <- stats4(r)
  
  qa_rows[[fn]] <- data.frame(
    file = fn,
    exact_2km_geometry = geom_ok,
    nrow = nrow(r),
    ncol = ncol(r),
    res_x = res(r)[1],
    res_y = res(r)[2],
    min = s["min"],
    median = s["median"],
    mean = s["mean"],
    max = s["max"]
  )
}

qa_df <- do.call(rbind, qa_rows)

write.csv(
  qa_df,
  file.path(QA_DIR, "static_model_raster_QA.csv"),
  row.names = FALSE
)

if (!all(qa_df$exact_2km_geometry)) {
  stop(
    "One or more ordinary model rasters do not exactly match the 2 km grid.",
    call. = FALSE
  )
}

# Confirm all 365 JULES-derived time-step rasters exist.
hydro_counts <- c(
  runoff = length(list.files(RUNOFF_DIR_OUT, pattern = "^runoff_[0-9]+\\.tif$")),
  quickflow = length(list.files(QUICKFLOW_DIR_OUT, pattern = "^quickflow_[0-9]+\\.tif$")),
  evap = length(list.files(EVAP_DIR_OUT, pattern = "^evap_[0-9]+\\.tif$")),
  sub_surf_roff = length(list.files(JULES_SUBSURF_DIR_OUT, pattern = "^sub_surf_roff_[0-9]+\\.tif$")),
  drain = length(list.files(JULES_DRAIN_DIR_OUT, pattern = "^drain_[0-9]+\\.tif$"))
)

if (any(hydro_counts != N_TIMESTEPS)) {
  stop(
    "JULES-derived daily raster count is incomplete: ",
    paste(names(hydro_counts), hydro_counts, sep = "=", collapse = ", "),
    call. = FALSE
  )
}

write.csv(
  data.frame(variable = names(hydro_counts), n_files = as.integer(hydro_counts)),
  file.path(QA_DIR, "JULES_daily_output_file_counts.csv"),
  row.names = FALSE
)

nfdata_ready <- routing_dem_ready && matched_hydro_ready

summary_lines <- c(
  "Colne P-FASE 2 km / CEH-GEAR / JULES / LCM2024 preprocessing",
  "====================================================",
  paste("Simulation year:", SIM_YEAR),
  paste("Target CRS:", TARGET_CRS),
  paste("Target resolution:", TARGET_RES, "m"),
  paste("Valid model cells:", n_valid),
  paste("CEH-GEAR source:", GEAR_NC),
  paste("LCM source:", LCM2024),
  paste("JULES hydrology source:", JULES_OUTPUT_DIR),
  paste("HydroSHEDS-derived 2 km routing DEM ready:", routing_dem_ready),
  paste("JULES 2019 runoff + quickflow + evap ready:", matched_hydro_ready),
  paste("Ready for nfdata create:", nfdata_ready),
  "",
  "Rainfall:",
  "  CEH-GEAR 1 km daily rainfall was spatially MEAN-aggregated to 2 km",
  "  and converted from mm/day to m/day.",
  "",
  "Land cover:",
  "  LCM2024 band 1 was clipped to Colne and retained at 25 m.",
  "  nfdata will map it to the 12 NanoFASE land-use fractions per 2 km cell.",
  "",
  if (!matched_hydro_ready) {
    paste0(
      "ACTION REQUIRED: JULES ", SIM_YEAR,
      " hydrology preprocessing did not complete."
    )
  } else {
    paste0(
      "JULES hydrological forcing is complete: runoff=runoff, ",
      "quickflow=surf_roff, evap=positive fqw_gb; sub_surf_roff and drain ",
      "retained as separate QA/future-PFAS products."
    )
  },
  if (!routing_dem_ready) {
    "ACTION REQUIRED: HydroSHEDS source DEM was not processed successfully."
  } else {
    "2 km DEM is ready for nfdata D8 routing derivation."
  }
)

writeLines(
  summary_lines,
  file.path(QA_DIR, "preprocessing_summary.txt")
)

cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
