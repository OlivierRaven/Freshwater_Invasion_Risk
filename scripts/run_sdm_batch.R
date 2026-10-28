# Runs analysis.qmd's "species-sdms-batch-all" chunk unattended, outside
# Quarto's render (which isn't suited to a multi-hour/multi-day background
# job). This script is a straight copy of that chunk's code plus the minimal
# setup it depends on (predictor stack, occ_sdm, species list) -- if you
# change the modeling logic, change it in analysis.qmd first and copy the
# change here, so the two never drift.
#
# Usage: Rscript scripts/run_sdm_batch.R  (safe to interrupt/resume any time
# -- already-successful species are skipped via data/derived/sdm_batch_results.csv)

setwd("/Users/olivierraven/Documents/Education:work/Freshwater_Invasion_Risk")
suppressMessages({
  library(tidyverse)
  library(terra)
  library(biomod2)
  library(randomForest)
  library(R.utils)
  library(dismo)
})

data_root <- Sys.getenv("FIW_DATA_ROOT", "/Volumes/Disk_OR/Freshwater_Invasion_Risk_data")
if (!dir.exists(data_root)) stop("Data drive not found at: ", data_root)
hydroatlas_dir <- file.path(data_root, "hydroatlas")
futurestreams_dir <- file.path(data_root, "futurestreams")

# ---- predictor stack (mirrors analysis.qmd "species-sdms-prototype") ----
keep11 <- c("UPLAND_SKM", "CATCH_SKM", "ORD_CLAS", "slp_dg_cav", "ele_mt_cav",
            "Q_mean", "Q_bfi_mean", "Q_zfw", "WT_mean", "WT_min", "WT_max")

hydro_vars <- c("UPLAND_SKM", "CATCH_SKM", "ORD_CLAS", "slp_dg_cav", "ele_mt_cav")
hydro_paths <- file.path(hydroatlas_dir, "rasterized", paste0(hydro_vars, ".tif"))
sdm_hydro_stack <- rast(hydro_paths)
names(sdm_hydro_stack) <- hydro_vars

fs_hist_files <- list.files(futurestreams_dir, pattern = "E2O_hist.*[.]nc$", full.names = TRUE)
fs_need_hyphen <- c("Q-mean", "Q-zfw", "WT-mean", "WT-min", "WT-max")
fs_need_clean  <- c("Q_mean", "Q_zfw", "WT_mean", "WT_min", "WT_max")
fs_hist_keep <- fs_hist_files[sapply(fs_need_hyphen, function(v) grep(paste0("^", v, "_"), basename(fs_hist_files)))]
sdm_fs_stack <- rast(fs_hist_keep)
names(sdm_fs_stack) <- fs_need_clean

bfi_file <- grep("^Q-bfi_", basename(fs_hist_files), value = TRUE)
Q_bfi_mean_r <- mean(rast(file.path(futurestreams_dir, bfi_file)))
names(Q_bfi_mean_r) <- "Q_bfi_mean"

ext(sdm_hydro_stack) <- ext(sdm_fs_stack)
sdm_predictor_stack <- c(sdm_hydro_stack, sdm_fs_stack, Q_bfi_mean_r)[[keep11]]
sdm_predictor_stack[["WT_max"]][sdm_predictor_stack[["WT_max"]] > 350] <- NA

occ_sdm <- read_csv("data/derived/occurrences_with_predictors.csv", show_col_types = FALSE) |>
  rename(Q_mean = `Q-mean`, Q_bfi_mean = `Q-bfi_mean`, Q_zfw = `Q-zfw`,
         WT_mean = `WT-mean`, WT_min = `WT-min`, WT_max = `WT-max`)

bulk_target_species <- read_csv("data/derived/bulk_353_species_target_list.csv", show_col_types = FALSE)

message(sprintf("Predictor stack ready (%d layers). %d species in target list.",
                 terra::nlyr(sdm_predictor_stack), nrow(bulk_target_species)))

# ---- batch loop (mirrors analysis.qmd "species-sdms-batch-all") ----
sdm_all_dir <- file.path(data_root, "outputs", "sdm_all")
dir.create(sdm_all_dir, showWarnings = FALSE, recursive = TRUE)

sdm_results_path <- "data/derived/sdm_batch_results.csv"
sdm_eval_path <- "data/derived/sdm_batch_evaluation.csv"
sdm_varimp_path <- "data/derived/sdm_batch_variable_importance.csv"

fit_and_save_one_sdm <- function(species_name, n_background = 10000, seed = 42) {
  clean_name <- gsub("[ -]+", "_", species_name)  # also strips hyphens (e.g. Alisma plantago-aquatica) -- same formula-parsing bug as pitfall 5, caught in the first full-batch smoke test when resp.name itself carried a hyphen

  presence <- occ_sdm |>
    filter(species_requested == species_name, if_all(all_of(keep11), ~ !is.na(.x)), WT_max <= 350) |>
    dplyr::select(decimalLongitude, decimalLatitude)

  set.seed(seed)
  background <- occ_sdm |>
    filter(species_requested != species_name) |>
    slice_sample(n = n_background) |>
    dplyr::select(decimalLongitude, decimalLatitude)

  all_pts <- bind_rows(presence, background)
  resp <- c(rep(1, nrow(presence)), rep(NA, nrow(background)))
  all_pts_vect <- vect(all_pts, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

  expl_values <- map(keep11, function(v) {
    e <- terra::extract(sdm_predictor_stack[[v]], all_pts_vect)[, 2]
    if (length(e) != nrow(all_pts)) stop("Row mismatch extracting ", v)
    e
  }) |> set_names(keep11) |> as_tibble()

  bm_data <- BIOMOD_FormatingData(
    resp.var = resp,
    resp.xy = as.matrix(all_pts[, c("decimalLongitude", "decimalLatitude")]),
    expl.var = expl_values,
    resp.name = clean_name,
    dir.name = sdm_all_dir,
    PA.nb.rep = 0
  )

  bm_models <- BIOMOD_Modeling(
    bm.format = bm_data,
    models = c("GLM", "RF"),
    CV.strategy = "kfold",
    CV.k = 2,
    CV.nb.rep = 1,
    OPT.strategy = "bigboss",
    var.import = 3,
    metric.eval = c("TSS", "AUCroc"),
    seed.val = seed
  )

  bm_proj <- BIOMOD_Projection(
    bm.mod = bm_models,
    proj.name = "current",
    new.env = sdm_predictor_stack,
    models.chosen = "all",
    build.clamping.mask = FALSE
  )

  suitability <- mean(get_predictions(bm_proj))
  # Flat filename, not a per-species subfolder -- biomod2 sanitizes resp.name
  # into its OWN internal folder name (e.g. replacing "_" with "."), which
  # doesn't match our clean_name, so a subfolder path built from clean_name
  # doesn't exist and writeRaster fails (caught in the first full-batch
  # smoke test, 2026-10-28). Writing flat into sdm_all_dir (guaranteed to
  # exist, created once before the loop) sidesteps that entirely.
  writeRaster(suitability, file.path(sdm_all_dir, paste0(clean_name, "_suitability_current.tif")), overwrite = TRUE)

  write_csv(get_evaluations(bm_models) |> mutate(species = species_name),
            sdm_eval_path, append = file.exists(sdm_eval_path))
  write_csv(get_variables_importance(bm_models) |> mutate(species = species_name),
            sdm_varimp_path, append = file.exists(sdm_varimp_path))

  rm(bm_data, bm_models, bm_proj, suitability); invisible(gc(verbose = FALSE))
  nrow(presence)
}

species_done <- if (file.exists(sdm_results_path)) {
  read_csv(sdm_results_path, show_col_types = FALSE) |> filter(status == "success") |> pull(species)
} else character(0)

species_to_run <- setdiff(bulk_target_species$species_requested, species_done)
message(sprintf("%d species already done, %d remaining", length(species_done), length(species_to_run)))

for (i in seq_along(species_to_run)) {
  sp <- species_to_run[i]
  t0 <- Sys.time()
  result <- tryCatch(
    {
      n_pres <- fit_and_save_one_sdm(sp)
      list(status = "success", error = NA_character_, n_presence = n_pres)
    },
    error = function(e) {
      message("FAILED: ", sp, " -- ", conditionMessage(e))
      list(status = "failed", error = conditionMessage(e), n_presence = NA_integer_)
    }
  )
  runtime_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  write_csv(
    tibble(species = sp, status = result$status, n_presence = result$n_presence,
           runtime_min = runtime_min, error = result$error,
           completed_at = as.character(Sys.time())),
    sdm_results_path, append = file.exists(sdm_results_path)
  )
  message(sprintf("[%d/%d] %s: %s (%.1f min)", i, length(species_to_run), sp, result$status, runtime_min))
}

message("=== BATCH COMPLETE ===")
