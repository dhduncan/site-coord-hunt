# Check real inputs before using them. Run locally, in RStudio or from the repo root:
#
#   Rscript data-raw/check_inputs.R
#
# The report it prints holds counts and summaries only: no coordinates, record ids, place
# names or unit ids, and counts below `suppress_below` are shown as "<N". Read it against
# your data-sharing agreement before sharing it. The data itself never needs to leave your
# machine. (Set details = TRUE to see full error messages in your own console: they can
# contain ids or names, so do not share that output.)
#
# Edit the settings, run, and review the report. The check does not change any file.

suppressPackageStartupMessages({library(terra); library(sf)})
for (f in c("dates", "scoring", "candidates", "data_io", "check_inputs")) source(file.path("R", paste0(f, ".R")))

# ---- settings -----------------------------------------------------------------------
data_dir  <- "path/to/your/data"       # folder with covariates.tif and records.csv (outside the repo)
layers_dir <- NULL                      # optional folder of other rasters (land cover, dwellings, travel time, ...)
polygons <- list(                       # the approximate-location layers, as in build_candidates.R
  list(file = "path/to/parishes.gpkg",            id_col = "PARISH_ID", name_col = "PARISH_NAME",   prefix = "parish"),
  list(file = "path/to/facility_catchments.gpkg", id_col = "FAC_ID",    name_col = "FACILITY_NAME", prefix = "facility")
)
suppress_below <- 5                     # counts below this are shown as "<N"; 0 turns suppression off
time_scoring   <- TRUE                  # time one scoring run on the real mask
details        <- FALSE

# ---- run ------------------------------------------------------------------------------
if (!dir.exists(data_dir)) stop("Edit `data_dir` (and the other settings) in data-raw/check_inputs.R first.", call. = FALSE)
polygons <- Filter(function(l) file.exists(l$file), polygons)
check_inputs(data_dir, layers_dir, if (length(polygons)) polygons else NULL, suppress_below, time_scoring, details)
