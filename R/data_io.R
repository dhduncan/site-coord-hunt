# The only file that reads or writes project data. Swap this for a database
# backend later without touching the app or scoring code.

CORRECTION_COLS <- c("action_id", "record_id", "cell_id", "lon", "lat", "weight",
                     "weight_mode", "method", "source", "confidence", "user",
                     "timestamp")
FLAG_COLS <- c("action", "record_id", "user", "timestamp", "note")

# ---- inputs -----------------------------------------------------------------

#' Read the project's inputs.
#'   covariates.tif  fixed model mask and covariates (cell IDs are stable)
#'   candidates.csv  unit_id, cell_id, weight   (built by data-raw/build_candidates.R)
#'   units.csv       unit_id, unit_name, layer, n_cells, mask_share
#'   records.csv     record_id, unit_id (or admin_name), lon, lat, n_specimens, date
#'   admin.gpkg      unit outlines for the map, with a unit_id column (optional)
load_inputs <- function(dir = "data") {
  need <- c("covariates.tif", "candidates.csv", "units.csv", "records.csv")
  missing <- need[!file.exists(file.path(dir, need))]
  if (length(missing))
    stop("Missing in '", dir, "': ", toString(missing),
         ". candidates.csv and units.csv come from data-raw/build_candidates.R.", call. = FALSE)
  covs <- terra::rast(file.path(dir, "covariates.tif"))
  units <- utils::read.csv(file.path(dir, "units.csv"), stringsAsFactors = FALSE,
                           colClasses = c(unit_id = "character", unit_name = "character", layer = "character"))
  candidates <- utils::read.csv(file.path(dir, "candidates.csv"), stringsAsFactors = FALSE,
                                colClasses = c(unit_id = "character"))
  check_candidates(candidates, units, covs)
  rec_file <- file.path(dir, "records.csv")
  cls <- c(record_id = "character", admin_name = "character", unit_id = "character",
           unit_layer = "character", date = "character")
  cls <- cls[names(cls) %in% names(utils::read.csv(rec_file, nrows = 0))]
  records <- utils::read.csv(rec_file, stringsAsFactors = FALSE, na.strings = c("", "NA"), colClasses = cls)
  if (is.null(records$date)) records$date <- NA_character_
  assigned <- assign_unit_ids(records, units)
  records <- assigned$records
  records$unit_name <- units$unit_name[match(records$unit_id, units$unit_id)]
  admin <- NULL
  if (file.exists(file.path(dir, "admin.gpkg"))) {
    admin <- sf::st_read(file.path(dir, "admin.gpkg"), quiet = TRUE)
    if (!("unit_id" %in% names(admin))) {
      warning("admin.gpkg has no unit_id column; unit outlines will not be drawn.", call. = FALSE)
      admin <- NULL
    }
  }
  list(covs = covs, units = units, candidates = candidates, records = records,
       admin = admin, record_problems = assigned$report)
}

# ---- corrections: one row per marked cell, grouped by action_id -------------

empty_corrections <- function() {
  data.frame(action_id = character(), record_id = character(), cell_id = numeric(),
             lon = numeric(), lat = numeric(), weight = numeric(),
             weight_mode = character(), method = character(), source = character(),
             confidence = character(), user = character(), timestamp = character(),
             stringsAsFactors = FALSE)
}

new_action_id <- function() {
  sprintf("%s-%04d", format(Sys.time(), "%Y%m%dT%H%M%S"), sample.int(9999, 1))
}

read_corrections <- function(dir = "data") {
  f <- file.path(dir, "corrections.csv")
  if (!file.exists(f)) return(empty_corrections())
  cls <- c(record_id = "character", action_id = "character")
  cls <- cls[names(cls) %in% names(utils::read.csv(f, nrows = 0))]
  d <- utils::read.csv(f, stringsAsFactors = FALSE, na.strings = c("", "NA"), colClasses = cls)
  # Files written by the first prototype have one row per record and no action
  # id or weight: each row is its own single-cell action.
  if (is.null(d$action_id)) d$action_id <- paste0("legacy-", seq_len(nrow(d)))
  if (is.null(d$weight)) d$weight <- 1
  if (is.null(d$weight_mode)) d$weight_mode <- "even"
  for (col in setdiff(CORRECTION_COLS, names(d))) d[[col]] <- NA
  # Earlier versions stored the cell centre as lon/lat for cell-level marks. Only
  # an exact point has real coordinates, so blank them; the centre is derived
  # from cell_id wherever it is needed.
  cell_level <- d$method %in% c("cell", "cells")
  d$lon[cell_level] <- NA_real_
  d$lat[cell_level] <- NA_real_
  d[, CORRECTION_COLS, drop = FALSE]
}

#' Append one save (one action): one row per marked cell. Weights are
#' proportions and must sum to 1. Only an exact point (method "point") carries
#' lon/lat; a marked cell is identified by cell_id alone, so a cell centre is never
#' stored as if it were a measured location.
append_correction <- function(rows, dir = "data") {
  stopifnot(all(CORRECTION_COLS %in% names(rows)), nrow(rows) >= 1,
            length(unique(rows$action_id)) == 1,
            all(rows$method %in% c("cell", "cells", "point")),
            all(!is.na(rows$cell_id)),
            all(rows$weight > 0), abs(sum(rows$weight) - 1) < 1e-6)
  is_point <- rows$method == "point"
  if (any(is_point & (is.na(rows$lon) | is.na(rows$lat))))
    stop("an exact point needs lon and lat", call. = FALSE)
  if (any(!is_point & (!is.na(rows$lon) | !is.na(rows$lat))))
    stop("a marked cell must not carry lon/lat; they are derived from cell_id", call. = FALSE)
  rows <- rows[, CORRECTION_COLS, drop = FALSE]
  f <- file.path(dir, "corrections.csv")
  utils::write.table(rows, f, sep = ",", row.names = FALSE, na = "",
                     col.names = !file.exists(f), append = file.exists(f),
                     qmethod = "double")
  invisible(rows)
}

#' Rows of the most recent action for each record. The log is append-only, so
#' "latest" is the last action for that record in file order.
latest_action_rows <- function(corrections) {
  if (nrow(corrections) == 0) return(corrections)
  last_id <- tapply(corrections$action_id, corrections$record_id, function(a) a[length(a)])
  corrections[paste(corrections$record_id, corrections$action_id) %in%
                paste(names(last_id), last_id), , drop = FALSE]
}

# ---- flags: located sites flagged as possibly wrong -------------------------

empty_flags <- function() {
  data.frame(action = character(), record_id = character(), user = character(),
             timestamp = character(), note = character(), stringsAsFactors = FALSE)
}

read_flags <- function(dir = "data") {
  f <- file.path(dir, "flags.csv")
  if (!file.exists(f)) return(empty_flags())
  d <- utils::read.csv(f, stringsAsFactors = FALSE, na.strings = c("", "NA"),
                       colClasses = c(record_id = "character", note = "character"))
  d[, FLAG_COLS, drop = FALSE]
}

#' Append a flag (action = "flag", by an entomologist) or a dismissal
#' (action = "dismiss", by a modeller). Flags never change a record's location.
append_flag <- function(row, dir = "data") {
  stopifnot(all(FLAG_COLS %in% names(row)), row$action %in% c("flag", "dismiss"))
  row <- row[, FLAG_COLS, drop = FALSE]
  f <- file.path(dir, "flags.csv")
  utils::write.table(row, f, sep = ",", row.names = FALSE, na = "",
                     col.names = !file.exists(f), append = file.exists(f),
                     qmethod = "double")
  invisible(row)
}

#' Record ids whose latest flag action is "flag" (not since dismissed).
current_flags <- function(flags) {
  if (nrow(flags) == 0) return(character(0))
  last <- tapply(flags$action, flags$record_id, function(a) a[length(a)])
  names(last)[last == "flag"]
}

# ---- record states ----------------------------------------------------------

#' Work out what is known about every record.
#'
#' status: "located"  - one location inside the mask (a cell or a point), or a
#'                      marked-cell set where one cell holds at least
#'                      `located_share` of the weight (that cell is used);
#'         "narrowed" - several marked cells, none heavy enough;
#'         "coarse"   - only the location ID is known.
#' A coarse record whose unit has exactly one candidate cell counts as located there
#' (`basis` "single candidate cell"): nothing is left for anyone to resolve.
#' Returns list(records = records plus `status`, the effective lon/lat, `exact`
#'                        (TRUE where lon/lat are a real point, FALSE where they are the
#'                        centre of the cell, derived from cell_id), and `basis` (what
#'                        made it located: record coordinates, entomologist, single
#'                        candidate cell),
#'              cells   = marked cells of narrowed records with weights and the
#'                        derived cell centres).
record_states <- function(records, corrections, covs, located_share = 0.8, candidates = NULL) {
  records$status <- "coarse"
  records$exact  <- FALSE
  records$basis  <- NA_character_
  has_xy <- !is.na(records$lon) & !is.na(records$lat)
  if (any(has_xy)) {
    v <- terra::extract(covs[[1]], as.matrix(records[has_xy, c("lon", "lat")]))
    records$status[has_xy][!is.na(v[[1]])] <- "located"
    records$exact[has_xy][!is.na(v[[1]])]  <- TRUE     # coordinates that came with the record
    records$basis[has_xy][!is.na(v[[1]])]  <- "record coordinates"
  }
  cells <- data.frame(record_id = character(), cell_id = numeric(),
                      lon = numeric(), lat = numeric(), weight = numeric(),
                      stringsAsFactors = FALSE)
  last <- latest_action_rows(corrections)
  for (id in unique(last$record_id)) {
    i <- match(id, records$record_id)
    if (is.na(i)) next
    r <- last[last$record_id == id & !is.na(last$cell_id), , drop = FALSE]
    if (nrow(r) == 0) next
    # Cell centres come from the cell ID; only an exact point keeps its own coordinates.
    xy <- terra::xyFromCell(covs[[1]], r$cell_id)
    is_point <- r$method %in% "point" & !is.na(r$lon) & !is.na(r$lat)
    r$lon <- ifelse(is_point, r$lon, xy[, 1])
    r$lat <- ifelse(is_point, r$lat, xy[, 2])
    heavy <- which.max(r$weight)
    unique_heavy <- sum(r$weight == r$weight[heavy]) == 1L
    records$basis[i] <- "entomologist"
    if (unique_heavy && r$weight[heavy] >= located_share - 1e-9) {
      records$status[i] <- "located"
      records$lon[i] <- r$lon[heavy]; records$lat[i] <- r$lat[heavy]
      records$exact[i] <- is_point[heavy]
    } else {
      records$status[i] <- "narrowed"
      cells <- rbind(cells, r[, c("record_id", "cell_id", "lon", "lat", "weight")])
    }
  }
  # A unit whose polygon covers exactly one mask cell has nothing left to resolve: its
  # records are located at that cell. (An entomologist's own mark always takes precedence.)
  if (!is.null(candidates) && nrow(candidates)) {
    n_cand <- table(candidates$unit_id)
    single <- names(n_cand)[n_cand == 1]
    auto <- records$status == "coarse" & records$unit_id %in% single
    if (any(auto)) {
      cell <- candidates$cell_id[match(records$unit_id[auto], candidates$unit_id)]
      xy <- terra::xyFromCell(covs[[1]], cell)
      records$status[auto] <- "located"
      records$lon[auto] <- xy[, 1]; records$lat[auto] <- xy[, 2]
      records$exact[auto] <- FALSE
      records$basis[auto] <- "single candidate cell"
    }
  }
  list(records = records, cells = cells)
}

#' lon/lat of the sites that define the reference set: located records, minus
#' flagged ones if the modeller chooses to exclude them until reviewed.
reference_xy <- function(states, flagged_ids, exclude_flagged = TRUE) {
  r <- states$records
  keep <- r$status == "located"
  if (exclude_flagged) keep <- keep & !(r$record_id %in% flagged_ids)
  r[keep, c("lon", "lat"), drop = FALSE]
}
