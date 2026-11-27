# A check of real inputs before they are used, built to be run locally on data that
# cannot leave the machine and then shared as text.
#
# The report contains counts and summaries only: never a coordinate, a record id, a
# place name or a unit id. Counts below `suppress_below` are shown as "<N" (a common
# rule in data-sharing agreements). It is still the reader's job to check the report
# against their own agreement before sharing it. Set `details = TRUE` to see full error
# messages in your own console, which may contain ids or names: do not share that output.

# ---- a small report object ---------------------------------------------------

new_report <- function(suppress_below = 5, details = FALSE) {
  r <- new.env()
  r$lines <- character(); r$problems <- character(); r$warnings <- character()
  r$suppress <- suppress_below; r$details <- details
  r
}
rp_head    <- function(r, title) r$lines <- c(r$lines, "", paste0("== ", title, " =="))
rp_line    <- function(r, ...) r$lines <- c(r$lines, paste0("  ", ...))
rp_problem <- function(r, ...) { m <- paste0(...); r$problems <- c(r$problems, m); r$lines <- c(r$lines, paste0("  PROBLEM: ", m)) }
rp_warn    <- function(r, ...) { m <- paste0(...); r$warnings <- c(r$warnings, m); r$lines <- c(r$lines, paste0("  warning: ", m)) }
#' A count, with small non-zero counts suppressed.
rp_n <- function(r, n) {
  if (is.na(n)) return("unknown")
  if (r$suppress > 0 && n > 0 && n < r$suppress) paste0("<", r$suppress) else format(n, big.mark = ",", trim = TRUE)
}
rp_pct <- function(r, n, total) {
  if (total == 0) return("")
  if (r$suppress > 0 && n > 0 && n < r$suppress) return("")
  sprintf(" (%.1f%%)", 100 * n / total)
}
#' Run a check; if it fails, say so without the error text (which could hold ids or names).
rp_try <- function(r, what, expr) {
  tryCatch(expr, error = function(e) {
    rp_problem(r, what, " could not be checked",
               if (r$details) paste0(": ", conditionMessage(e)) else " (run with details = TRUE locally to see why)")
  })
}

# ---- the mask and covariates -------------------------------------------------

check_mask <- function(r, covs) {
  rp_head(r, "Mask and covariates")
  rp_line(r, "layers: ", terra::nlyr(covs), " (", paste(names(covs), collapse = ", "), ")")
  rp_line(r, "cells: ", format(terra::ncell(covs), big.mark = ","), "; extent (xmin, xmax, ymin, ymax): ",
          paste(round(as.vector(terra::ext(covs)), 3), collapse = ", "))
  if (!terra::is.lonlat(covs)) {
    rp_problem(r, "the mask is not in longitude/latitude coordinates. The app reads record coordinates, ",
               "draws cells and sizes them as degrees, so a projected mask would give wrong results. ",
               "Reproject the mask and every layer to WGS84 (EPSG:4326).")
  } else {
    d <- terra::crs(covs, describe = TRUE)
    wgs84 <- identical(as.character(d$code), "4326") || grepl("WGS 84", d$name, fixed = TRUE)
    if (wgs84) rp_line(r, "coordinate system: WGS84 longitude/latitude") else
      rp_warn(r, "geographic coordinates, but not recognised as WGS84 (", d$name, "). Record coordinates are assumed to be WGS84.")
  }
  rs <- terra::res(covs)
  mid <- mean(as.vector(terra::ext(covs))[3:4])
  km <- c(rs[1] * 111.32 * cos(mid * pi / 180), rs[2] * 110.57)
  rp_line(r, "cell size: ", signif(rs[1], 4), " x ", signif(rs[2], 4), " degrees (about ",
          signif(km[1], 3), " x ", signif(km[2], 3), " km at the middle of the extent)")
  if (abs(rs[1] - rs[2]) / rs[1] > 1e-6) rp_warn(r, "cells are not square in degrees")

  m <- !is.na(terra::values(covs[[1]], mat = FALSE))
  rp_line(r, "inside the mask (first layer not NA): ", format(sum(m), big.mark = ","), " cells (",
          sprintf("%.1f", 100 * mean(m)), "% of the raster)")
  rp_line(r, "memory for all layers in memory, roughly: ", format(signif(terra::ncell(covs) * terra::nlyr(covs) * 8 / 1e6, 2), big.mark = ","), " MB")

  # alignment of the masks across layers
  na_mat <- is.na(terra::values(covs, mat = TRUE))
  inconsistent <- sum(rowSums(na_mat) %in% seq_len(ncol(na_mat) - 1))
  if (ncol(na_mat) > 1 && inconsistent > 0)
    rp_problem(r, "the layers are not masked identically: ", rp_n(r, inconsistent),
               " cells are NA in some layers but not others. The mask is taken from the first layer, so ",
               "covariate gaps inside it would break scoring.")
  else rp_line(r, "all layers are masked identically")

  s <- terra::global(covs, c("min", "max", "mean", "sd"), na.rm = TRUE)
  for (i in seq_len(terra::nlyr(covs)))
    rp_line(r, "  ", names(covs)[i], ": min ", signif(s$min[i], 4), ", max ", signif(s$max[i], 4),
            ", mean ", signif(s$mean[i], 4), ", sd ", signif(s$sd[i], 4))
  flat <- which(is.na(s$sd) | s$sd == 0)
  if (length(flat)) rp_problem(r, "constant covariate(s): ", paste(names(covs)[flat], collapse = ", "),
                               ". A layer with no variation breaks the distance calculation.")
  if (terra::nlyr(covs) > 1 && sum(m) > 10) {
    smp <- terra::spatSample(covs, min(50000, sum(m)), method = "regular", na.rm = TRUE, as.df = TRUE)
    smp <- stats::na.omit(smp)
    smp <- smp[, vapply(smp, function(x) isTRUE(stats::sd(x) > 0), TRUE), drop = FALSE]
    if (ncol(smp) > 1 && nrow(smp) > 10) {
      cc <- stats::cor(smp); diag(cc) <- 0
      hi <- which(abs(cc) > 0.95, arr.ind = TRUE); hi <- hi[hi[, 1] < hi[, 2], , drop = FALSE]
      if (nrow(hi)) rp_warn(r, "highly correlated covariates (|r| > 0.95): ",
                            paste(sprintf("%s & %s", colnames(cc)[hi[, 1]], colnames(cc)[hi[, 2]]), collapse = "; "),
                            ". Near-duplicate layers make the distance calculation unstable.")
    }
  }
  invisible(m)
}

# ---- other rasters (land cover, dwellings, travel time, ...) -----------------

check_layers <- function(r, layers_dir, covs) {
  rp_head(r, "Other raster layers")
  files <- list.files(layers_dir, pattern = "\\.(tif|tiff|vrt)$", ignore.case = TRUE, full.names = TRUE)
  files <- files[basename(files) != "covariates.tif"]
  if (!length(files)) { rp_line(r, "none found in the layers folder"); return(invisible()) }
  m <- !is.na(terra::values(covs[[1]], mat = FALSE))
  for (f in files) rp_try(r, paste0("layer ", basename(f)), {
    x <- terra::rast(f)
    why <- character()
    if (!terra::same.crs(x, covs)) why <- c(why, "coordinate system differs")
    if (!isTRUE(all.equal(as.vector(terra::ext(x)), as.vector(terra::ext(covs))))) why <- c(why, "extent differs")
    if (!isTRUE(all.equal(terra::res(x), terra::res(covs)))) why <- c(why, "cell size differs")
    if (length(why)) {
      rp_problem(r, basename(f), " is not aligned to the mask (", paste(why, collapse = ", "), "). ",
                 "Resample or reproject it to the mask first.")
    } else {
      v <- terra::values(x[[1]], mat = FALSE)
      na_in <- sum(is.na(v) & m)
      u <- unique(stats::na.omit(v[seq_len(min(length(v), 500000))]))
      kind <- if (length(u) <= 50 && all(u == round(u))) paste0("looks categorical, ", length(u), " distinct values") else
        paste0("continuous, range ", signif(min(v, na.rm = TRUE), 4), " to ", signif(max(v, na.rm = TRUE), 4))
      rp_line(r, basename(f), ": aligned; ", kind, "; NA inside the mask: ", rp_n(r, na_in), " cells")
      if (na_in > 0) rp_warn(r, basename(f), " has gaps inside the mask; decide how they should be filled before using it as a prior.")
    }
  })
}

# ---- polygon layers ------------------------------------------------------------

#' `layers`: list of list(file, id_col, name_col, prefix), as in data-raw/build_candidates.R.
#' Returns the combined conversion (candidates, units, report) or NULL.
check_polygons <- function(r, layers, covs) {
  rp_head(r, "Polygon layers (approximate locations)")
  built <- list()
  for (i in seq_along(layers)) {
    l <- layers[[i]]; tag <- paste0("layer ", i, " (", basename(l$file), ")")
    if (!file.exists(l$file)) { rp_problem(r, tag, ": file not found"); next }
    rp_try(r, tag, {
      p <- sf::st_read(l$file, quiet = TRUE)
      rp_line(r, tag, ": ", rp_n(r, nrow(p)), " polygons")
      ok <- TRUE
      if (crs_missing(p)) { rp_problem(r, tag, " has no coordinate system"); ok <- FALSE }
      if (!(l$id_col %in% names(p))) { rp_problem(r, tag, " has no id column '", l$id_col, "'"); ok <- FALSE }
      if (!(l$name_col %in% names(p))) { rp_problem(r, tag, " has no name column '", l$name_col, "'"); ok <- FALSE }
      if (ok) {
        id <- as.character(p[[l$id_col]])
        if (anyNA(id) || any(!nzchar(id))) { rp_problem(r, tag, ": ", rp_n(r, sum(is.na(id) | !nzchar(id))), " polygons have no id"); ok <- FALSE }
        if (anyDuplicated(id)) { rp_problem(r, tag, ": ", rp_n(r, sum(duplicated(id))), " duplicate ids"); ok <- FALSE }
        rp_line(r, "  empty geometries: ", rp_n(r, sum(sf::st_is_empty(p))), "; invalid: ", rp_n(r, sum(!sf::st_is_valid(p))))
      }
      if (ok) built[[length(built) + 1]] <- build_candidates(p, covs, l$id_col, l$name_col, if (is.null(l$prefix)) "" else l$prefix)
    })
  }
  if (!length(built)) return(invisible(NULL))
  out <- tryCatch(do.call(combine_candidates, built), error = function(e) {
    rp_problem(r, "unit ids repeat across layers; give each layer its own prefix"); NULL })
  if (is.null(out)) return(invisible(NULL))
  u <- out$units
  rp_line(r, "units: ", rp_n(r, nrow(u)), "; candidate cells per unit: min ", min(u$n_cells), ", median ",
          stats::median(u$n_cells), ", max ", max(u$n_cells))
  rp_line(r, "units with no cells in the mask: ", rp_n(r, sum(u$n_cells == 0)), "; covering a single cell (their records will be located automatically): ",
          rp_n(r, sum(u$n_cells == 1)), "; mostly outside the mask: ", rp_n(r, sum(u$n_cells > 0 & !is.na(u$mask_share) & u$mask_share < 0.5)))
  if (any(u$n_cells == 0)) rp_warn(r, "units with no cells cannot be ranked; their records will be listed as a problem in the app")
  invisible(out)
}

# ---- records ---------------------------------------------------------------------

check_records <- function(r, records, covs, units = NULL) {
  rp_head(r, "Records")
  n <- nrow(records)
  rp_line(r, "records: ", rp_n(r, n), "; columns: ", paste(names(records), collapse = ", "))
  if (!("record_id" %in% names(records))) { rp_problem(r, "no record_id column"); return(invisible()) }
  if (anyDuplicated(records$record_id)) rp_problem(r, rp_n(r, sum(duplicated(records$record_id))), " duplicate record ids")

  # coordinates
  lon <- suppressWarnings(as.numeric(records$lon)); lat <- suppressWarnings(as.numeric(records$lat))
  has <- !is.na(lon) & !is.na(lat)
  rp_line(r, "with coordinates: ", rp_n(r, sum(has)), rp_pct(r, sum(has), n), "; without (coarse): ", rp_n(r, sum(!has)))
  bad_range <- has & (abs(lon) > 180 | abs(lat) > 90)
  if (any(bad_range)) rp_problem(r, rp_n(r, sum(bad_range)), " records have coordinates outside the valid longitude/latitude range")
  in_ext <- has & !bad_range & lon >= terra::xmin(covs) & lon <= terra::xmax(covs) & lat >= terra::ymin(covs) & lat <= terra::ymax(covs)
  cell <- rep(NA_real_, n)
  if (any(in_ext)) cell[in_ext] <- terra::cellFromXY(covs[[1]], cbind(lon[in_ext], lat[in_ext]))
  val <- rep(NA_real_, n); if (any(in_ext)) val[in_ext] <- terra::values(covs[[1]], mat = FALSE)[cell[in_ext]]
  inside <- in_ext & !is.na(val)
  rp_line(r, "coordinates inside the mask: ", rp_n(r, sum(inside)), rp_pct(r, sum(inside), sum(has)))
  rp_line(r, "  outside the raster extent: ", rp_n(r, sum(has & !in_ext)), "; inside the extent but outside the mask: ", rp_n(r, sum(in_ext & is.na(val))))
  out <- has & !inside
  if (any(out)) {
    rp_warn(r, rp_n(r, sum(out)), " records with coordinates fall outside the mask and will count as coarse")
    sw <- out & abs(lat) <= 180 & abs(lon) <= 90 &
      lat >= terra::xmin(covs) & lat <= terra::xmax(covs) & lon >= terra::ymin(covs) & lon <= terra::ymax(covs)
    sw_in <- rep(FALSE, n)
    if (any(sw)) sw_in[sw] <- !is.na(terra::extract(covs[[1]], cbind(lat[sw], lon[sw]))[[1]])
    if (any(sw_in)) rp_warn(r, rp_n(r, sum(sw_in)), " of them would be inside the mask if longitude and latitude were swapped")
  }

  # units
  if (!is.null(units)) {
    if ("unit_id" %in% names(records)) {
      has_u <- !is.na(records$unit_id)
      rp_line(r, "with a unit_id: ", rp_n(r, sum(has_u)), "; unit_id not among the units: ", rp_n(r, sum(has_u & !(records$unit_id %in% units$unit_id))))
    } else if ("admin_name" %in% names(records)) {
      am <- assign_unit_ids(records, units)
      named <- sum(!is.na(records$admin_name))
      miss <- am$report
      rp_line(r, "named by admin_name: ", rp_n(r, named), "; matched to exactly one unit: ", rp_n(r, sum(!is.na(am$records$unit_id))))
      rp_line(r, "  matching no unit: ", rp_n(r, sum(miss$problem == "name matches no unit")),
              "; matching more than one unit: ", rp_n(r, sum(miss$problem == "name matches more than one unit")))
      if (nrow(miss)) rp_warn(r, rp_n(r, nrow(miss)), " records cannot be matched to a unit by name; add a unit_id or a unit_layer column")
      records$unit_id <- am$records$unit_id
    } else rp_problem(r, "records have neither a unit_id nor an admin_name column")
    if ("unit_id" %in% names(records)) {
      located_with_unit <- sum(inside & !is.na(records$unit_id))
      rp_line(r, "located records that also have a unit (needed to check a prior on probable location): ", rp_n(r, located_with_unit))
    }
  }

  # dates
  if ("date" %in% names(records)) {
    d <- trimws(as.character(records$date)); d[is.na(d)] <- ""
    full <- grepl("^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", d); mon <- grepl("^[0-9]{4}-[0-9]{1,2}$", d); yr <- grepl("^[0-9]{4}$", d)
    none <- d == ""; other <- !(full | mon | yr | none)
    rp_line(r, "dates: full ", rp_n(r, sum(full)), ", month ", rp_n(r, sum(mon)), ", year ", rp_n(r, sum(yr)),
            ", missing ", rp_n(r, sum(none)), ", other format ", rp_n(r, sum(other)))
    if (any(other)) {
      pat <- gsub("[0-9]", "9", gsub("[A-Za-z]", "A", d[other]))
      top <- utils::head(sort(table(pat), decreasing = TRUE), 3)
      rp_warn(r, "dates in another format will be treated as undated. Patterns seen (9 = digit, A = letter): ",
              paste(names(top), collapse = ", "))
    }
    p <- parse_partial_date(d)
    if (any(!is.na(p$year))) {
      rp_line(r, "years covered: ", min(p$year, na.rm = TRUE), " to ", max(p$year, na.rm = TRUE))
      if (min(p$year, na.rm = TRUE) < 1900 || max(p$year, na.rm = TRUE) > as.integer(format(Sys.Date(), "%Y")) + 1)
        rp_warn(r, "some years look implausible")
    }
  } else rp_warn(r, "no date column: the entomologist view will show no date span")
  invisible(inside)
}

# ---- how long scoring takes ------------------------------------------------------

check_timing <- function(r, covs, records, cand = NULL) {
  rp_head(r, "Scoring time")
  lon <- suppressWarnings(as.numeric(records$lon)); lat <- suppressWarnings(as.numeric(records$lat))
  has <- !is.na(lon) & !is.na(lat)
  if (!any(has)) { rp_line(r, "no records with coordinates, so no reference set to time"); return(invisible()) }
  ref <- data.frame(lon = lon[has], lat = lat[has])
  ref <- ref[!is.na(terra::extract(covs[[1]], as.matrix(ref))[[1]]), , drop = FALSE]
  rp_line(r, "reference sites inside the mask: ", rp_n(r, nrow(ref)), "; cells scored: ", format(sum(!is.na(terra::values(covs[[1]], mat = FALSE))), big.mark = ","))
  if (nrow(ref) < 8) { rp_warn(r, "fewer than 8 located records inside the mask: the app cannot score yet"); return(invisible()) }
  t <- system.time(sc <- score_cells(covs, ref))[["elapsed"]]
  rp_line(r, "MESS and SHAPE for the whole mask: ", sprintf("%.1f", t), " seconds. The app does this again after every save.")
  if (t > 10) rp_warn(r, "scoring is slow enough to make saving feel sluggish; see the scoring-speed note in the issues")
  if (!is.null(cand) && !is.null(records$unit_id)) {
    t2 <- system.time({
      st <- record_states(records[, intersect(c("record_id", "unit_id", "lon", "lat", "date"), names(records))], empty_corrections(), covs, candidates = cand$candidates)
      rs <- record_summary(cand$candidates, sc, st)
    })[["elapsed"]]
    rp_line(r, "ranking ", rp_n(r, nrow(rs)), " unlocated records over their candidate cells: ", sprintf("%.1f", t2), " seconds")
  }
  invisible(t)
}

# ---- the whole check ---------------------------------------------------------------

#' Check real inputs. `dir` holds covariates.tif and records.csv. `layers_dir` (optional)
#' holds other rasters. `polygons` is a list of list(file, id_col, name_col, prefix).
#' Returns the report (invisibly); the text is printed.
check_inputs <- function(dir, layers_dir = NULL, polygons = NULL, suppress_below = 5,
                         time_scoring = TRUE, details = FALSE) {
  r <- new_report(suppress_below, details)
  r$lines <- c("DATA CHECK",
               paste0("This report holds counts and summaries only: no coordinates, record ids, place names or unit ids. ",
                      "Counts below ", suppress_below, " are shown as '<", suppress_below, "'. Read it against your data-sharing agreement before sharing it."),
               paste0("R ", getRversion(), "; terra ", utils::packageVersion("terra"), "; sf ", utils::packageVersion("sf"),
                      if (requireNamespace("flexsdm", quietly = TRUE)) paste0("; flexsdm ", utils::packageVersion("flexsdm")) else "; flexsdm not installed"))
  covs <- NULL; records <- NULL; cand <- NULL
  f_cov <- file.path(dir, "covariates.tif"); f_rec <- file.path(dir, "records.csv")
  if (!file.exists(f_cov)) rp_problem(r, "covariates.tif not found in the data folder") else
    rp_try(r, "the mask", { covs <- terra::rast(f_cov); check_mask(r, covs) })
  if (!is.null(layers_dir) && !is.null(covs)) rp_try(r, "the other layers", check_layers(r, layers_dir, covs))
  if (!is.null(polygons) && !is.null(covs)) rp_try(r, "the polygon layers", cand <- check_polygons(r, polygons, covs))
  if (!file.exists(f_rec)) rp_problem(r, "records.csv not found in the data folder") else if (!is.null(covs)) {
    rp_try(r, "the records", {
      records <- utils::read.csv(f_rec, stringsAsFactors = FALSE, na.strings = c("", "NA"),
                                 colClasses = c(record_id = "character", unit_id = "character", admin_name = "character", date = "character")[
                                   c("record_id", "unit_id", "admin_name", "date") %in% names(utils::read.csv(f_rec, nrows = 0))])
      check_records(r, records, covs, if (!is.null(cand)) cand$units)
    })
  }
  if (time_scoring && !is.null(covs) && !is.null(records))
    rp_try(r, "scoring time", {
      recs <- records
      if (!is.null(cand) && is.null(recs$unit_id) && "admin_name" %in% names(recs)) recs$unit_id <- assign_unit_ids(recs, cand$units)$records$unit_id
      check_timing(r, covs, recs, cand)
    })
  rp_head(r, "Summary")
  rp_line(r, length(r$problems), " problem(s), ", length(r$warnings), " warning(s)")
  for (p in r$problems) rp_line(r, "PROBLEM: ", p)
  for (w in r$warnings) rp_line(r, "warning: ", w)
  cat(r$lines, sep = "\n"); cat("\n")
  invisible(r)
}
