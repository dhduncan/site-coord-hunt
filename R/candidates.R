# Candidate cells: turning approximate locations (polygons) into model cells.
#
# A coarse record is attributed to a unit (a parish, a facility catchment, ...). The
# model works on a raster mask, so every unit needs a list of the mask cells it could
# be in, and how likely each is. That list is the candidate table:
#
#   unit_id, cell_id, weight
#
# - A cell can appear under several units (a facility catchment inside a parish), so
#   the table is many-to-many.
# - weight is the share of the polygon's area that falls in the cell, i.e. the chance
#   the true location is there if it is equally likely anywhere in the polygon. It is
#   normalised over mask cells only and sums to 1 for each unit.
# - cell_id is the cell number in the mask raster (fixed: every resource is aligned
#   to it).

#' TRUE if a layer has no usable coordinate system. A missing one is NA, but formats such as
#' GeoPackage store "Undefined Cartesian SRS" instead, which reads back as a defined system.
crs_missing <- function(x) {
  cr <- sf::st_crs(x)
  is.na(cr) || isTRUE(grepl("Undefined Cartesian|unknown unit", cr$input))
}

#' Convert one polygon layer to candidate cells.
#'
#' `polygons`: sf layer. `mask`: SpatRaster; NA cells are outside the mask.
#' `id_col`, `name_col`: columns holding the unit's id and its display name.
#' `prefix`: put in front of every id (e.g. "parish") so ids stay unique across layers.
#' Returns list(candidates, units, report). `report` lists units that could not be
#' converted properly (no cells in the mask, mostly outside it, invalid geometry);
#' nothing is dropped without being listed.
build_candidates <- function(polygons, mask, id_col, name_col = id_col, prefix = "") {
  stopifnot(inherits(polygons, "sf"), id_col %in% names(polygons), name_col %in% names(polygons))
  if (crs_missing(polygons))
    stop("The polygon layer has no coordinate system; set it before converting.", call. = FALSE)
  ids <- as.character(polygons[[id_col]])
  if (anyNA(ids) || any(!nzchar(ids))) stop("Every polygon needs an id.", call. = FALSE)
  unit_id <- if (nzchar(prefix)) paste0(prefix, ":", ids) else ids
  if (anyDuplicated(unit_id))
    stop("Duplicate unit ids in the polygon layer: ",
         toString(unique(unit_id[duplicated(unit_id)])), call. = FALSE)

  report <- data.frame(unit_id = character(), problem = character(), detail = character(),
                       stringsAsFactors = FALSE)
  add_problem <- function(id, problem, detail = "")
    report[nrow(report) + 1, ] <<- c(id, problem, detail)

  polygons <- sf::st_transform(polygons, terra::crs(mask))
  bad <- !sf::st_is_valid(polygons)
  if (any(bad)) {
    polygons[bad, ] <- sf::st_make_valid(polygons[bad, ])
    for (u in unit_id[bad]) add_problem(u, "invalid geometry", "repaired with st_make_valid()")
  }
  empty <- sf::st_is_empty(polygons)

  m1 <- mask[[1]]
  keep <- which(!empty)                 # terra cannot convert an empty geometry
  e <- terra::extract(m1, terra::vect(polygons[keep, ]), cells = TRUE, exact = TRUE)
  e$ID <- keep[e$ID]                    # back to row numbers in `polygons`
  e <- e[!is.na(e$cell), , drop = FALSE]
  area_all <- terra::cellSize(m1, unit = "m")
  e$area <- terra::values(area_all, mat = FALSE)[e$cell] * e$fraction  # polygon area in the cell
  in_mask <- !is.na(e[[names(m1)]])
  # Same area measure as terra::cellSize(), so the share inside the mask is consistent.
  poly_area <- numeric(nrow(polygons))
  poly_area[keep] <- terra::expanse(terra::vect(polygons[keep, ]), unit = "m")
  inside <- tapply(e$area[in_mask], factor(e$ID[in_mask], levels = seq_len(nrow(polygons))), sum)
  inside[is.na(inside)] <- 0
  e <- e[in_mask, , drop = FALSE]
  e$unit_id <- unit_id[e$ID]
  tot <- ave(e$area, e$unit_id, FUN = sum)
  cand <- data.frame(unit_id = e$unit_id, cell_id = e$cell, weight = e$area / tot,
                     stringsAsFactors = FALSE)
  cand <- cand[order(cand$unit_id, cand$cell_id), , drop = FALSE]
  rownames(cand) <- NULL

  n_cells <- as.integer(table(factor(cand$unit_id, levels = unit_id)))
  share   <- ifelse(poly_area > 0, pmin(1, as.numeric(inside) / poly_area), NA_real_)
  for (k in which(n_cells == 0))
    add_problem(unit_id[k], "no cells in the mask",
                if (empty[k]) "the geometry is empty" else "the polygon lies outside the mask")
  for (k in which(n_cells > 0 & !is.na(share) & share < 0.5))
    add_problem(unit_id[k], "mostly outside the mask",
                sprintf("only %.0f%% of its area is inside; candidates cover that part only",
                        100 * share[k]))

  units <- data.frame(unit_id = unit_id, unit_name = as.character(polygons[[name_col]]),
                      layer = if (nzchar(prefix)) prefix else NA_character_,
                      n_cells = n_cells, mask_share = round(share, 4),
                      stringsAsFactors = FALSE)
  list(candidates = cand, units = units, report = report)
}

#' Combine the results for several layers. Ids must be unique across layers.
combine_candidates <- function(...) {
  parts <- list(...)
  units <- do.call(rbind, lapply(parts, `[[`, "units"))
  if (anyDuplicated(units$unit_id))
    stop("Unit ids repeat across layers: ", toString(unique(units$unit_id[duplicated(units$unit_id)])),
         ". Give each layer its own prefix.", call. = FALSE)
  list(candidates = do.call(rbind, lapply(parts, `[[`, "candidates")),
       units = units,
       report = do.call(rbind, lapply(parts, `[[`, "report")))
}

#' Give every record a unit_id.
#'
#' If the records already have `unit_id`, it is checked against `units`. Otherwise the
#' record's `admin_name` (the name from the source data) is matched to a unit name, within
#' `unit_layer` if the records carry one. A name that matches no unit, or more than one,
#' is reported and left NA: it is never guessed.
#' Returns list(records, report).
assign_unit_ids <- function(records, units) {
  report <- data.frame(record_id = character(), problem = character(), detail = character(),
                       stringsAsFactors = FALSE)
  if (!is.null(records$unit_id)) {
    unknown <- !is.na(records$unit_id) & !(records$unit_id %in% units$unit_id)
    if (any(unknown))
      report <- rbind(report, data.frame(record_id = records$record_id[unknown],
                                         problem = "unknown unit_id",
                                         detail = records$unit_id[unknown], stringsAsFactors = FALSE))
    return(list(records = records, report = report))
  }
  nm <- records$admin_name
  records$unit_id <- NA_character_
  for (i in which(!is.na(nm))) {
    pool <- units[units$unit_name == nm[i], , drop = FALSE]
    if (!is.null(records$unit_layer) && !is.na(records$unit_layer[i]))
      pool <- pool[pool$layer %in% records$unit_layer[i], , drop = FALSE]
    if (nrow(pool) == 1) {
      records$unit_id[i] <- pool$unit_id
    } else {
      report[nrow(report) + 1, ] <- c(records$record_id[i],
        if (nrow(pool) == 0) "name matches no unit" else "name matches more than one unit",
        if (nrow(pool) == 0) nm[i] else paste(nm[i], "->", toString(pool$unit_id)))
    }
  }
  list(records = records, report = report)
}

#' Units that records point to but the app cannot use: ids missing from `units`, and
#' units with no candidate cells. Shown to the modeller, never silently dropped.
unit_problems <- function(records, units, candidates) {
  used <- unique(records$unit_id[!is.na(records$unit_id)])
  list(unknown = setdiff(used, units$unit_id),
       no_cells = setdiff(intersect(used, units$unit_id), unique(candidates$unit_id)))
}

#' Refuse a candidate table that does not fit the mask or the units: it would give
#' silently wrong results. Cell IDs are fixed (the mask never changes), so a mismatch
#' means the table was built for a different raster.
check_candidates <- function(candidates, units, covs) {
  need <- c("unit_id", "cell_id", "weight")
  if (!all(need %in% names(candidates)))
    stop("candidates.csv needs columns: ", toString(need), call. = FALSE)
  if (anyDuplicated(units$unit_id)) stop("units.csv has repeated unit_id values.", call. = FALSE)
  if (nrow(candidates) == 0) return(invisible(TRUE))
  if (any(candidates$cell_id < 1 | candidates$cell_id > terra::ncell(covs) | is.na(candidates$cell_id)))
    stop("candidates.csv has cell IDs outside the mask raster; was it built for a different raster?", call. = FALSE)
  if (anyNA(terra::values(covs[[1]], mat = FALSE)[candidates$cell_id]))
    stop("candidates.csv has cells that are outside the mask; was it built for a different raster?", call. = FALSE)
  if (anyDuplicated(candidates[, c("unit_id", "cell_id")]))
    stop("candidates.csv lists the same cell twice for one unit.", call. = FALSE)
  if (any(candidates$weight <= 0)) stop("candidates.csv has non-positive weights.", call. = FALSE)
  sums <- tapply(candidates$weight, candidates$unit_id, sum)
  if (any(abs(sums - 1) > 1e-6))
    stop("Weights must sum to 1 for each unit; not so for: ", toString(names(sums)[abs(sums - 1) > 1e-6]), call. = FALSE)
  extra <- setdiff(candidates$unit_id, units$unit_id)
  if (length(extra)) stop("candidates.csv has units missing from units.csv: ", toString(extra), call. = FALSE)
  invisible(TRUE)
}

#' Convert polygon layers and write the files the app reads.
#'
#' `layers`: list of list(polygons = <sf>, id_col, name_col, prefix). `mask`: SpatRaster.
#' Writes candidates.csv, units.csv and admin.gpkg (unit outlines for the map, with a
#' unit_id column) into `dir`. Returns the report of units that could not be converted
#' properly, which is also printed.
write_candidate_files <- function(layers, mask, dir = "data") {
  built <- lapply(layers, function(l)
    build_candidates(l$polygons, mask, l$id_col, l$name_col, l$prefix))
  out <- do.call(combine_candidates, built)
  outlines <- do.call(rbind, lapply(seq_along(layers), function(i) {
    pol <- sf::st_transform(layers[[i]]$polygons, terra::crs(mask))
    sf::st_sf(unit_id = built[[i]]$units$unit_id, unit_name = built[[i]]$units$unit_name,
              layer = built[[i]]$units$layer, geometry = sf::st_geometry(pol))
  }))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(out$candidates, file.path(dir, "candidates.csv"), row.names = FALSE)
  utils::write.csv(out$units, file.path(dir, "units.csv"), row.names = FALSE, na = "")
  sf::st_write(outlines, file.path(dir, "admin.gpkg"), delete_dsn = TRUE, quiet = TRUE)
  message(sprintf("Wrote %d candidate rows for %d units to %s/", nrow(out$candidates), nrow(out$units), dir))
  if (nrow(out$report)) {
    message("Units that need attention:")
    print(out$report, row.names = FALSE)
  }
  invisible(out$report)
}
