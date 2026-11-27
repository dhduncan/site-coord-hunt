# Run from repo root: Rscript tests/test_candidates.R
suppressPackageStartupMessages({library(terra); library(sf)})
source("R/candidates.R")
chk <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)
fails <- function(expr) inherits(tryCatch(expr, error = function(e) e), "error")
near <- function(a, b, tol = 1e-4) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol))

# A planar 4 x 4 raster of 1 m cells, so areas can be worked out by hand.
# Cell numbers run along rows from the top left: row 1 = cells 1-4 (y 3-4), row 4 = cells 13-16 (y 0-1).
mk_mask <- function(na_cells = integer()) {
  r <- rast(nrows = 4, ncols = 4, xmin = 0, xmax = 4, ymin = 0, ymax = 4, crs = "EPSG:32636")
  values(r) <- 1; r[na_cells] <- NA; r
}
sq <- function(x0, x1, y0, y1) st_polygon(list(rbind(c(x0,y0), c(x1,y0), c(x1,y1), c(x0,y1), c(x0,y0))))
layer <- function(ids, geoms, crs = 32636) st_sf(id = ids, nm = paste("Unit", ids), geometry = st_sfc(geoms, crs = crs))
W <- function(b, u) setNames(b$candidates$weight[b$candidates$unit_id == u], b$candidates$cell_id[b$candidates$unit_id == u])

# ---- weights are the share of the polygon's area in each cell ------------------
mask <- mk_mask()
b <- build_candidates(layer(c("tiny", "straddle", "wide"),
        list(sq(0.2, 0.8, 0.2, 0.8),    # inside one cell, no cell centre inside it
             sq(0.5, 1.5, 0.5, 1.5),    # a quarter in each of four cells
             sq(0.5, 2.5, 0, 1))),      # half, whole, half along the bottom row
      mask, "id", "nm", prefix = "x")
chk(near(W(b, "x:tiny"), 1) && names(W(b, "x:tiny")) == "13", "a polygon inside one cell has that cell at weight 1")
chk(near(W(b, "x:straddle"), rep(0.25, 4)) && setequal(names(W(b, "x:straddle")), c("9", "10", "13", "14")),
    "a polygon straddling four cells gives each a quarter")
chk(near(W(b, "x:wide"), c(0.25, 0.5, 0.25)) && identical(names(W(b, "x:wide")), c("13", "14", "15")),
    "weights follow the area share: half, whole, half")
chk(all(abs(tapply(b$candidates$weight, b$candidates$unit_id, sum) - 1) < 1e-9), "weights sum to 1 for every unit")
chk(setequal(b$units$unit_id, c("x:tiny", "x:straddle", "x:wide")) && identical(b$units$n_cells[b$units$unit_id == "x:wide"], 3L),
    "units table: ids are prefixed and cells counted")
chk(nrow(b$report) == 0, "a clean conversion reports nothing")
chk(all(abs(b$units$mask_share - 1) < 1e-6), "a polygon fully inside the mask has share 1")

# ---- cells can belong to several units (overlap) -------------------------------
ov <- build_candidates(layer(c("a", "b"), list(sq(0, 2, 0, 2), sq(1, 3, 1, 3))), mask, "id", "nm")
chk(length(intersect(ov$candidates$cell_id[ov$candidates$unit_id == "a"], ov$candidates$cell_id[ov$candidates$unit_id == "b"])) > 0,
    "overlapping units share cells")
chk(all(abs(tapply(ov$candidates$weight, ov$candidates$unit_id, sum) - 1) < 1e-9), "each overlapping unit still sums to 1")

# ---- the mask: weights are renormalised over cells inside it --------------------
holey <- mk_mask(na_cells = c(13, 14))           # bottom row, two left cells outside the mask
h <- build_candidates(layer(c("half", "gone", "away"),
       list(sq(0, 4, 0, 1),                      # the whole bottom row: half of it is outside the mask
            sq(0, 2, 0, 1),                      # entirely over the masked-out cells
            sq(10, 11, 10, 11))),                # entirely outside the raster
     holey, "id", "nm")
chk(near(W(h, "half"), c(0.5, 0.5)) && setequal(names(W(h, "half")), c("15", "16")), "weights renormalise over mask cells")
chk(abs(h$units$mask_share[h$units$unit_id == "half"] - 0.5) < 0.01, "mask_share says how much of the polygon is inside")
chk(h$units$n_cells[h$units$unit_id == "gone"] == 0 && h$units$n_cells[h$units$unit_id == "away"] == 0, "units with no cells are kept in the units table")
chk(!any(h$candidates$unit_id %in% c("gone", "away")), "...and have no candidate rows")
chk(all(c("gone", "away") %in% h$report$unit_id) && all(h$report$problem[h$report$unit_id %in% c("gone", "away")] == "no cells in the mask"),
    "units with no cells in the mask are reported, not dropped silently")
m <- build_candidates(layer("mostly", list(sq(0, 4, 0, 1))), mk_mask(na_cells = c(13, 14, 15)), "id", "nm")
chk(m$report$problem == "mostly outside the mask", "a unit mostly outside the mask is reported")

# ---- awkward input -------------------------------------------------------------
e <- build_candidates(layer(c("ok", "void"), list(sq(0, 1, 0, 1), st_polygon())), mask, "id", "nm")
chk(e$units$n_cells[e$units$unit_id == "void"] == 0 && "void" %in% e$report$unit_id, "an empty geometry is reported, not a crash")
bow <- st_polygon(list(rbind(c(0,0), c(2,2), c(2,0), c(0,2), c(0,0))))          # self-intersecting bowtie
bw <- build_candidates(layer("bow", list(bow)), mask, "id", "nm")
chk("invalid geometry" %in% bw$report$problem && bw$units$n_cells > 0, "an invalid geometry is repaired and reported")
chk(fails(build_candidates(layer(c("d", "d"), list(sq(0,1,0,1), sq(1,2,1,2))), mask, "id", "nm")), "duplicate ids are refused")
chk(fails(build_candidates(layer(c("a", NA), list(sq(0,1,0,1), sq(1,2,1,2))), mask, "id", "nm")), "a missing id is refused")
nocrs <- layer("a", list(sq(0,1,0,1))); st_crs(nocrs) <- NA
chk(fails(build_candidates(nocrs, mask, "id", "nm")), "a polygon layer with no CRS is refused")
# GeoPackage stores a missing CRS as "Undefined Cartesian", which must be refused too
gp <- tempfile(fileext = ".gpkg"); suppressWarnings(st_write(nocrs, gp, quiet = TRUE))
chk(fails(build_candidates(st_read(gp, quiet = TRUE), mask, "id", "nm")), "a layer stored with an undefined CRS is refused")
# a layer in another coordinate system is transformed to the mask's
ll <- st_transform(layer(c("a", "b"), list(sq(0, 2, 0, 2), sq(1, 3, 1, 3))), 4326)
tr <- build_candidates(ll, mask, "id", "nm")
chk(setequal(tr$candidates$cell_id[tr$candidates$unit_id == "a"], ov$candidates$cell_id[ov$candidates$unit_id == "a"]),
    "a layer in another CRS gives the same cells")

# ---- combining layers: ids must stay unique ------------------------------------
l1 <- build_candidates(layer("1", list(sq(0,1,0,1))), mask, "id", "nm", prefix = "parish")
l2 <- build_candidates(layer("1", list(sq(2,3,2,3))), mask, "id", "nm", prefix = "facility")
cmb <- combine_candidates(l1, l2)
chk(setequal(cmb$units$unit_id, c("parish:1", "facility:1")), "prefixes keep ids unique across layers")
chk(fails(combine_candidates(build_candidates(layer("1", list(sq(0,1,0,1))), mask, "id", "nm"),
                             build_candidates(layer("1", list(sq(2,3,2,3))), mask, "id", "nm"))),
    "the same id in two layers is refused")

# ---- matching names to units ---------------------------------------------------
units <- data.frame(unit_id = c("parish:1", "parish:2", "facility:7", "parish:9"),
                    unit_name = c("Central", "Kasese", "Central", "Kasese"),
                    layer = c("parish", "parish", "facility", "parish"), stringsAsFactors = FALSE)
rec <- data.frame(record_id = paste0("r", 1:5), admin_name = c("Kasese", "Central", "Nowhere", NA, "Central"),
                  unit_layer = c(NA, NA, NA, NA, "facility"), stringsAsFactors = FALSE)
am <- assign_unit_ids(rec, units)
chk(identical(am$records$unit_id, c(NA, NA, NA, NA, "facility:7")), "a name is matched only when exactly one unit fits (layer given for r5)")
chk(setequal(am$report$record_id, c("r1", "r2", "r3")) &&
    all(am$report$problem[am$report$record_id %in% c("r1", "r2")] == "name matches more than one unit") &&
    am$report$problem[am$report$record_id == "r3"] == "name matches no unit",
    "ambiguous and unmatched names are reported, never guessed")
ok <- assign_unit_ids(data.frame(record_id = c("a", "b"), unit_id = c("parish:1", "parish:99"), stringsAsFactors = FALSE), units)
chk(nrow(ok$report) == 1 && ok$report$record_id == "b" && ok$report$problem == "unknown unit_id", "an unknown unit_id is reported")
chk(identical(ok$records$unit_id, c("parish:1", "parish:99")), "a unit_id is never changed")

# ---- the candidate table must fit the mask -------------------------------------
good <- b$candidates; gu <- b$units
chk(isTRUE(check_candidates(good, gu, mask)), "a table built for this mask passes")
chk(fails(check_candidates(transform(good, cell_id = cell_id + 100), gu, mask)), "cell ids beyond the raster are refused")
chk(fails(check_candidates(good, gu, holey)), "cells outside the mask are refused (wrong raster)")
chk(fails(check_candidates(transform(good, weight = weight * 0.9), gu, mask)), "weights that do not sum to 1 are refused")
chk(fails(check_candidates(rbind(good, good[1, ]), gu, mask)), "the same cell twice for a unit is refused")
chk(fails(check_candidates(good, gu[gu$unit_id != "x:wide", ], mask)), "units missing from units.csv are refused")
chk(fails(check_candidates(good[, c("unit_id", "cell_id")], gu, mask)), "a missing column is refused")

# ---- writing and reading back --------------------------------------------------
dir <- tempfile(); dir.create(dir)
rep <- write_candidate_files(list(list(polygons = layer(c("tiny", "wide"), list(sq(0.2,0.8,0.2,0.8), sq(0.5,2.5,0,1))),
                                       id_col = "id", name_col = "nm", prefix = "p")), mask, dir)
chk(all(file.exists(file.path(dir, c("candidates.csv", "units.csv", "admin.gpkg")))), "all three files are written")
back <- read.csv(file.path(dir, "candidates.csv"), colClasses = c(unit_id = "character"))
chk(isTRUE(check_candidates(back, read.csv(file.path(dir, "units.csv"), colClasses = c(unit_id = "character")), mask)), "written files pass the check")
out <- st_read(file.path(dir, "admin.gpkg"), quiet = TRUE)
chk(setequal(out$unit_id, c("p:tiny", "p:wide")) && all(c("unit_name", "layer") %in% names(out)), "outlines carry unit_id for the map")

cat("All candidate tests passed\n")
