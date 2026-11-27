# Run from repo root: Rscript tests/test_check_inputs.R
suppressPackageStartupMessages({library(terra); library(sf)})
for (f in c("dates", "scoring", "candidates", "data_io", "check_inputs")) source(file.path("R", paste0(f, ".R")))
chk <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)

# ---- helpers ------------------------------------------------------------------
# A data folder copied from the simulated data, with optional replacements.
mkdir <- function(records = NULL, covs = NULL) {
  d <- tempfile(); dir.create(d)
  file.copy(c("data/covariates.tif", "data/records.csv"), d)
  if (!is.null(records)) write.csv(records, file.path(d, "records.csv"), row.names = FALSE, na = "")
  if (!is.null(covs)) writeRaster(covs, file.path(d, "covariates.tif"), overwrite = TRUE)
  d
}
# The simulated polygon layers, split into two files.
pdir <- tempfile(); dir.create(pdir)
outl <- st_read("data/admin.gpkg", quiet = TRUE)
for (ly in c("parish", "facility")) st_write(outl[outl$layer == ly, ], file.path(pdir, paste0(ly, ".gpkg")), quiet = TRUE)
polys <- list(list(file = file.path(pdir, "parish.gpkg"),   id_col = "unit_id", name_col = "unit_name", prefix = ""),
              list(file = file.path(pdir, "facility.gpkg"), id_col = "unit_id", name_col = "unit_name", prefix = ""))
run <- function(dir, ..., polygons = polys) {
  txt <- capture.output(r <- check_inputs(dir, polygons = polygons, ...))
  list(text = paste(txt, collapse = "\n"), r = r)
}
has <- function(x, pattern) grepl(pattern, x$text)
recs <- read.csv("data/records.csv", stringsAsFactors = FALSE, colClasses = c(unit_id = "character"))
covs <- rast("data/covariates.tif")

# ---- clean data: no problems, no warnings, and every section present -----------------
ok <- run("data")
chk(length(ok$r$problems) == 0, paste("clean data should have no problems:", toString(ok$r$problems)))
chk(length(ok$r$warnings) == 0, paste("clean data should have no warnings:", toString(ok$r$warnings)))
for (s in c("Mask and covariates", "Polygon layers", "Records", "Scoring time", "Summary"))
  chk(has(ok, s), paste("missing section:", s))
chk(has(ok, "WGS84 longitude/latitude") && has(ok, "all layers are masked identically"), "mask lines")
chk(has(ok, "MESS and SHAPE for the whole mask"), "scoring time reported")

# ---- privacy: nothing identifying can appear -------------------------------------------
leak <- function(x, d = recs) {
  xy <- c(d$lon, d$lat); xy <- xy[!is.na(xy) & abs(xy - round(xy)) > 0.01]       # ignore round numbers that occur in text anyway
  c(record_id = any(vapply(d$record_id, function(i) grepl(i, x$text, fixed = TRUE), TRUE)),
    unit_id   = any(vapply(unique(d$unit_id), function(i) grepl(i, x$text, fixed = TRUE), TRUE)),
    unit_name = grepl("Parish|Facility", x$text),
    coords    = any(vapply(xy, function(v) grepl(sprintf("%.3f", v), x$text, fixed = TRUE), TRUE)))
}
chk(!any(leak(ok)), paste("clean report leaks:", toString(names(which(leak(ok))))))
# ...including when things go wrong, and when the error text could carry ids or names
dup <- outl[outl$layer == "parish", ]; dup$unit_id[2] <- dup$unit_id[1]
dupf <- file.path(pdir, "dup.gpkg"); st_write(dup, dupf, quiet = TRUE, delete_dsn = TRUE)
bad <- run("data", polygons = list(list(file = dupf, id_col = "unit_id", name_col = "unit_name", prefix = "")))
chk(has(bad, "duplicate ids") && length(bad$r$problems) > 0, "duplicate polygon ids are a problem")
chk(!any(leak(bad)), paste("report leaks when polygon ids are duplicated:", toString(names(which(leak(bad))))))
broken <- mkdir(); writeLines("not,a,raster", file.path(broken, "covariates.tif"))
bk <- run(broken)
chk(has(bk, "could not be checked") && !any(leak(bk)), "an unreadable mask is reported")
chk(!has(bk, "could not be checked:") && !has(bk, "GDAL|not recognized|supported file format"), "...without the error text, which could hold ids or names")
chk(has(run(broken, details = TRUE), "could not be checked: "), "details = TRUE shows the error for the user's own console")

# ---- the mask --------------------------------------------------------------------------
proj <- mkdir(covs = terra::project(covs, "EPSG:32636"))
pj <- run(proj, polygons = NULL)
chk(any(grepl("not in longitude/latitude", pj$r$problems)), "a projected mask is a problem")
nonsq <- rast(nrows = 30, ncols = 60, nlyrs = 3, xmin = 30, xmax = 33, ymin = -3, ymax = 0, crs = "EPSG:4326")
values(nonsq) <- matrix(rnorm(30 * 60 * 3), ncol = 3); names(nonsq) <- names(covs)     # 0.05 x 0.1 degree cells
chk(any(grepl("not square", run(mkdir(covs = nonsq), polygons = NULL)$r$warnings)), "non-square cells are noted")

const <- covs; const[[2]] <- const[[2]] * 0 + 1
chk(any(grepl("constant covariate", run(mkdir(covs = const), polygons = NULL)$r$problems)), "a constant covariate is a problem")
dupl <- c(covs, covs[[1]] * 2); names(dupl) <- c(names(covs), "twice_elevation")
chk(any(grepl("highly correlated", run(mkdir(covs = dupl), polygons = NULL)$r$warnings)), "near-duplicate layers are a warning")
gap <- covs; gap[[3]][500] <- NA; gap[[1]][500] <- 1                      # a covariate gap inside the mask
chk(any(grepl("not masked identically", run(mkdir(covs = gap), polygons = NULL)$r$problems)), "inconsistent masks are a problem")

# ---- other raster layers ---------------------------------------------------------------
ld <- tempfile(); dir.create(ld)
good <- covs[[1]]; values(good) <- sample(1:6, ncell(good), replace = TRUE); good <- mask(good, covs[[1]]); names(good) <- "landcover"
writeRaster(good, file.path(ld, "landcover.tif"))
fine <- rast(nrows = 120, ncols = 120, xmin = 30, xmax = 33, ymin = -3, ymax = 0, crs = "EPSG:4326"); values(fine) <- 1
writeRaster(fine, file.path(ld, "dwellings_fine.tif"))
off <- rast(nrows = 60, ncols = 60, xmin = 31, xmax = 34, ymin = -3, ymax = 0, crs = "EPSG:4326"); values(off) <- 1
writeRaster(off, file.path(ld, "shifted.tif"))
ly <- run("data", layers_dir = ld, polygons = NULL)
chk(has(ly, "landcover.tif: aligned; looks categorical, 6 distinct values"), "an aligned categorical layer is recognised")
chk(any(grepl("dwellings_fine.tif is not aligned.*cell size differs", ly$r$problems)), "a different cell size is a problem")
chk(any(grepl("shifted.tif is not aligned.*extent differs", ly$r$problems)), "a different extent is a problem")
chk(!any(leak(ly)), "layer report leaks nothing")

# ---- records ----------------------------------------------------------------------------
extra <- data.frame(record_id = sprintf("X%02d", 1:10), unit_id = "parish:01",
                    lon = c(0.123, 0.123, NA, NA, NA, NA, NA, NA, NA, NA), lat = c(0.456, 0.456, NA, NA, NA, NA, NA, NA, NA, NA),
                    n_specimens = 1, date = c("12/05/2019", "12/05/2019", "12/05/2019", "12/05/2019", "12/05/2019", "12/05/2019", "2020", "2020-03", "2020-03-04", NA),
                    stringsAsFactors = FALSE)
extra$lon[3] <- -1.234; extra$lat[3] <- 31.567                            # swapped longitude and latitude
r2 <- rbind(recs, extra)
rr <- run(mkdir(records = r2), polygons = polys)
chk(has(rr, "dates: .*other format 6") && has(rr, "Patterns seen .*: 99/99/9999"), "an unrecognised date format is flagged with its pattern, not its values")
chk(has(rr, "outside the raster extent: <5"), "small counts are suppressed by default")
chk(has(rr, "<5 of them would be inside the mask if longitude and latitude were swapped"), "swapped coordinates are detected")
rr0 <- run(mkdir(records = r2), polygons = polys, suppress_below = 0)
# three records lie outside the raster extent: two near (0, 0) and the swapped one
chk(has(rr0, "outside the raster extent: 3") && has(rr0, "1 of them would be inside the mask if longitude and latitude were swapped"), "suppression can be turned off")
chk(!has(rr, "outside the raster extent: 3"), "the exact small count is not shown when suppressed")
dups <- recs; dups$record_id[2] <- dups$record_id[1]
dd <- run(mkdir(records = dups), polygons = NULL)
chk(any(grepl("duplicate record ids", dd$r$problems)), "duplicate record ids are a problem")
chk(!any(leak(dd)), "duplicate record ids are counted, never listed")
chk(!any(leak(rr, r2)), "record problems leak nothing")

# names only: matched, ambiguous and unknown, as counts
named <- recs[, c("record_id", "lon", "lat", "n_specimens", "date")]
un <- read.csv("data/units.csv", stringsAsFactors = FALSE)
named$admin_name <- un$unit_name[match(recs$unit_id, un$unit_id)]
named$admin_name[1:6] <- "Nowhere"
nm <- run(mkdir(records = named), polygons = polys)
chk(has(nm, "named by admin_name: 100; matched to exactly one unit: 94"), "name matching is summarised as counts")
chk(has(nm, "matching no unit: 6"), "unknown names counted")
chk(!any(leak(nm)) && !grepl("Nowhere", nm$text), "name matching leaks no names")

# ---- polygon layers ---------------------------------------------------------------------
nocrs <- outl[outl$layer == "parish", ]; st_crs(nocrs) <- NA
# a shapefile with no .prj (reads back as NA), and a GeoPackage written without a CRS (reads back as "Undefined Cartesian")
nf <- file.path(pdir, "nocrs.shp"); suppressWarnings(st_write(nocrs[, c("unit_id", "unit_name")], nf, quiet = TRUE, delete_dsn = TRUE))
unlink(sub("\\.shp$", ".prj", nf))
chk(is.na(st_crs(st_read(nf, quiet = TRUE))), "(setup) the shapefile has no CRS")
ng <- file.path(pdir, "nocrs.gpkg"); suppressWarnings(st_write(nocrs[, c("unit_id", "unit_name")], ng, quiet = TRUE, delete_dsn = TRUE))
for (f in c(nf, ng)) {
  np <- run("data", polygons = list(list(file = f, id_col = "unit_id", name_col = "unit_name", prefix = "")))
  chk(any(grepl("has no coordinate system", np$r$problems)), paste("a polygon layer with no CRS is a problem:", basename(f)))
}
chk(any(grepl("file not found", run("data", polygons = list(list(file = "nope.gpkg", id_col = "a", name_col = "b")))$r$problems)), "a missing polygon file is a problem")
chk(any(grepl("no id column", run("data", polygons = list(list(file = polys[[1]]$file, id_col = "missing", name_col = "unit_name", prefix = "")))$r$problems)),
    "a missing id column is a problem")
# a polygon outside the mask: reported as a unit with no cells, and no ids shown
far <- outl[outl$layer == "parish", ][1, ]; far$unit_id <- "zzfar9"
st_geometry(far) <- st_sfc(st_polygon(list(rbind(c(10, 10), c(11, 10), c(11, 11), c(10, 11), c(10, 10)))), crs = 4326)
farf <- file.path(pdir, "far.gpkg"); st_write(rbind(outl[outl$layer == "parish", ], far), farf, quiet = TRUE, delete_dsn = TRUE)
fr <- run("data", polygons = list(list(file = farf, id_col = "unit_id", name_col = "unit_name", prefix = "")))
chk(has(fr, "units with no cells in the mask: <5") && any(grepl("units with no cells cannot be ranked", fr$r$warnings)), "units with no cells are reported")
chk(!has(fr, "zzfar9"), "...as a count, never by id")

# ---- scoring time -----------------------------------------------------------------------
few <- recs; few$lon[!is.na(few$lon)][-(1:3)] <- NA; few$lat[!is.na(few$lat)][-(1:3)] <- NA
chk(any(grepl("fewer than 8 located records", run(mkdir(records = few), polygons = NULL)$r$warnings)), "too few located records is a warning")
chk(!has(run("data", polygons = NULL, time_scoring = FALSE), "Scoring time"), "timing can be switched off")

cat("All input-check tests passed\n")
