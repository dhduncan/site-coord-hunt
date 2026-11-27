# Run from repo root: Rscript tests/test_data_io.R
suppressPackageStartupMessages({library(terra); library(sf)})
source("R/dates.R"); source("R/scoring.R"); source("R/candidates.R"); source("R/data_io.R")
chk <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)

inp <- load_inputs("data")
covs <- inp$covs; cand <- inp$candidates
single_units <- names(which(table(cand$unit_id) == 1))
n_orig <- sum(!is.na(inp$records$lon))                          # records that came with coordinates
auto_ids <- inp$records$record_id[inp$records$unit_id %in% single_units & is.na(inp$records$lon)]
n_ref0 <- n_orig + length(auto_ids)                              # located before anyone marks anything
chk(length(auto_ids) > 0, "the simulated data should include units that cover a single cell")
st0 <- record_states(inp$records, empty_corrections(), covs, candidates = cand)
chk(sum(st0$records$status == "located") == n_ref0, "located = records with coordinates + single-cell units")
chk(sum(st0$records$status == "coarse") == nrow(inp$records) - n_ref0, "the rest are coarse")
pr <- unit_problems(inp$records, inp$units, cand)
chk(length(pr$unknown) == 0 && length(pr$no_cells) == 0, "unexpected unit problems")
chk(nrow(inp$record_problems) == 0, "unexpected record problems")

coarse <- st0$records$record_id[st0$records$status == "coarse"]
a <- coarse[1]; b <- coarse[2]
cell_xy <- function(cell) { xy <- xyFromCell(covs[[1]], cell); list(lon = xy[1], lat = xy[2]) }
# Only an exact point carries lon/lat (here: the cell centre, a valid point inside
# the cell). Cell-level marks are identified by cell_id alone.
act <- function(id, cells, w, mode = "manual", method = "cells") {
  xy <- lapply(cells, cell_xy)
  exact <- method == "point"
  data.frame(action_id = new_action_id(), record_id = id, cell_id = cells,
             lon = if (exact) sapply(xy, `[[`, "lon") else NA_real_,
             lat = if (exact) sapply(xy, `[[`, "lat") else NA_real_,
             weight = w, weight_mode = mode, method = method, source = "t",
             confidence = "high", user = "t", timestamp = "x")
}
dir <- tempfile(); dir.create(dir)

# several cells: narrowed until one cell is heavy enough
valid <- which(!is.na(values(covs[[1]], mat = FALSE)))
cs <- valid[c(100, 500, 900)]
append_correction(act(a, cs[1:2], c(0.5, 0.5), "even"), dir)
st <- record_states(inp$records, read_corrections(dir), covs, candidates = cand)
chk(st$records$status[st$records$record_id == a] == "narrowed", "50/50 should be narrowed")
chk(nrow(st$cells) == 2 && all(st$cells$weight == 0.5), "marked cells not kept")
chk(!anyNA(st$cells$lon) && !anyNA(st$cells$lat), "marked cells should get derived centres")

# A threshold must not let which.max() choose among tied cells.
tie_dir <- tempfile(); dir.create(tie_dir)
append_correction(act(a, cs[1:2], c(0.5, 0.5), "even"), tie_dir)
tie <- record_states(inp$records, read_corrections(tie_dir), covs,
                     located_share = 0.5, candidates = cand)
chk(tie$records$status[tie$records$record_id == a] == "narrowed",
    "a tied 50/50 split must stay narrowed at a 50% threshold")
chk(nrow(tie$cells) == 2, "both tied cells must remain available")

three_way_dir <- tempfile(); dir.create(three_way_dir)
append_correction(act(a, cs, rep(1 / 3, 3), "even"), three_way_dir)
three_way <- record_states(inp$records, read_corrections(three_way_dir), covs,
                           located_share = 0.3, candidates = cand)
chk(three_way$records$status[three_way$records$record_id == a] == "narrowed",
    "three equal cells must stay narrowed at a lower threshold")

dominant_dir <- tempfile(); dir.create(dominant_dir)
append_correction(act(b, cs, c(0.5, 0.3, 0.2), "manual"), dominant_dir)
dominant <- record_states(inp$records, read_corrections(dominant_dir), covs,
                          located_share = 0.5, candidates = cand)
chk(dominant$records$status[dominant$records$record_id == b] == "located",
    "a unique cell at the threshold should still count as located")

append_correction(act(b, cs[1:2], c(0.85, 0.15)), dir)
st <- record_states(inp$records, read_corrections(dir), covs, located_share = 0.8, candidates = cand)
rb <- st$records[st$records$record_id == b, ]
chk(rb$status == "located", "85% on one cell should count as located at 0.8")
chk(isTRUE(all.equal(c(rb$lon, rb$lat), unlist(cell_xy(cs[1])), check.attributes = FALSE)),
    "located record should use the heaviest cell")
st9 <- record_states(inp$records, read_corrections(dir), covs, located_share = 0.9, candidates = cand)
chk(st9$records$status[st9$records$record_id == b] == "narrowed", "threshold must be adjustable")

# a later save replaces the earlier one in full, in both directions
append_correction(act(a, cs[3], 1, "even", "cell"), dir)
st <- record_states(inp$records, read_corrections(dir), covs, candidates = cand)
chk(st$records$status[st$records$record_id == a] == "located", "single cell should replace narrowed")
chk(!(a %in% st$cells$record_id), "old marked cells should be gone")
append_correction(act(a, cs[1:3], rep(1/3, 3), "even"), dir)
st <- record_states(inp$records, read_corrections(dir), covs, candidates = cand)
chk(st$records$status[st$records$record_id == a] == "narrowed", "narrowing again should replace located")
chk(sum(st$cells$record_id == a) == 3, "expected 3 marked cells")

# input validation and the original data
err <- tryCatch(append_correction(act(a, cs[1:2], c(0.5, 0.6)), dir), error = function(e) e)
chk(inherits(err, "error"), "weights that do not sum to 1 must be refused")
chk(is.na(inp$records$lon[inp$records$record_id == a]), "original records mutated")

# files from the first prototype (no action id or weight) still load
old <- tempfile(); dir.create(old)
writeLines(c('"record_id","lon","lat","cell_id","method","source","confidence","user","timestamp"',
             sprintf('"%s",%f,%f,%d,"cell","x","high","t","y"', a, cell_xy(cs[1])$lon, cell_xy(cs[1])$lat, cs[1])),
           file.path(old, "corrections.csv"))
lc <- read_corrections(old)
chk(lc$weight == 1 && nzchar(lc$action_id), "legacy corrections not upgraded")
chk(is.na(lc$lon) && is.na(lc$lat), "legacy cell-level centres should be blanked on read")
chk(record_states(inp$records, lc, covs, candidates = cand)$records$status[inp$records$record_id == a] == "located", "legacy row should locate")

# cell-level marks store no coordinates; the centre is derived from the cell id
dc <- tempfile(); dir.create(dc)
append_correction(act(a, cs[1], 1, "even", "cell"), dc)
raw <- utils::read.csv(file.path(dc, "corrections.csv"))
chk(is.na(raw$lon) && is.na(raw$lat), "a cell-only mark must be stored without coordinates")
stc <- record_states(inp$records, read_corrections(dc), covs, candidates = cand)
ra <- stc$records[stc$records$record_id == a, ]
chk(ra$status == "located" && !ra$exact, "a cell-level location is located but not exact")
chk(isTRUE(all.equal(c(ra$lon, ra$lat), unlist(cell_xy(cs[1])), check.attributes = FALSE)),
    "the centre should be derived from the cell id")
chk(nrow(reference_xy(stc, character(0))) == n_ref0 + 1, "a cell-level site still joins the reference set")
chk(all(st0$records$exact[st0$records$basis %in% "record coordinates"]), "coordinates that came with a record are exact")
chk(!any(st0$records$exact[st0$records$status != "located"]), "unlocated records are never exact")

# an exact point keeps its own coordinates, which differ from the cell centre
pt <- act(b, cs[2], 1, "even", "point"); pt$lon <- pt$lon + 0.01; pt$lat <- pt$lat - 0.01
append_correction(pt, dc)
stp <- record_states(inp$records, read_corrections(dc), covs, candidates = cand)
rp <- stp$records[stp$records$record_id == b, ]
chk(rp$status == "located" && rp$exact, "an exact point is located and exact")
chk(abs(rp$lon - pt$lon) < 1e-9 && abs(rp$lat - pt$lat) < 1e-9, "an exact point must keep its own coordinates")

# the rule is enforced on write
bad <- act(a, cs[1], 1, "even", "cell"); bad$lon <- 31; bad$lat <- -1
chk(inherits(tryCatch(append_correction(bad, dc), error = function(e) e), "error"), "a cell mark with coordinates must be refused")
bad <- act(a, cs[1], 1, "even", "point"); bad$lon <- NA_real_
chk(inherits(tryCatch(append_correction(bad, dc), error = function(e) e), "error"), "an exact point without coordinates must be refused")
bad <- act(a, cs[1], 1, "even", "cell"); bad$method <- "centroid"
chk(inherits(tryCatch(append_correction(bad, dc), error = function(e) e), "error"), "an unknown method must be refused")

# flags: flag, dismiss, flag again
fdir <- tempfile(); dir.create(fdir)
chk(length(current_flags(read_flags(fdir))) == 0, "no flags at start")
fl <- function(action, id) data.frame(action = action, record_id = id, user = "u", timestamp = "x", note = "")
append_flag(fl("flag", "R001"), fdir); append_flag(fl("flag", "R002"), fdir)
chk(setequal(current_flags(read_flags(fdir)), c("R001", "R002")), "two flags expected")
append_flag(fl("dismiss", "R001"), fdir)
chk(identical(current_flags(read_flags(fdir)), "R002"), "dismiss should clear a flag")
append_flag(fl("flag", "R001"), fdir)
chk(setequal(current_flags(read_flags(fdir)), c("R001", "R002")), "a new flag after dismissal counts")
chk(nrow(reference_xy(st0, c("R001", "R002"), TRUE)) == n_ref0 - 2, "flagged sites should leave the reference set")
chk(nrow(reference_xy(st0, c("R001", "R002"), FALSE)) == n_ref0, "exclusion must be optional")
chk(nrow(reference_xy(st0, character(0), TRUE)) == n_ref0, "no flags, nothing excluded")

# records that arrive with names (no ids) are matched to units; problems are reported, never guessed
dn <- tempfile(); dir.create(dn)
invisible(file.copy(file.path("data", c("covariates.tif", "candidates.csv", "admin.gpkg")), dn))
named <- inp$records[, c("record_id", "lon", "lat", "n_specimens", "date")]
named$admin_name <- inp$units$unit_name[match(inp$records$unit_id, inp$units$unit_id)]
named$admin_name[named$record_id == "R050"] <- "Nowhere"                  # matches no unit
un <- inp$units; un$unit_name[un$unit_id == "parish:02"] <- "Parish 01"    # now "Parish 01" names two units
utils::write.csv(named, file.path(dn, "records.csv"), row.names = FALSE, na = "")
utils::write.csv(un, file.path(dn, "units.csv"), row.names = FALSE, na = "")
ln <- load_inputs(dn)
chk(!anyNA(ln$records$unit_id[!(inp$records$unit_id %in% c("parish:01", "parish:02")) & named$record_id != "R050"]),
    "names that match exactly one unit are matched")
chk(is.na(ln$records$unit_id[ln$records$record_id == "R050"]), "an unknown name is left unmatched")
amb <- inp$records$record_id[inp$records$unit_id %in% c("parish:01", "parish:02")]
chk(all(is.na(ln$records$unit_id[ln$records$record_id %in% amb])), "an ambiguous name is left unmatched, not guessed")
chk(setequal(ln$record_problems$record_id, c("R050", amb)), "every unmatched record is reported")
sn <- record_states(ln$records, empty_corrections(), covs, candidates = ln$candidates)
chk(all(sn$records$status[sn$records$record_id %in% c("R050", amb) & is.na(sn$records$lon)] == "coarse"), "unmatched records stay coarse")
rsn <- record_summary(ln$candidates, score_cells(covs, reference_xy(sn, character(0))), sn)
chk(all(rsn$n_cells[rsn$record_id %in% c("R050", amb)] == 0), "unmatched records have no candidate cells and are not ranked")

# a unit that covers one cell locates its records automatically
ra <- st0$records[st0$records$record_id %in% auto_ids, ]
chk(all(ra$status == "located") && all(ra$basis == "single candidate cell") && !any(ra$exact),
    "single-cell units: located, basis recorded, not an exact point")
want <- xyFromCell(covs[[1]], cand$cell_id[match(ra$unit_id, cand$unit_id)])
chk(isTRUE(all.equal(cbind(ra$lon, ra$lat), want, check.attributes = FALSE)), "placed at that cell's centre")
chk(all(is.na(st0$records$basis[st0$records$status == "coarse"])), "coarse records have no basis")
chk(all(st0$records$basis[st0$records$record_id %in% inp$records$record_id[!is.na(inp$records$lon)]] == "record coordinates"),
    "records with coordinates say so")
chk(sum(record_states(inp$records, empty_corrections(), covs)$records$status == "located") == n_orig,
    "without the candidate table nothing is placed automatically")
da <- tempfile(); dir.create(da)
append_correction(act(auto_ids[1], cs[1], 1, "even", "cell"), da)
sta <- record_states(inp$records, read_corrections(da), covs, candidates = cand)
chk(sta$records$basis[sta$records$record_id == auto_ids[1]] == "entomologist" &&
    isTRUE(all.equal(sta$records$lon[sta$records$record_id == auto_ids[1]], xyFromCell(covs[[1]], cs[1])[1])),
    "an entomologist's mark takes precedence over automatic placement")

# ranking: narrowed records are scored over their marked cells only, weighted
st <- record_states(inp$records, read_corrections(dir), covs, located_share = 0.9, candidates = cand)
sc <- score_cells(covs, reference_xy(st, character(0)))
rs <- record_summary(cand, sc, st)
rb <- rs[rs$record_id == b, ]
cells_b <- st$cells[st$cells$record_id == b, ]
direct <- terra::extract(sc, cells_b$cell_id)
chk(rb$n_cells == 2, "narrowed record should be scored over its 2 marked cells")
chk(isTRUE(all.equal(rb$mean_shape, sum(cells_b$weight * direct$shape))), "mean_shape should be the weighted mean")
chk(isTRUE(all.equal(rb$max_shape, max(direct$shape))), "max_shape over marked cells")
u3 <- inp$records$unit_id[inp$records$record_id == coarse[3]]
chk(rs$n_cells[rs$record_id == coarse[3]] == sum(cand$unit_id == u3), "coarse record scored over all its candidate cells")
chk(!any(rs$record_id %in% st$records$record_id[st$records$status == "located"]), "located records must not be ranked")

# a coarse record's mean is weighted by the share of the polygon in each cell
fu <- names(which(table(cand$unit_id[grepl("^facility", cand$unit_id)]) > 4))[1]
fr <- st$records$record_id[st$records$unit_id == fu & st$records$status == "coarse"][1]
fc <- cand[cand$unit_id == fu, ]
chk(length(unique(round(fc$weight, 6))) > 1, "the test unit should have unequal weights")
fv <- terra::extract(sc, fc$cell_id)
chk(isTRUE(all.equal(rs$mean_shape[rs$record_id == fr], sum(fc$weight * fv$shape))), "coarse mean_shape should be weighted")
chk(isTRUE(all.equal(rs$max_shape[rs$record_id == fr], max(fv$shape))), "coarse max_shape over candidate cells")

u <- unit_summary(rs, st$records, inp$units)
chk(all(nzchar(u$dates)), "every location should have a date label")
chk(sum(u$n_records) == sum(st$records$status != "located"), "unit record counts should add up")
cat("All data_io tests passed\n")
