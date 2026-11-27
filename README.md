# site-coord-hunt

Identifying high-value sites for geolocation.

Surveillance records located only to a parish, village or facility are often
larger than the model's cell size. Re-locating all of them is expensive. This repo
holds a Shiny app that ranks those coarse records by how much novel environmental
space a precise location would add to the set of sites already located (MESS and
SHAPE), and lets an entomologist record the location quickly.

## Concept

One app with a role tab and two views:

- **Modeller**: scoring settings, cell and unit scores, and the prioritised work
  list. The modeller also decides which flagged sites to keep in the reference set.
- **Entomologist**: the work list and the location tools. They see **no value
  information**, so the location they record is not drawn toward cells that are
  useful to the model. For each location they see the number of records and the
  month-and-year time span, to help recall. They also see the sites already located,
  and can flag one as possibly wrong.

A location can be one cell, a dragged point, or several cells with an even split or
manual percentage weights when the entomologist is unsure.

The tab is a convenience, not access control.

## Documents

- [DECISIONS.md](DECISIONS.md): what has been decided and why. Start here.
- [DESIGN.md](DESIGN.md): the detail: metrics, scoring of multi-cell records, data formats and open questions.

## Run the prototype

Needs R with `shiny`, `leaflet`, `terra`, `sf`, `DT` and `flexsdm`. The app opens on the
Entomologist tab; switch to Modeller with the tab at the top.

```r
source("data-raw/simulate_data.R")  # once; writes simulated data/ files
shiny::runApp()
```

To use real data you need the fixed mask (`covariates.tif`), your records, and your
polygon layers converted to cells. Edit and run `data-raw/build_candidates.R` to write
`candidates.csv`, `units.csv` and `admin.gpkg` (the units' outlines); then put those
with `covariates.tif` and `records.csv` in a folder (formats in DESIGN.md) and run
`options(sch.data_dir = "path/to/folder"); shiny::runApp()`.
Geolocations are appended to `corrections.csv` and flags to `flags.csv` in that folder;
the original `records.csv` is never modified. Dates in `records.csv` may be a full date,
a month (`2019-04`), a year (`2019`) or blank.

## Checking real data first

Real data should be checked before it is used, and it never needs to leave the machine it
is held on. Keep it outside the repo, then edit the settings in `data-raw/check_inputs.R`
(your data folder, the polygon layers, optionally a folder of other rasters) and run it
from the repo root:

```bash
Rscript data-raw/check_inputs.R
```

It changes nothing. It prints a report on the mask (it must be WGS84 longitude/latitude),
whether every layer lines up with it, the polygon layers, the records (how many fall inside
the mask, how many match a unit, what the dates look like), and how long one scoring run
takes on the real mask. The report holds **counts and summaries only**: no coordinates,
record ids, place names or unit ids, and counts below 5 are shown as `<5`
(`suppress_below` changes that). Check it against your data-sharing agreement before sharing
it. `details = TRUE` shows full error messages in your own console; do not share that output.

## Tests

```bash
Rscript tests/test_dates.R
Rscript tests/test_check_inputs.R
Rscript tests/test_candidates.R
Rscript tests/test_scoring.R
Rscript tests/test_data_io.R
```

## Status

The prototype runs on simulated data.

| Built | Designed, not built |
|---|---|
| Role tab with Modeller and Entomologist views | Limiting each entomologist's list to their own area |
| MESS and SHAPE scoring against located sites, with a choice of covariates | Sequential "next best" selection |
| Work list by location with record count and month-and-year span | Re-locating a flagged site |
| Click a cell, drag a point, or mark several cells (even split or manual percentages) | Habitat restriction of candidate cells |
| Scoring of narrowed records over their marked cells | Facilities that are only a point and a radius (they need converting to polygons first) |
| Polygons converted to weighted candidate cells; overlapping units; automatic placement of single-cell units | |
| Unmatched or ambiguous names reported, never guessed | |
| A local check of real inputs that reports counts only | |
| Located sites shown to entomologists, flagging, modeller review | Multi-user storage |
| Append-only corrections and flags logs; recompute after each save | |
| Simulated data (with dates) and tests | |

The map's basemap tiles need internet access.
