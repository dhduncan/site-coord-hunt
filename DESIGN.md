# Design: prioritising coarse records for geolocation

Status: draft for discussion between spatial modellers and entomologists.
The reasons behind each decision are summarised in [DECISIONS.md](DECISIONS.md).

## Problem

Surveillance records are often located only to a parish, village, health
facility or other admin unit. Many of these units are larger than the cell size
of the model, so a record cannot be assigned to a single cell. Re-locating
every record is expensive. We want to rank the coarse records by how much a
precise location would help the model, and then make the re-location itself
quick.

## Assumed model objective

An environmental-niche / species distribution model fitted at the resolution of
a raster mask. For this objective the most valuable new site is one that adds
**novel environmental space** relative to the sites that are already located at
model resolution. A record in a conditions-poor corner of covariate space
constrains the response curves more than one that duplicates conditions we
already have. If the objective changes (e.g. spatial prediction error, range
edge delineation), the value metric changes and this document should be revisited.

## Definitions

- **Reference set**: records geolocated at model resolution (a cell, or a point
  that falls in one). Covariates are extracted at these cells.
- **Unit**: a coarse place that records can be attributed to (a parish, a facility
  catchment, a village), identified by a `unit_id`. Units can overlap: a facility
  catchment can lie inside a parish, and some records are attributed to one and some
  to the other.
- **Coarse record**: a record with only a unit as its location.
- **Candidate cells**: the mask cells a unit could be in, each with a weight (see
  "Converting polygons to cells").
- **Single-candidate unit**: a unit whose polygon covers exactly one mask cell. Its
  records are located at that cell automatically, because there is nothing left to
  resolve.
- **Narrowed record**: a coarse record for which the entomologist has marked one
  or more candidate cells, possibly with weights, without pinning it to one cell.
- **Novelty of a cell**: how far the covariates of that cell lie outside the
  reference set.

## Value metrics

Two complementary extrapolation metrics, both computed per cell against the
reference set:

1. **MESS** (Elith et al. 2010, Multivariate Environmental Similarity Surface).
   Univariate: for each covariate, how far a value is outside the range of the
   reference set, then the minimum across covariates. Negative = at least one
   covariate is outside the reference range. It cannot see unusual
   *combinations* of covariates that are individually in range.
2. **SHAPE** (Velazco et al. 2023, via `flexsdm::extra_eval()`). Distance-based
   (Mahalanobis by default): the distance from a cell to its nearest reference
   site, as a percentage of the mean distance of reference sites to their own
   centroid. 0 means the cell matches a reference site; it is unbounded, and
   values above 100 are farther from any reference site than the reference
   cloud's typical spread. Because it is distance-based it can flag unusual
   combinations of covariates that are each within range.

We show both and do not merge them into one score by default. The app uses
SHAPE as the primary ranking value and MESS (< 0) as the flag for "outside the
observed range of at least one covariate". We do not use flexsdm's `univar_comb`
layer: in flexsdm 1.4.2 its codes are the reverse of its help text (checked
against MESS on the simulated data), so range status is taken from MESS. Both are
reported per cell, and the modeller can change which one drives the ranking.

## From cells to records

A coarse record spans many cells, and we do not know which one is right. For
every coarse record we score each of its unit's candidate cells, then report:

- `max_novelty`: best-case payoff, i.e. if the true location is the most novel cell
- `mean_novelty`: expected payoff, the mean over candidate cells weighted by how
  likely each is (the share of the polygon in it)
- `range`: how much the answer depends on where the record really is
  (`max - min`)
- `n_cells`: number of candidate cells
- `n_records`: number of records sharing the unit (geolocating the unit once helps
  all of them if they are from the same site; flagged separately)

Default ranking is by `max_novelty`; records with a large range are the ones
where a good location helps most and a wrong guess hurts most. Where
records carry a **plausible-habitat constraint** (e.g. a land-cover mask) the
candidate cells can be restricted, which narrows the range.

## Sequential selection

Novelty is relative to the reference set, so every site that is geolocated
changes the scores of all other candidates. The app therefore:

- recomputes scores whenever the reference set changes, and
- offers a "next best" ranking in which a candidate unit that is similar in
  covariate space to an already-selected one is down-weighted (greedy
  selection that adds each chosen record's best cell to a *working* reference set
  and rescores).

This avoids sending effort to ten adjacent parishes that are all novel in the
same way.

## Two views

The same data serves two people with different jobs, so the app has two views,
selected by a **role tab** (Modeller / Entomologist) in a single Shiny app. They
share the scoring and storage code but show different things.

The tab is a convenience, not access control: anyone can switch to the Modeller tab.
That is adequate against unconscious drift toward high-value cells, which is the
concern. If hiding value ever had to be enforced, the entomologist view would need
a separate deployment. The Entomologist tab is the landing tab, so opening the app
does not show value information.

The entomologist's work list is ordered by the modeller's scores but carries none.
If scoring is not possible (for example, too few located sites), the list falls back to
A to Z so fieldwork is never blocked by a modelling problem.

| | Modeller view | Entomologist view |
|---|---|---|
| Question | Which records would add the most to the model? | Where did this sample really come from? |
| Controls | Covariate subset, MESS vs SHAPE, reference set rules, how many locations to prioritise | None of the scoring settings |
| Sees | Cell and unit scores, environmental-space diagnostics, ranked table, how the ranking shifts as sites are added, the corrections log with source and confidence | A work list for their area of responsibility, location polygons, the model grid, located sites, and the recall aids below |
| Produces | A prioritised work list of locations | Corrections: one or more cells, with weights, source and confidence |

### Entomologist view: no value information

Entomologists never see the value of a cell to the model: no novelty colours, no
SHAPE or MESS numbers, no cell-level ranking. If they could see it, they would be
drawn, even unconsciously, toward high-value cells, and the location recorded
would reflect what is useful to the model and not where the sample came from. That
would bias the very data the model is fitted to. Their only guide to a location is
what helps them remember it.

The work list is ordered at the level of the location only (that is the purpose of
the prioritisation) and carries no scores. For each location ID it shows what
helps recall:

- the **number of records** attributed to it
- the **time span** of those records as month and year, e.g. "Mar 2018 to Nov
  2021" ("undated" where there is no date, and year only where only a year is known)

The location's outline and the cells it could be in are shown (geometry only, with no
scores), the model grid is drawn for the cells in view (zoom in to see it on a large
study area), and so are the sites
that are already located (an exact point is drawn as a dot, a site known only to its
cell as that cell's square, so the map shows what was recorded). Showing located sites helps recall, and the hope is that
entomologists notice ones that look wrong. For that to be useful there is a **flag
as possibly wrong** action on a located site (a note is optional). Flags go to the
modeller view. By default flagged sites are left out of the reference set until the
modeller reviews and dismisses the flag (a setting the modeller can turn off), because
a misplaced reference site distorts the novelty scores of every other candidate.
Re-locating a flagged site is not yet supported: the flag only reaches the modeller. A neutral base map (satellite or terrain) is a recall aid and is fine.

### Marking a location

For a record in the list, the entomologist can:

1. **Click one cell**: the record is located at model resolution (`method = "cell"`).
   Only the cell is recorded: no coordinates are stored, and none are shown to the
   entomologist.
2. **Click to drop a marker and drag it** to the exact spot: a precise location
   (`method = "point"`), from which the cell is derived.
3. **Mark several cells** when unsure, then say how belief is spread over them:
   - **Evenly split**: every marked cell counts equally (exactly 1/n, so two cells are
     50/50 and three are one third each). This is also the way to say "somewhere in
     these cells", with no limit on how many.
   - **Manual weights**: the entomologist gives each cell a percentage. The total
     must be 100, and the interface keeps it so: the last cell takes the remainder,
     so the sum is always right. The one thing that can still go wrong is the other
     cells adding up to 100 or more; the interface then says so and blocks saving
     until every cell has at least 1%. Whole percentages only, to avoid claiming more
     precision than anyone has.

   Clicking one cell is the single-cell case (weight 1). Weights are stored as
   proportions that sum to 1 and shown as percentages, so an even split over three
   cells is stored as 1/3 each and not as 33, 33, 33. Overall doubt is expressed by the
   separate confidence field, not by the weights.

Each save records a source and a confidence level. A new save for the same record
replaces the earlier one in full; the earlier rows stay in the log.

### What a multi-cell record means for scoring

- **The candidate cells become the marked cells.** A narrowed record is scored over
  only its marked cells, and `mean_novelty` becomes a weighted mean using the
  entomologist's weights. `max_novelty` and `range` are taken over marked cells.
  So marking cells shrinks the polygon-wide uncertainty rather than ignoring it.
- **Narrowed is not located.** A record with several cells is not added to the
  reference set, because we do not know which cell's environment it represents. The
  modeller sets the rule for when a narrowed record counts as located: one cell
  holds at least a chosen share of the weight (default 80%, a slider in the modeller
  view), and that cell is then used as its location. Narrowed records below the rule stay
  on the work list with their reduced candidate set.
- **The model sees the uncertainty too.** Fitting the model with a record that
  has several weighted cells is a modelling choice (weighted presences, or
  averaging over cells) and is outside this repo, but the corrections file
  provides the cells and weights it needs.

The entomologist sees none of this scoring. They see only that a record is
located, narrowed or still coarse.

Corrections are never written over the original data; see the next section.

## Converting polygons to cells

Approximate locations arrive as polygons. The model works on the raster mask, so each
unit is converted once to its candidate cells (`data-raw/build_candidates.R`, using
`R/candidates.R`). The rules, all visible in the output:

- **Weight = the share of the polygon's area that falls in the cell.** This assumes the
  true location is equally likely anywhere within the polygon (uniform; ignoring
  habitat or access, which a prior on probable location could refine under the rule in
  `DECISIONS.md` 6d; none is built yet). It is not the fraction of
  the cell covered: a large polygon has many cells, each covered fully, and they should
  not all weigh the same as a sliver at the edge.
- **Normalised over mask cells only.** Weights sum to 1 for each unit. A polygon that
  extends past the mask (coast, border, lake) keeps only its in-mask part, and
  `mask_share` records how much that was.
- **Overlap is allowed.** A cell can be a candidate for many units; the table is
  many-to-many. Unit ids are unique across layers (a `prefix` per layer).
- **Small polygons are not dropped.** A polygon smaller than a cell still gets the
  cell(s) it touches, with weights in proportion to its area in each. If that is one
  cell, the unit is a single-candidate unit.
- **Problems are listed, never dropped silently:** units with no cells in the mask,
  units mostly outside it, and invalid geometries (repaired and reported). The Modeller
  tab also lists records whose unit is unknown or has no cells.
- **The mask is fixed**, so the cell IDs in the table stay valid.

## Data

All inputs are plain files so the repo runs without a database. The storage layer
(`R/data_io.R`) is the only code that touches files, so it can be swapped for a
database later.

- `covariates.tif`: covariate stack (multi-layer GeoTIFF), `NA` = outside the mask. The
  grid is the cell structure of this raster. The mask is fixed: it is the maximal
  extent of the project, and every resource is masked and aligned to it, so cell IDs
  are stable and can be stored and used directly.
- `candidates.csv`: `unit_id`, `cell_id`, `weight`. Which mask cells each unit could be
  in, and how likely each is. Built from the polygon layers by
  `data-raw/build_candidates.R`. The app refuses a table that does not fit the mask
  (cell IDs outside the raster, cells outside the mask, weights that do not sum to 1
  for a unit, a unit missing from `units.csv`).
- `units.csv`: `unit_id`, `unit_name`, `layer`, `n_cells`, `mask_share` (the share of
  the polygon inside the mask).
- `admin.gpkg`: the unit outlines, with a `unit_id` column, for drawing on the map
  only. Nothing is scored from it.
- `records.csv`: one row per record.

| Column | Meaning |
|---|---|
| `record_id` | unique id |
| `unit_id` | the unit the record is attributed to. Where the source data has only names, give `admin_name` instead (and `unit_layer` if the same name occurs in more than one layer); the app matches names to units, and reports any that match none or more than one, never guessing |
| `lon`, `lat` | precise coordinates; blank if not yet located |
| `n_specimens` | optional, count |
| `date` | recommended: shown in the entomologist view as month and year; partial dates (year only, or year and month) are accepted |

- `corrections.csv`: append-only log written by the app. One row per marked
  cell, grouped by `action_id` (one save = one action). Columns: `action_id`,
  `record_id`, `cell_id`, `lon`, `lat`, `weight` (proportion; sums to 1 per
  action), `weight_mode` (`even` or `manual`), `method` (`cell`, `point` or
  `cells`), `source`, `confidence`, `user`, `timestamp`. A single cell or point has
  `weight = 1`. **`lon` and `lat` are filled only for an exact point** (`method =
  point`); a marked cell is identified by `cell_id` alone, so a cell centre is never
  stored as if it were a measured location. Where the app needs a centre (to score a
  cell, or to draw it) it derives it from the cell ID, and the file is refused if a
  cell-level row carries coordinates. Older logs that stored centres are blanked on
  read. The current state of a record is its latest action; the original
  `records.csv` is never modified, so every edit is auditable and reversible. Files
  written by the first prototype (one row per record, no `action_id` or `weight`)
  are still read, each row counting as a single-cell action.

- `flags.csv`: append-only log of flag actions on located sites: `action` (`flag`
  by an entomologist, `dismiss` by a modeller), `record_id`, `user`, `timestamp`,
  `note`. A site is currently flagged if its latest action is `flag`. Flags never
  change a record's location.

A record counts as **located** when it has one location (a cell or a point) that
falls inside a mask cell. A record with several marked cells is **narrowed** (see
"What a multi-cell record means for scoring"). Names that do not match any polygon
are reported, not silently dropped.

## Simulated data

The first prototype runs on simulated covariates, admin units and records
(`data-raw/simulate_data.R`), so the repo needs no private data. The simulation
deliberately includes a few admin units that cover unusual environmental
conditions, so the ranking has something to find. There are two layers of units: 36
parishes that line up with the grid, and 8 facility catchments that do not and that
overlap the parishes. Two catchments are smaller than a cell, so their records are
placed automatically. Record dates are simulated with the
mix of precision real records have (full dates, month only, year only, and a few
missing).

## Open questions

1. Real names may match the polygon layers imperfectly. Unmatched and ambiguous names
   are now reported, but facilities that exist only as a point (with a catchment
   radius) still need converting to polygons, for example by buffering.
2. Should novelty be computed against all covariates in the model, or a chosen subset?
3. Do records from the same unit and date describe one site, or several?
4. What is an acceptable georeferencing precision (cell centre vs. point), and does
   the entomologist's own certainty need to feed the model as uncertainty?
5. The work list is a live queue that re-ranks as sites are located (working
   default). Is that right, or does the modeller need to publish a fixed list?
6. A narrowed record counts as located once one cell holds at least 80% of the
   weight (working default, a placeholder). What share is right?
7. How should a flagged site get re-located? Today a flag only reaches the modeller.
8. If the final model is not niche-based (e.g. a spatial smoother), the value
   metric needs to change (see "Assumed model objective").

## References

- Elith J, Kearney M, Phillips S (2010). The art of modelling range-shifting
  species. *Methods in Ecology and Evolution* 1:330-342.
- Velazco SJE, et al. (2023). SHAPE extrapolation metric. *Ecography*,
  doi:10.1111/ecog.06992 (as cited in the `flexsdm::extra_eval` documentation).
- Velazco SJE, et al. (2022). flexsdm: An R package for supporting a
  comprehensive and flexible species distribution modelling workflow. *Methods in
  Ecology and Evolution*.
