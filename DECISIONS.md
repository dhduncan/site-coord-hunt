# Decisions and their reasons

A short log of what has been decided so far and why, so the reasoning can be
challenged. [DESIGN.md](DESIGN.md) has the detail.

Status: **decided** (you chose it), **assumed** (we are working from it and it can
change), **working default** (accepted for now, easy to change), **open**.

## The question and how we answer it

| # | Decision | Reason | Status |
|---|---|---|---|
| 1 | Value of a record = how much novel environmental space a precise location would add, relative to the sites already located at model resolution. | For an environmental-niche model, a site in conditions we already sample adds little; one in novel conditions constrains the response curves. | assumed; revisit if the model objective changes |
| 2 | The model mask is a raster file, and its cells are the grid. | It ties cells directly to the covariate stack. | decided |
| 2a | The mask is fixed: the maximal extent of the project, with every resource masked and aligned to it. | Cell IDs are then stable, so they can be stored and used directly without checking for a changing grid. | decided |
| 3 | Novelty is measured with MESS and SHAPE (Velazco et al.), shown side by side. SHAPE ranks; MESS < 0 flags "outside the observed range of at least one covariate". | MESS cannot see unusual combinations of covariates that are each in range; SHAPE can. | decided |
| 4 | A coarse record is scored cell by cell over its unit's candidate cells, then summarised (max, weighted mean, range). | The spread inside a polygon is the point: a polygon with one very novel cell and many ordinary ones is different from a uniformly ordinary one. | decided |
| 5 | Scores are recomputed whenever a site is located. A "next best" ranking that down-weights similar candidates is planned. | Each located site changes the novelty of every other candidate; without this, effort goes to adjacent parishes that are all novel in the same way. | recompute: built; next best: not built |
| 5a | The modeller chooses which covariates feed the scores. | Not every covariate in a data pack belongs in the model. | built |
| 6 | Approximate locations are polygons, converted once to a candidate table (`unit_id`, `cell_id`, `weight`). Units can overlap, and records point to a unit by id. Where only names exist, they are matched to units; unknown or ambiguous names are reported, never guessed. | Matching by name silently merged two polygons that shared a name, and the polygon layer would be re-read on every score. A table is many-to-many, so a facility catchment inside a parish works, and it makes the overlap rule explicit. | decided; built |
| 6a | Weight = the share of the polygon's area in each cell (uniform within the polygon), normalised over mask cells. | The best we can say without more information is that the sample is equally likely anywhere in the polygon. Habitat could refine it later. | decided; built |
| 6b | A unit whose polygon covers exactly one mask cell has its records located at that cell automatically, and they join the reference set. | Nothing is left to resolve; asking an entomologist to click the only possible cell wastes their time. An entomologist's own mark still takes precedence. | decided; built |
| 6c | The model grid is drawn only for the cells in view, and only when that is a manageable number. | Drawing every cell of a real study area would swamp the browser (about 18 MB for 90,000 cells). | decided; built |
| 6d | **Circumspect use of priors on probable location.** Candidate weights start as the area share of the polygon (uniform within it). A prior on where a sample probably came from (land cover, dwellings, accessibility) is a single per-area layer, built from independent information and never from the niche covariates. It is invoked only if, on records whose location we know, it gives the known cell more probability on average (a higher log score) than a uniform prior does, judged on records not used to build it. Otherwise weights stay uniform within the polygon. Layers are not added because they seem plausible. | Novelty is highest where the species, and the sampling, is least likely: a summit is almost certainly novel and almost certainly not where the sample came from. A prior can correct that, but each one is another assumption, so it must earn its place. A prior fitted from the niche covariates would be circular: it would rate unfamiliar conditions as unlikely, which discounts exactly the novelty we are looking for. | agreed in principle; not built |
| 6e | With a prior in place, ranking is by expected novelty (the weighted mean), with the best cell kept as a diagnostic, and the ranking under a uniform prior is shown alongside for sensitivity. | The best cell ignores the weights, so it rewards improbable cells. Records that would expand the niche most may have a small prior chance, so the user should be able to see which locations depend on the prior. | agreed in principle; not built |

## Two people, one app

| # | Decision | Reason | Status |
|---|---|---|---|
| 7 | One app with a role tab (Modeller / Entomologist). No access control between the tabs. | One codebase, shared scoring and storage. Hiding value is aimed at unconscious drift, not at someone determined to look. | decided |
| 8 | Entomologists never see cell value (no novelty colours, scores or cell-level ranking). | If they could, their choice of cell would drift toward cells useful to the model, so the recorded location would no longer be where the sample came from, and the data the model is fitted to would be biased. | decided |
| 9 | The work list is ordered by location, with no scores. It shows the number of records and the month-and-year span for each location. | Ordering by location is what the prioritisation is for. Count and dates help recall where a site was. | decided |
| 10 | Entomologists see the sites already located, and can flag one as possibly wrong. Flags go to the modeller view, which can exclude flagged sites from the reference set until reviewed. | A person who knows the area can spot a misplaced site, and a wrongly placed reference site distorts novelty scores for everything else. A flag needs somewhere to go to be useful. | decided; built |
| 10a | Flagged sites leave the reference set by default until the modeller dismisses the flag; the modeller can turn the exclusion off. Re-locating a flagged site is not yet possible. | Keeping a doubtful site in the reference set quietly skews every score; leaving it out only makes the ranking slightly more cautious. | working default |
| 10b | The Entomologist tab is the landing tab. | Anyone opening the app starts on the view that shows no value information. | working default |

## How a location is recorded

| # | Decision | Reason | Status |
|---|---|---|---|
| 11 | Three ways to mark a location: click one cell, click to drop a marker and drag it to the exact spot, or mark several cells. A location chosen outside the record's own unit shows a warning but can still be saved. | Cell resolution is all the model needs; a point stores more when it is known; several cells express honest doubt. The recorded name may itself be wrong, so the warning does not block. | decided; built |
| 12 | With several cells, the entomologist chooses an even split or manual percentages. Manual weights always total 100, because the last cell takes the remainder. Weights are stored as proportions. | This replaced a fixed vote budget of 9, which cannot split evenly in two (5/4 reads as a lean) and caps the number of cells. Even split is exact, and auto-balancing means a bad sum never has to be rejected. Proportions avoid 33+33+33=99. | decided |
| 13 | Overall doubt is the separate confidence field, not the weights. | Weights are normalised to sum to 1, so they cannot carry "I'm not sure about any of this". | decided |
| 14 | A record with several marked cells is narrowed, not located, and does not join the reference set until one cell holds at least 80% of the weight (a slider for the modeller). A narrowed record is scored over its marked cells with a weighted mean. | We do not know which cell's environment it represents. Scoring over the marked cells shrinks the uncertainty without pretending it is gone. | working default (the 80% share is a placeholder); built |
| 14a | A mark on a cell stores only the cell ID. Coordinates are stored only for an exact point; centres are derived from the cell ID when needed. A cell-level site is drawn as its cell square and an exact point as a dot. | A dot at a cell centre looks like a point someone placed, and centres in the log look like measured coordinates. Scoring needs only the cell. | decided; built |
| 15 | Corrections go to an append-only log; the original records file is never modified. A later save for a record replaces the earlier one in full. | Every edit is auditable and reversible, and storage can be swapped without touching the app. | decided (storage: no preference expressed) |
| 16 | The work list is a live queue that re-ranks as sites are located, not a file the modeller publishes. If scoring is not possible it falls back to A to Z. | Matches decision 5. A published file would go stale after the first save, and a modelling problem should never stop fieldwork. | working default; built |

## How we build it

| # | Decision | Reason | Status |
|---|---|---|---|
| 17 | Prototype on simulated data first. | No private data needed to run the repo end to end. The simulation is biased so some units have conditions no located site shares, giving the ranking something to find. | decided |
| 18 | Real data stays outside the repo and on the machine it is held on, under its data-sharing agreement. It is checked and run locally; what comes back to the developer is a report of counts and summaries, never coordinates, record ids, place names or unit ids. Counts below 5 are shown as "<5" by default. | The Uganda data is under a strict data-sharing agreement. Whether it may be processed in a cloud session is a question for the agreement, so it is not attached to one. The report is built so that it cannot carry identifying detail, and a test checks that. | decided; built |

## What building taught us

- **SHAPE is not scaled 0 to 1.** `flexsdm::extra_eval()` returns a percentage with no upper bound. Colours on the map are capped at the 98th percentile so cells inside one unit stay distinguishable.
- **Shiny does not render outputs on a hidden tab.** Drawing commands sent to a map on the other tab were silently lost until both maps were built at start and refitted when their tab is shown.
- **flexsdm's `univar_comb` layer is reversed in version 1.4.2.** Its help text says 1 is univariate and 2 is combinatorial, but the code gives 2 where any covariate is outside the reference range. We take "outside range" from MESS instead.

## Not yet decided

The ranking: a unit is currently ranked by its best candidate cell. That ignores how likely the cell is, so it favours improbable cells (decision 6e), and it favours large polygons (a district will outrank a village; nested levels make that worse). The size effect is noted, not yet addressed.

See "Open questions" in [DESIGN.md](DESIGN.md): facilities that are only a point and a radius,
covariate subset, whether records from one unit and date are one site, acceptable
precision, and what changes if the model is not niche-based.
