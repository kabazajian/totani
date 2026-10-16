# Two manuscript figures: owner run and review

This implements the approved template overview and regional NFW-normalization
figure in the existing analysis. It preserves the likelihood, templates, masks,
response, fit freedoms and manuscript text. The supplied
`Totani_two_figures_Codex_brief.md` was the task brief; this checkout had no
`CODEX_BRIEF.md` or adopted `CODEX_BRIEF_ARCHITECT_POLICY.md`.

The checkout inspected was `b177828`, remote `kabazajian/totani`, branch
`figures/template-and-regional-consistency`. The current editorial source is
`Totani_Halo_v2_copyedited.tex`, supplied separately. The repository article is
older; only focused insertion snippets are supplied. The article and bibliography
are unchanged. The unavailable sibling `gce-fisher/CONVENTIONS.md` and `CLAUDE.md`
could not be read here; owner exports use the vendored loaders and their actual
coordinate rules. No LFS/submodule assets were substituted for numerical arrays.

## Run on Kev's machine

Use the existing analysis Python (the README's gce-fisher `.venv`, or the adopted
environment with NumPy, SciPy, Astropy and Matplotlib >=3.6). Python must be >=3.11.
No fermitools command, photon reduction, new GALPROP archive, GPU, cluster, or
new simulation ensemble is required. If installing into an existing suitable
environment, `python -m pip install -e '.[paper-figures]'` supplies the plotting
extra; the renderer alone needs only NumPy and Matplotlib.

From this review branch, one command completes export, checks and rendering:

```bash
bash scripts/run_paper_figures.sh \
  --data-root "$TOTANI_DATA_ROOT" \
  --gce-fisher-root "$GCE_FISHER_ROOT" \
  --erosita-root "$EROSITA_DATA_ROOT" \
  --fermi-data-root "$FERMI_DATA" \
  --allow-nfw-constructor --recover-profiles --workers 1
```

Set these variables to existing directories first. `TOTANI_DATA_ROOT` is the
owner analysis project root containing `runs/` and `results/fits/ranking_cells/`,
not a photon archive. `GCE_FISHER_ROOT` contains the flat bubble delivery under
`data/external/`. `EROSITA_DATA_ROOT` contains the normalized full-sky FITS, or its
`eROSITA_fullsky_lat_handoff/` subdirectory. `FERMI_DATA` contains `catalog/`.
Set `PYTHON` to the analysis interpreter if `python3` is not that interpreter.
When this branch is checked out in Kev's existing analysis tree, `--data-root .`
is sufficient for the runs and caches. Paths also accept CLI flags or JSON
configuration; no Python source edits are needed. No machine-specific paths or
credentials are committed.

**The two explicit opt-ins matter.** `--allow-nfw-constructor` calls the existing
adopted J-map constructor if no verified intensity FITS has been supplied. It
uses the original 960x960 geometry, parameters, integration and chunking, and
does not build a response or fit. `--recover-profiles` reproduces the existing
regional fits because their exporter discarded the regional coefficients and
profiles. It performs 54 four-amplitude bin fits and 54 common-bin fits, followed
by 24 one-dimensional joint profiles at bin 5. Each profile includes the other
three regional amplitudes and all foregrounds. Several bracket/root evaluations
per profile are needed; runtime depends on the owner CPU. All fit bounds and
optimizer options come from `ns_halo_spectra.joint`. The recovery must reproduce
the saved common coefficients and likelihood differences (within 0.05 in
minus twice log likelihood) before it can be used. A disagreement stops that
export and must be investigated. It never starts a calibration ensemble.

Existing caches are read only. The wrapper caps BLAS/OpenMP threads at one;
`--workers` controls only the independent profile tasks, with no nested pool.
Recovered profiles are written under `results/paper_figures/`, with input identities,
solver identity and source revision. Subsequent runs omit both opt-ins when a
verified NFW intensity FITS is configured and saved profiles exist. Profiles can
be reused without a new recovery; omit `--recover-profiles`. A default `prepare`
never fits or regenerates J.

The wrapper processes the two figures independently. Use `--figure templates`
or `--figure regional` to finish one while the other lacks inputs. Missing assets
produce exit status 2 and precise requirements; they never produce mock manuscript
figures. The underlying commands can also be run separately:

```bash
python scripts/paper_figures.py inspect --config configs/paper_figures.json
python scripts/paper_figures.py prepare --config configs/paper_figures.json
python scripts/paper_figures.py check --config configs/paper_figures.json
python scripts/paper_figures.py render --config configs/paper_figures.json
```

`prepare` exports; `render` reads only the bundle and makes no network request.
`check` validates scientific definitions and bundle consistency, rather than
calling a likelihood. A rendering machine can use this checkout plus the compact
bundle without LAT/GALPROP/eROSITA deliveries or the analysis environment.

## Required inputs and current availability

| Asset and unit | Owner location / producer | Available here |
|---|---|---|
| Actual Galactic CAR WCS, degrees | `data_root/runs/fff2de5bc7e0/ccube.fits` | Missing |
| NFW input J, inherited `GeV^2 cm^-5 sr^-1` convention | `totani.halo.nfw_rho2_on_target(actual_WCS)`, explicit opt-in; optional verified FITS via `nfw_intensity_path` in JSON | Constructor available; map missing |
| Flat Fermi bubbles 0/1 intensity shape | `gce_fisher_root/data/external/FermiBubblesPublicFiles/bubbleTemp_wcs.fits` | Missing |
| Structured Fermi bubbles input, arbitrary delivered units | `external/FermiBubbles_structured_fig12_inpainted_ALLSKY_Macias.fits` | Present, real partial export checked |
| eROSITA N/S brightness weights, each originally unit solid-angle integral | `erosita_root/erosita_brightness_fullsky_normalized.fits`, extensions `T_N`, `T_S`, `U_N`, `U_S` | Missing |
| Extended-source exclusion geometry, degrees | `fermi_data_root/catalog/gll_psc_v35.fit`, `ExtendedSources` | Missing |
| Common outside-bubbles mask and folded J unit counts | `data_root/results/fits/ranking_cells/cells2/library.npz`, `cells2_psfuc/library.npz`, `cells1_psfuc/library.npz` | Arrays missing; saved summary present |
| Final common coefficients and regional LR, inherited coefficient / `-2 log L` | `results/fits/quadrant_calibration_cells2_psfuc_2026-09-29.json` | Present |
| Four-amplitude estimates and 68% profiles | `results/paper_figures/recovered_profiles.json`; explicit original-solver recovery | Missing; not saved by old exporter |
| Original pre-PSF regional response cubes | `data_root/results/fits/ranking_cells/cells2_psfuc/halo_quadrants.npz` | Missing |
| Baseline gas1 diffuse response cubes | Same cache, `row__SL_Z6_R30_T150_C2__gas1_ics1.npz` | Missing |
| Twelve brightness pieces | Same cache, `erosita_flex.npz` | Missing |
| ULTRACLEAN PSF reference | `external/psf/`, original run `psf_P8R3_ULTRACLEAN_V3_evtype3.fits` | Tracked reference present; cached response required |
| Calibration, 200 independent skies per configuration | Active final quadrant JSON above | Present; SO reference diffuse row |
| Existing 1-degree comparison | `results/fits/quadrant_calibration_cells1_psfuc_2026-09-29.json` | Present; included as diagnostics |

Use the existing final caches from Kev's analysis; do not rebuild a full analysis
just to export figures. `quadrant_calibration.build_quadrants` originally produced
the quadrant cache and `erosita_flex_ns.build_flex` produced the flexible pieces.
If these owner caches have been lost, this bounded handoff reports their absence
instead of silently rerunning their parent analysis. `inspect` checks actual files,
including whether a supposed FITS asset is an LFS pointer. A verified optional J
FITS must match the actual WCS and declare the adopted constructor's parameter
and integration provenance; the adapter reports the required header keys.

## Scientific conventions and saved-result checks

The active `_psfuc` files were rerun after the October 1 exact-position
point-source rebuild, despite retaining September filenames. This is documented
in `validation/STATE_2026-09-28.md` sections 0c/0d and
`results/fits/psfuc_chain3_2026-10-01.sh`. The export rejects `pixelsrc` archive paths,
CLEAN-proxy and response-free selections. It does not choose a result by date
or maximize a displayed statistic. Other existing diffuse rows and 1-degree
comparisons remain diagnostics in the bundle.

The regional basis uses `H[e,r]=R_e[J*m_r]`, and its sum is checked against full
J on **every response cell before the fitted-cell cut**. Source masks act in
the same order as the analysis. The basis includes PSF spillover, so observed
quadrant photon counts are unsuitable for a diagonal count-to-amplitude
conversion. The original solver normalizes each signal by its absolute count
sum and already returns `alpha/c`; the adapter does not divide twice. Points,
profile roots and saved scans share this inherited physical normalization.
The common line is directly fitted to the same data. It is never a regional
average or a per-panel denominator. Estimates are signed on one linear scale.

Standard labels are `WN`, `WS`, `EN`, `ES`, with W=`l<0` and E=`l>=0`.
The rotated partition is exactly N: `b>|l|`, S: `b<-|l|`,
E: `l>=0 and l>=|b|`, W: `l<0 and l<=-|b|`; diagonal ties go to E/W.
The compact schematics depict these exact coordinate formulas. They show the
plane exclusion but omit extended-source holes, which Figure A shows. They
are labeled as partition schematics, not substitutes for the original pixel masks.

The points illustrate bin 5 (actual center about 20.8 GeV); `Q6` uses bins 3–8
over approximately 5.6–130 GeV. The bundle stores the actual repository energy
edges, nine-bin common coefficients, likelihoods, per-bin statistics, all
four-amplitude coefficients, and illustrated-bin profile scans. No `Q+` or
maximum-over-partitions probability is substituted for a fixed-partition test.

| Baseline configuration | Standard `Q6`, exceedances | Rotated `Q6`, exceedances |
|---|---|---|
| Structured + uniform eROSITA N/S | 155.9, 0/200 | 223.1, 0/200 |
| Flat + uniform eROSITA N/S | 161.8, 0/200 | 208.8, 0/200 |
| Flat + twelve brightness pieces | 15.7, 107/200 | 174.6, 0/200 |

These were read from the active saved outputs and agree with the current
manuscript at its quoted precision. The manuscript's `p=0.5–0.9` is a model-set
range: baseline raw exceedances are 107/200=0.535, while the old add-one
estimator is 108/201=0.537. The figure reports the actual count. All existing
null simulations use `SO_Z6_R30_T150_C5 gas1_ics1`, whereas the displayed data use
`SL_Z6_R30_T150_C2 gas1_ics1`. Annotations/caption explicitly identify this
**reference null**, not an exactly baseline-matched probability. Null arrays,
unique sky identities and finite statistics are checked before counting.
Zero exceedances are `0/200`, not evidence for an extreme Gaussian tail; if a
bound is discussed, the one-sided finite-simulation 95% bound is about 1.5%.

The common outside mask is chosen on the unconvolved 2-degree flat support,
then carried to the four contained 1-degree cells. It is not reselected at fine
resolution. Saved summaries give 2190 coarse and 8760 fine cells, with identical
folded J count weight `0.5247083836046436` (rounded 52.5%). This is **response-folded
template count weight at bin 5, not sky area**. Real array equality and WCS
registration still require owner caches and are checked during preparation.

The eROSITA delivery's N/S weights are independently normalized by the original
analysis. The assembled panel therefore states **template weights**, rather than
calibrated relative X-ray brightness. No display normalization is applied per
hemisphere or native/reconstructed part. The provenance masks use delivery
metadata plus the documented analysis East-reconstruction rule; they are not
inferred from screen position.

## Bundle and figure products

Outputs, profiles and proofs under `results/paper_figures/` are ignored by Git.
Only code/configuration/tests/docs belong on this review branch. Share the compact
bundle privately with Oscar if he is to render it elsewhere; keep unreleased
bundles and generated figures out of the remote.

The exchange format is schema version 1:

| File | Contents |
|---|---|
| `manifest.json` | Scientific/test provenance, source revision, per-figure status, relative file paths and SHA256 identities |
| `spatial.json` | Actual WCS, source identities, intensity units, original and whole-map display scales, eROSITA provenance, common-mask counts and folded-count definition |
| `spatial.npz` | Numeric maps, coordinate edges, Boolean observation/exclusion/support/provenance masks, original coarse and contained fine common masks; no pickle |
| `regional.json` | Both exact partition definitions and schematics, energy edges, common and four physical amplitudes, joint profile intervals/scans, likelihoods, nuisance freedoms/bounds, basis closure, calibration null values and counts, existing diagnostics |
| `render_manifest.json` in figure directory | Product identities associated with each input bundle figure; production proof refuses stale figures |

Every displayed point, color scale and caption statistic is recoverable from
these files. `null` means absent, never measured zero. Production rejects incomplete
figures and test provenance. `np.load(..., allow_pickle=False)` is used throughout.
The default render produces `template_overview.pdf` and
`regional_nfw_normalizations.pdf` with 600-dpi image layers and vector text,
boundaries and error bars; each also has a 300-dpi PNG and grayscale PNG.

Figures use the actual REVTeX width, 510 TeX points (7.05687 inches), all labels
10 points, white background, LaTeX text in Computer Modern Roman to match the
REVTeX manuscript (a TeX installation must be on PATH), restrained colors and
equal map aspect ratios. No tight
bounding box rescales the placed figure. Display uses one recorded max-absolute
scale per whole template, with the same rule for flat/structured maps. Original
negative values remain; excluded, zero and missing pixels have distinct encodings.
Boundaries use actual pixel edges without smoothing or interpolated outlines.

## Manuscript insertion and proof

Use `docs/paper_figure_insertions.tex`: block A at the end of Section II, block B
near `sec:quadrantfinal` in Section V. It contains the two `figure*` blocks,
captions using existing citation keys and one introductory sentence per figure.
Keep Table I, existing figures, approved title/abstract/prose/results and separate
`refs.bib` in place. The proof helper copies the supplied current manuscript to a
separate generated directory; it never writes the original.

Append this to the owner wrapper command for a proof in Kev's existing REVTeX
environment:

```bash
--manuscript /path/to/current/Totani_Halo_v2_copyedited.tex
```

Or after rendering:

```bash
python scripts/paper_figure_proof.py \
  --manuscript /path/to/current/Totani_Halo_v2_copyedited.tex
```

If the editorial anchors have changed, use the focused snippets manually; the
helper rejects an ambiguous insertion. Review the two figure pages at full placed
size and in grayscale, including captions and float placement. The automatically
written `proof_report.json` records input identity, preservation of approved text
and compiler warnings; it does not claim human visual inspection.

The four presentation references in the brief were viewed at the requested PDF
pages: [Bartels et al., Figure 1](https://arxiv.org/pdf/1711.04778),
[Calore et al., Figures 15/16](https://arxiv.org/pdf/1409.0042),
[Leane and Slatyer, Figures 1/2](https://arxiv.org/pdf/2002.12370), and
[Su et al., Figure 4](https://arxiv.org/pdf/1005.5480).
Only matched framing, direct region labels and restrained overlays informed the
design; none of their measurements is imported.

## Verification here and remaining owner work

All 44 focused tests passed locally. They cover WCS handedness, map/mask registration, exact common-region
coverage, response-basis closure including leakage and unequal weights, physical
coefficient conversion, genuine signed profile endpoints, original joint-solver
API, fixed-partition statistic/calibration selection, and rejection of test data.
Tests use very small synthetic inputs and never run the real likelihood.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q \
  tests/test_paper_figures_bundle.py tests/test_paper_figures_spatial.py \
  tests/test_paper_figures_regional.py tests/test_paper_figure_proof.py
```

Layout fixtures are deliberately separate and visibly marked
`TEST DATA - NOT SCIENTIFIC RESULTS`:

```bash
python scripts/paper_figures.py fixture
python scripts/paper_figures.py check --test-data
python scripts/paper_figures.py render --test-data
```

Their directory is `results/figure_tests/`, never the production figure directory.
Production refusal and default no-refit behavior were exercised here. Real saved
statistics/common-mask summaries and the delivered structured map were checked;
full masks, regional amplitudes and real scientific figures remain owner work.
The actual current editorial manuscript was compiled with watermarked fixtures
and inspected at final size and in grayscale. The revised layouts fit the page
with captions, embedded vector fonts and 600-dpi raster layers. This proves layout,
not the missing scientific inputs. Once real exports arrive, repeat the production
checks and visual manuscript proof because real dynamic ranges may differ.

| Figure | Inputs here | Adapter tested | Real render | Scientific checks | Manuscript proof | Review branch |
|---|---|---|---|---|---|---|
| Spatial templates | Structured only; remaining inputs listed above missing | Real partial export + synthetic full export | Pending owner | Synthetic geometry/conversion; saved common summary verified | Watermarked layout only | Code/docs on feature branch |
| Regional NFW | Final common/LR/calibration available; caches/four profiles missing | Real summary export + actual solver on tiny fixture | Pending owner | Final statistics verified; physical/profile tests pass | Watermarked layout only | Code/docs on feature branch |

Message Oscar can send to Kev (not sent automatically):

> The review branch adds the two requested figure exporters/renderers, scientific
> checks and insertion snippets without changing the manuscript or fits. Please
> run the command in `docs/paper_figures.md` in your existing analysis environment.
> The explicit recovery option retrieves the discarded regional coefficients and
> joint profile intervals using the same solver and checks them against your final
> outputs. The local proofs use labeled test data; real figures need your cached
> inputs. The calibration labels disclose that the existing null ensemble uses
> the SO reference model rather than the plotted SL baseline. Please review the
> final PDFs at manuscript size and share the compact bundle privately.
