# totani — foreground dependence of the 20 GeV Galactic halo gamma-ray excess

Analysis code for

> K. N. Abazajian, Y. Li, and O. Macias, *Foreground dependence of the 20 GeV
> Galactic halo gamma-ray excess and its proposed connection to the
> LUX-ZEPLIN event* (2026), arXiv:[to be added].

Totani (arXiv:2507.07209) reported a halo-like excess peaking near 20 GeV over
|ℓ| ≤ 60°, 10° ≤ |b| ≤ 60° and fitted it with dark-matter annihilation. This
repository re-measures that excess with fifteen years of Fermi-LAT data, first
under the original foreground assumptions and 10° binning, then on 2° and 1°
cells with the instrument response, against 281 Galactic diffuse-emission
configurations, three descriptions of the Fermi bubbles, four of Loop I, and
uniform and X-ray-brightness-weighted eROSITA-bubble templates, with simulated
skies to calibrate what foreground mismodeling can produce. Every number,
table and figure in the paper is produced by a script here; the map from
paper item to script and output file is [docs/RUNBOOK.md](docs/RUNBOOK.md).

## Layout

| | |
|---|---|
| `src/totani/` | the library: data reduction onto the 0.125° grid, templates, cell likelihood, selection, PSF folding, simulated skies, figure export |
| `src/totani/_vendor/` | reduction and template helpers vendored from the authors' `gce-fisher` code, byte-identical and sha256-pinned (`_vendor/UPSTREAM.md`) |
| `scripts/` | one driver per analysis step and per paper item |
| `tests/` | `python -m pytest tests/` (units, grids, likelihood, statistics, vendored-code hashes) |
| `results/` | the tables (`.md`) and fit outputs (`results/fits/*.json`) that the paper quotes, and the chains (`results/fits/*.sh`) that produced them |
| `docs/RUNBOOK.md` | every step from public inputs to each paper item |
| `docs/paper_figures.md` | the two manuscript figures |
| `external/` | READMEs giving the URL and checksum of each external input |

## Installation

Python ≥ 3.11 with NumPy, SciPy and Astropy:

```bash
pip install -e .                   # library
pip install -e '.[paper-figures]'  # + matplotlib, for the figures
pip install pytest && python -m pytest tests/
```

Fermitools 2.5.3 is needed only to rebuild the counts and exposure cubes from
the photon files (runbook step 1); the Zenodo data record below lets you skip
that step.

## Inputs

Public inputs, each with its URL and checksum in the README of its directory:

| input | source | location |
|---|---|---|
| LAT photon and spacecraft files, 4FGL-DR4, IRFs | [Fermi Science Support Center](https://fermi.gsfc.nasa.gov/ssc/) | runbook step 1 |
| 2012 GALPROP grid of Ackermann et al. (2012) | [galprop.stanford.edu/PaperIISuppMaterial](https://galprop.stanford.edu/PaperIISuppMaterial/) | `external/galprop_paperII/` |
| 2014 flat Fermi-bubbles template | FSSC public files (Ackermann et al. 2014) | `$GCE_FISHER_ROOT/data/external/FermiBubblesPublicFiles/` |
| Haslam 408 MHz (Remazeilles et al. 2015) | [LAMBDA](https://lambda.gsfc.nasa.gov/) | `external/haslam/` |
| PPPC4DMID spectra (Cirelli et al. 2011) | [marcocirelli.net/PPPC4DMID](http://www.marcocirelli.net/PPPC4DMID.html) | `$GCE_FISHER_ROOT/data/external/PPPC4DMID/` |
| Cholis, Hu & Zhong (2026) bubble templates | [Zenodo 23199004](https://zenodo.org/records/23199004) | `external/cholis2026/` |
| LAT Galactic model `gll_iem_v07` (reference only) | FSSC | `external/fssc/` |

`GCE_FISHER_ROOT` is the data root for the templates shared with the authors'
`gce-fisher` analysis (default: a `gce-fisher` directory beside this one); its
expected layout is given in the runbook.

Available from the authors on request:

- the Galp21 configurations: six GALPROP v57 models computed for this work
  (GALPROP, https://galprop.stanford.edu/; Porter, Jóhannesson & Moskalenko
  2022). GALPROP outputs are not redistributed here, and the pipeline runs
  without them: every script that lists Galp21 rows skips them when their
  cubes are absent (`totani.ranking.available`);
- the A2020 diffuse templates of Abazajian et al. (2020) and the solar
  inverse-Compton template (Pohl et al. 2022);
- the structured Fermi-bubbles map and the eROSITA-bubbles brightness
  templates (O. Macias);
- the HESTIA particle data and halo profiles (Muru et al. 2025, from their
  authors).

## Analysis-ready inputs (Zenodo)

The Zenodo data record [DOI to be added] holds the counts cube, exposure, PSF
table and configuration of the data selection, the 10° cell caches of the
structure templates, Totani's residual-bubble constructions, the 2° and 1°
template libraries, and the simulation truths. Unpack it into the repository
root and start at runbook step 3. It contains no GALPROP-derived file: the
diffuse-model caches are rebuilt from the public 2012 grid by runbook step 2.

## Citation and acknowledgment

Please cite the paper above (`CITATION.cff`). The analysis uses GALPROP
(https://galprop.stanford.edu/; Ackermann et al. 2012; Porter, Jóhannesson &
Moskalenko 2017, 2022; Jóhannesson, Porter & Moskalenko 2018), Fermi-LAT data
and tools, and the external inputs listed above; cite them as their READMEs
state.

## License

MIT (see `LICENSE`).
