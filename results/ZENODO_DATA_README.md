# Analysis-ready inputs for "Foreground dependence of the 20 GeV Galactic halo gamma-ray excess and its proposed connection to the LUX-ZEPLIN event"

K. N. Abazajian, Y. Li, and O. Macias (2026), arXiv:[to be added].
Analysis code: https://github.com/kabazajian/totani (software record: [DOI to be added]).

This record holds the inputs that the analysis code reads once the Fermi-LAT
photon reduction is done, so that the fits, selections, simulations, tables
and figures of the paper can be rerun without Fermitools or the photon files.
Unpack the three archives into the root of a checkout of the code and follow
`docs/RUNBOOK.md` from step 3.

| archive | contents | size |
|---|---|---|
| `totani_runs.tar` | `runs/fff2de5bc7e0/`: the binned counts cube (`ccube.fits`; 960 × 960 pixels of 0.125° in plate carrée Galactic coordinates centred on the Galactic Center, 13 energy bins from 1.16 GeV to 1.06 TeV, of which the nine up to 130 GeV are fitted), the exposure cube (`expcube.fits`), the energy bounds, the data-selection configuration (`config.json`: Pass 8 ULTRACLEAN front+back, weeks 9–791, 2008-08-04 to 2023-08-03, zenith < 100°, standard good-time filter), and the exposure-weighted point-spread function table (`psf_P8R3_ULTRACLEAN_V3_evtype3.fits` from `gtpsf`, with its three-Gaussian fit in `.json`) | 112 MB |
| `totani_ranking_cells.tar` | `results/fits/ranking_cells/`: the 10° cell caches of the shared and structure templates (`shared__base.npz`, `structure__*.npz`: counts, point sources, isotropic, Fermi-bubble, Loop I, eROSITA-bubble, solar inverse-Compton and reference-model templates summed on Totani's 10° cells); Totani's residual-bubble constructions for each 2012-grid configuration (`totani_residual__*.npz`); the 2° and 1° template libraries with and without the instrument response (`cells2/`, `cells1/`, `cells2_psfuc/`, `cells1_psfuc/`: `library.npz`, the quadrant pieces of the NFW template, the eROSITA pieces, the Cholis et al. components, and the simulation truths); and the 10° simulation truths (`stage_e_truth__*.npz`) | 256 MB |
| `totani_stage_e.tar` | `results/fits/stage_e/`: the 10° simulated skies and their fits | 37 MB |

`ARCHIVES.sha256` gives the SHA-256 of each archive and `CONTENTS.sha256` that
of every file inside them.

## What is not here, and why

- **No GALPROP-derived file.** The per-model diffuse caches (`model__*.npz` at
  10°, `row__*.npz` at 2° and 1°, and the pixel-level simulation truths
  `stage_e_pix2/`) are rebuilt from the public 2012 GALPROP grid of Ackermann
  et al. (2012) by runbook step 2 of the code; the Galp21 configurations
  (six GALPROP v57 models computed for this work with GALPROP,
  https://galprop.stanford.edu/, Porter, Jóhannesson & Moskalenko 2022) are
  available from the authors on request and are not redistributed.
- **The photon and spacecraft files** (Fermi Science Support Center) and the
  **4FGL-DR4 catalog** are public and large; only their reduction is here.
- **The superseded CLEAN-proxy response** libraries (`cells*_psf/`) are
  omitted; every number in the paper uses the ULTRACLEAN response
  (`*_psfuc`).

## Provenance of the templates in the libraries

The structure templates and libraries were built by the released code from:
the 2014 flat Fermi-bubbles template (Ackermann et al. 2014, FSSC public
files); the structured Fermi-bubbles map of the LAT 2017 analysis as inpainted
by O. Macias; the eROSITA-bubbles brightness and uniform templates of
O. Macias (from the public eROSITA image, Predehl et al. 2020; the eastern
half reconstructed); the Haslam 408 MHz map (Remazeilles et al. 2015) and the
Wolleben (2007) shells for Loop I; the Cholis, Hu & Zhong (2026) templates
(Zenodo 23199004, CC BY 4.0); the A2020 gas, dust and inverse-Compton
templates (Abazajian et al. 2020; Pohl et al. 2022) and their solar
inverse-Compton template; the LAT `gll_iem_v07` model (reference only); the
4FGL-DR4 catalog; and the NFW annihilation profile at Totani's parameters.
The cell caches contain these templates summed on cells, integrated over the
analysis energy bins and, in the `_psfuc` libraries, folded with the LAT
response; the templates themselves are not here.

## License and citation

CC BY 4.0. Please cite the paper above. Fermi-LAT data: NASA/DOE Fermi
mission and the Fermi Science Support Center. GALPROP:
https://galprop.stanford.edu/ (Ackermann et al. 2012; Porter, Jóhannesson &
Moskalenko 2017, 2022; Jóhannesson, Porter & Moskalenko 2018). eROSITA
image: SRG/eROSITA, MPE/IKI (Predehl et al. 2020).
