# Validation against absolute references — the pipeline passes, after one real bug

> 🚨 **PARTIALLY SUPERSEDED by `validation/audit_2026-09-04.md` (same day).**
> The VALIDATION section stands (IGRB, point-source order-unity). The
> HALO-FIT TABLE is void: the fit included the disk cells Totani's
> procedure excludes, the alternation estimator was later shown to stall
> (injection slope 0.67 at bin 0), and 29–45% of the point-source flux
> belonged to masked sources that should not have been in the template.

**Date:** 2026-09-04. Supersedes `halo_fit_first_pass_2026-09-04.md` (same
day). Full per-bin output: `results/fits/halo_fit_v2_2026-09-04.json`.

## What was validated, and how (rule-1a compliant — no numbers off figures)

The two components with ABSOLUTE physical predictions anchor the whole chain:

**1. Point sources (the catalog is the reference — no external number).** Our
template is 4FGL-DR4 flux × exposure, so a correct chain must fit at f ≈ 1.
Result: **f = 0.65–0.96 across all 9 bins.** Order unity — excluding any gross
factor error (solid angle would be ~10⁵, keV/MeV ~10⁶, ΔΩ per-sr ~2×10⁵). The
sub-unity residual (worst at low energy) is real but modest; candidate causes
(PSF spill into masked regions, diffuse degeneracy, catalog epoch) are noted,
none of them a units error.

**2. Isotropic vs the measured IGRB** (Ackermann et al. 2015, machine-readable
table `IGRB_apj504089t3_mrt.txt` — a published measurement, anchorable).
After the solid-angle fix, our fitted isotropic amplitude is a physical
intensity, and its per-bin integrated flux against the IGRB's three
foreground models:

| bin | E (GeV) | ours (cm⁻²s⁻¹sr⁻¹) | IGRB-A | ratio vs A / B / C |
|---|---|---|---|---|
| 0 | 1.51 | 1.22e-07 | 1.37e-07 | 0.89 / 0.65 / 0.96 |
| 2 | 4.31 | 4.02e-08 | 3.04e-08 | 1.32 / 0.85 / 1.25 |
| 5 | 20.8 | 6.05e-09 | 4.19e-09 | 1.45 / 1.00 / 1.28 |
| 8 | 100 | 5.41e-10 | 3.97e-10 | 1.36 / 1.03 / 1.28 |

**0.9–1.5× the measurement across two decades of energy, sitting essentially
inside the IGRB's own A–C foreground-model envelope**, with the modest excess
over model A carrying the expected sign (our isotropic also absorbs residual
CR contamination and unmodeled smooth emission). The absolute units chain —
reduction, gtltcube/gtexpcube2, flat-spectrum bin integration, solid angle,
cell binning, Poisson fit — is validated end to end.

## The bug the validation caught

🚨 **Every diffuse template omitted the pixel solid angle.** An intensity map's
counts are value × exposure × ΔΩ(pixel), and ΔΩ ∝ cos b varies by 2× across
|b| ≤ 60° — a latitude-dependent SHAPE error no free normalization absorbs.
The tell: the fitted isotropic jumped FIVE ORDERS OF MAGNITUDE between
adjacent bins (5.1e-24 → 2.1e-18), while point sources — already integrated
fluxes, needing no ΔΩ — fit at order unity. `templates.pixel_solid_angle`
(exact CAR formula, verified against the analytic ROI total to 7 digits) is now
applied to every intensity-like template including the halo J.

After the fix: background-only halo-region residuals improved in every bin
(bin 0 mean|σ| 6.94 → 5.95; bin 5: 1.59 → 1.34), and the isotropic spectrum
became smooth and monotonic.

## The corrected halo fit

| bin | E (GeV) | halo counts | fraction | TS | SED proxy N×E |
|---|---|---|---|---|---|
| 0 | 1.51 | 49,874 | 2.7% | 249.8 | 7.5e4 |
| 1 | 2.55 | 40,205 | 4.2% | 322.4 | 1.03e5 |
| 2 | 4.31 | 25,455 | 5.5% | 279.4 | 1.10e5 |
| 3 | 7.28 | 13,401 | 6.4% | 190.8 | 9.8e4 |
| 4 | 12.29 | 7,619 | 8.0% | 152.1 | 9.4e4 |
| 5 | 20.76 | 4,188 | 9.7% | 121.4 | 8.7e4 |
| 6 | 35.06 | 1,979 | 9.7% | 62.0 | 6.9e4 |
| 7 | 59.22 | 900 | 10.1% | 28.1 | 5.3e4 |
| 8 | 100.02 | 306 | 8.4% | 7.8 | 3.1e4 |

Qualitatively unchanged from v1: a positive halo-shaped component in every
bin, SED proxy peaking at **2.5–4.3 GeV, not at 20 GeV**, TS formally large.

## Standing caveats (weakened, not lifted)

- The TS values are still not significances: halo-region residuals of
  ~6σ/cell remain at low energy, so the null model is still wrong there. The
  20 GeV bin's residuals (1.34 mean) are now close to Poisson-clean, which
  strengthens the *local* reading at the energies Totani's signal lives —
  but the systematics program (quadrant split, PS-masked, bubble variants,
  diffuse-model variants) has not been run.
- The point-source f < 1 trend with energy is unexplained at the tens-of-
  percent level.
- One diffuse family; the categorical-spread marginalization this project
  exists for still needs the |b| ≤ 60 Galp21 delivery.
