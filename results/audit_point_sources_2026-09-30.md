# Point-source audit (gll_psc_v35.fit, catalog VERSION v35; 2026-09-30)

## 1. Catalog

- rows: 7195; extended sources listed: 82
- SpectrumType values: ['LogParabola', 'PLSuperExpCutoff', 'PowerLaw']; handled: ['LogParabola', 'PLSuperExpCutoff', 'PowerLaw']; unhandled: none
- sources in |l|,|b| <= 60: 2651; with a non-finite or non-positive expected count in any bin: 0

## 2. Accounting

- kept: 2411 (in the fitted |b| >= 10 region: 1388)
- dropped because inside an extended-source mask circle: 240 (in |b| >= 10: 37)
- dropped for any other reason: 0
- masked sources at |b| >= 10: 4FGL J0024.0-7204, 4FGL J0040.7-7157, 4FGL J0046.7-7048, 4FGL J0058.0-7245e, 4FGL J0059.7-7210, 4FGL J0101.8-7543, 4FGL J0102.9-7051, 4FGL J0111.4-7304, 4FGL J0112.8-7506, 4FGL J1303.6-4622, 4FGL J1304.3-4353, 4FGL J1306.8-4035, 4FGL J1307.6-4259, 4FGL J1313.2-4437, 4FGL J1315.0-4236, 4FGL J1318.1-3937, 4FGL J1318.9-4556, 4FGL J1322.0-4624, 4FGL J1324.0-4330e, 4FGL J1325.4-4706, 4FGL J1325.5-4300, 4FGL J1326.6-4729, 4FGL J1328.5-4727, 4FGL J1330.2-3946, 4FGL J1333.2-4707, 4FGL J1334.1-4154, 4FGL J1334.6-4435, 4FGL J1334.8-3856, 4FGL J1335.5-4546, 4FGL J1335.6-4052, 4FGL J1336.5-4655, 4FGL J1336.9-4611, 4FGL J1338.0-4159, 4FGL J1338.8-4509, 4FGL J1626.9-2431e, 4FGL J1722.8-0418, 4FGL J1723.5-0501e

## 3. Cached templates against the recomputed kept-source totals (per bin, ratio)

| template | 1.51 | 2.55 | 4.31 | 7.28 | 12.3 | 20.8 | 35.1 | 59.2 | 100 |
|---|---|---|---|---|---|---|---|---|---|
| 10° shared (all 144 cells) | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 2° unconvolved | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 1° unconvolved | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| PSF fold of the kept map, flux left on the grid / before | 0.999767 | 0.999919 | 0.999921 | 0.999858 | 0.999777 | 0.999700 | 0.999630 | 0.999580 | 0.999545 |
| 2° PSF-folded library / independent fold | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 2° PSF-folded, max cell deviation from independent fold | 4.08e-07 | 7.08e-08 | 2.01e-09 | 1.55e-07 | 3.22e-08 | 7.76e-07 | 6.58e-07 | 8.50e-08 | 9.14e-08 |
| 1° PSF-folded library / independent fold | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 | 1.000000 |
| 1° PSF-folded, max cell deviation from independent fold | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |

## 4. Placement in the 10° template

- max |rebuilt − cached| / max cell, over all bins and cells: 0.00e+00

## 5. What is not modeled

- sources within 3° outside the |l|,|b| <= 60 grid (their PSF can spill in at the lowest energies): 165; at 1.5 GeV the 95% radius is 2.1°, so only the outermost cells are affected, and those are at the ROI edge, not near the halo
- fraction of the catalog flux in |b| >= 10 that is modeled (the rest sits inside the extended-source masks, whose pixels are removed from data and model alike): 0.9593, 0.9599, 0.9603, 0.9604, 0.9608, 0.9617, 0.9625, 0.9622, 0.9607

## 6. The template is in every fit

- `fine.base_stack` includes 'point_sources': True
- `selection.base_stack` (10° programme) includes 'point_sources': True
