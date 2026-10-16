# PSF class check: the table used vs ULTRACLEAN (2026-09-30)

The 2°/1° PSF folds (`src/totani/psf.py`) use gce-fisher's gtpsf table `runs/1bbcd8ac5ddb/psf.fits`:
P8R3_CLEAN_V3, evtype 3, Galactic Center direction, built with that run's livetime cube
(weeks 9–741, zenith < 90°). The data are P8R3_ULTRACLEAN_V3, evtype 3, weeks 9–791, zenith < 100°.
gtpsf was rerun with this analysis's own livetime cube (`runs/_ltcubes/ltcube_788ab00a4d3b.fits`)
for both classes; containment radii per analysis bin (E^-2.4 × exposure weighted, as in psf.py):

| bin (GeV) | R68 used | R68 CLEAN, our ltcube | R68 ULTRACLEAN | R95 used | R95 CLEAN | R95 ULTRACLEAN | ULTRACLEAN/used R68, R95 |
|---|---|---|---|---|---|---|---|
| 1.5 | 0.610 | 0.591 | 0.570 | 2.077 | 2.038 | 1.884 | 0.934, 0.907 |
| 2.6 | 0.417 | 0.395 | 0.383 | 1.584 | 1.519 | 1.405 | 0.917, 0.887 |
| 4.3 | 0.271 | 0.270 | 0.264 | 1.126 | 1.131 | 1.055 | 0.976, 0.937 |
| 7.3 | 0.197 | 0.193 | 0.189 | 0.858 | 0.847 | 0.801 | 0.960, 0.934 |
| 12.3 | 0.143 | 0.146 | 0.144 | 0.647 | 0.664 | 0.634 | 1.007, 0.980 |
| 20.8 | 0.118 | 0.119 | 0.117 | 0.555 | 0.562 | 0.537 | 0.991, 0.968 |
| 35.1 | 0.104 | 0.106 | 0.105 | 0.495 | 0.508 | 0.485 | 1.005, 0.979 |
| 59.2 | 0.100 | 0.101 | 0.100 | 0.455 | 0.464 | 0.445 | 0.999, 0.977 |
| 100 | 0.099 | 0.100 | 0.099 | 0.417 | 0.428 | 0.413 | 0.998, 0.991 |

Reading: in the response-valid bins (5.6–130 GeV) the ULTRACLEAN 68% radius is within 4% of the one
used and the 95% radius within 7% (the used kernel slightly wider). At 1.5–4.3 GeV the used kernel is
3–9% wider at 68% and 6–13% wider at 95%. Half of the low-energy difference is the livetime cube
(CLEAN with our cube is already 3–5% narrower at 1.5–2.6 GeV), half the event class.
