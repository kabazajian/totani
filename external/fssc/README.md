# FSSC downloads (gitignored)

`gll_iem_v07.fits` — the LAT team's Galactic interstellar emission model for
4FGL, downloaded 2026-09-19 from
https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/gll_iem_v07.fits
All-sky GLON-CAR 0.125 deg (2880 x 1441), 28 log-energy planes 50 MeV - 814 GeV,
dN/dE in MeV^-1 cm^-2 s^-1 sr^-1. CRPIX2 is an integer (721.0), so its
latitude pixel centres are offset half a pixel from ours; `totani.structures.gll_iem_v07`
averages the two bracketing rows. Reference row only in the ranking (Kev, 2026-09-19).
