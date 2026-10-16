# Cholis, Hu & Zhong (2026) Fermi-bubble templates

Source: Zenodo record 23199004 (CC BY 4.0), ZENODO.zip, md5 a32330cae216f6f47bbb719e095c108b,
accompanying arXiv:2610.10695. Downloaded 2026-10-11. The three bubble components and the two Loop I shells used by
`scripts/cholis_bubbles.py` are kept here (the FITS files are gitignored):

| file | sha256 |
|---|---|
| Bubbles_Edges_Alternative_E_flux_140x140_E_50-814009_MeV_Cartesian.fits | 6b45628dd131e33200f0b867d3cf1b65a866e3a22fa654d85de0d8696dacf1a0 |
| BubblesInner_with_Boundaries_*.fits | 8979e0e4e637809b95ed9a39d7cb78071ad71828e7b199da2bb6a07f15f72bb4 |
| Cocoons_EmissStoN_10_*.fits | e6dfc973810525c3bc99f5933165298e7f588fa7bf440aa54533065353807cba |
| LoopI_Bright_Shell_flux_140x140_E_50-814009_MeV_Cartesian.fits | 9284661fa5005fbb6e45bc628a60e85a6ad1d9da5ef95f81084375db5c440b15 |
| LoopI_Dim_Shell_flux_140x140_E_50-814009_MeV_Cartesian.fits | 3b95aa4236f59c12b2c5249b6f46c43a18545bf8eb7ef2ac62f91173e3aeb0dc |

The files carry no WCS. Grid inferred from content: 560x560x38, FITS axis 1 = Galactic l
increasing with column, axis 2 = b, 0.25-deg pixels with centers at -69.875 ... +69.875
(+-0.125 deg ambiguity in the pixel-center convention). Evidence: the flat-bubbles footprint
matches the 2014 Fermi-LAT template best with this orientation (overlap/union 0.741 vs 0.638
mirrored), and the bright Loop I shell's limb (North Polar Spur) lies at l = +17 to +22 deg.
Values are E^2 dPhi/dE (GeV cm^-2 s^-1 sr^-1); spatial shapes are energy independent.
Southern cocoon centroid (l, b) = (+7.6, -23.0) deg; Inner Bubbles are zero at |b| < 21 deg.
