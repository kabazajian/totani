# Does the halo's S/N at 20.8 GeV follow the eROSITA lobes? PSF-folded 2° cells, flat mask (2026-09-27)

Halo fitted as four quadrant pieces (W = l<0, native eROSITA; E = l>0, reconstructed). Medians over rows (data: 33-row subset + Totani's baseline; reference: 20 flat-bubble skies of the true row with a spherical halo, unif lobes true). Δ(−2lnL) H0 is the 9-bin background-only fit relative to the unif default (negative = better; 92 per extra template is the BIC cost). Quadrant LR: common amplitude vs 4 quadrants at 20.8 GeV, χ²₃. Lobe S/N and E/W: fitted lobe counts at 20.8 GeV.

## Data

| eROSITA model | lobe templates | Δ(−2lnL) H0 vs unif | halo total | S/N | S/N West | S/N East | quadrant LR p | lobe counts N / S | lobe E/W |
|---|---|---|---|---|---|---|---|---|---|
| none | 0 | +3061 | 6,079 | 1.12 | 0.94 | 1.35 | 0.0015 | 0 / 0 | nan |
| unif | 2 | +0 | 5,750 | 1.37 | 1.26 | 1.50 | 0.00047 | 1,241 / 280 | 0.75 |
| bright | 2 | +959 | 5,045 | 1.23 | 1.03 | 1.47 | 0.0018 | 989 / 607 | 0.78 |
| unif_outer | 2 | -78 | 5,773 | 1.27 | 1.17 | 1.37 | 0.0018 | 992 / 158 | 0.75 |
| bright_outer | 2 | +1039 | 5,661 | 1.24 | 1.07 | 1.43 | 0.0014 | 807 / 256 | 0.82 |
| unif_EW | 4 | -178 | 5,759 | 1.36 | 1.20 | 1.51 | 1.5e-05 | 1,219 / 403 | 0.28 |
| bright_EW | 4 | +732 | 5,232 | 1.35 | 1.08 | 1.60 | 0.014 | 1,064 / 480 | 0.26 |
| unif_flex | 12 | -1196 | 5,678 | 1.29 | 1.07 | 1.54 | 0.0063 | 1,303 / 587 | 0.50 |
| bright_flex | 12 | -1278 | 4,860 | 1.09 | 0.98 | 1.21 | 0.67 | 1,001 / 925 | 0.71 |
| unif_mirrorW | 2 | +1484 | 5,800 | 1.27 | 1.05 | 1.56 | 0.00076 | 967 / 615 | 1.00 |
| bright_mirrorW | 2 | +1792 | 5,176 | 1.11 | 0.88 | 1.40 | 0.0013 | 840 / 846 | 0.96 |
| unif_Wonly | 2 | +1153 | 5,973 | 1.28 | 1.21 | 1.33 | 2.2e-06 | 775 / 368 | 0.00 |
| bright_Wonly | 2 | +2370 | 5,666 | 1.23 | 1.06 | 1.37 | 1.4e-05 | 632 / 444 | 0.00 |

## Reference skies (spherical halo injected)

| eROSITA model | halo total | S/N [16–84%] | S/N West [16–84%] | S/N East [16–84%] |
|---|---|---|---|---|
| none | 6,940 | 1.20 [1.06, 1.25] | 1.24 [1.11, 1.36] | 1.08 [1.03, 1.18] |
| unif | 6,440 | 1.11 [0.99, 1.17] | 1.17 [1.04, 1.31] | 0.97 [0.89, 1.08] |
| bright | 5,732 | 1.00 [0.91, 1.08] | 1.07 [0.86, 1.16] | 0.94 [0.84, 1.05] |
| unif_outer | 6,371 | 1.13 [1.00, 1.18] | 1.20 [1.09, 1.36] | 0.99 [0.92, 1.10] |
| bright_outer | 6,320 | 1.14 [1.00, 1.21] | 1.18 [1.03, 1.31] | 1.04 [0.98, 1.16] |
| unif_EW | 6,451 | 1.11 [0.99, 1.16] | 1.13 [1.03, 1.29] | 1.03 [0.88, 1.11] |
| bright_EW | 5,634 | 1.01 [0.85, 1.08] | 1.05 [0.95, 1.26] | 0.88 [0.75, 1.01] |
| unif_flex | 6,445 | 1.09 [0.98, 1.21] | 1.08 [1.00, 1.32] | 1.06 [0.96, 1.19] |
| bright_flex | 6,240 | 1.10 [0.98, 1.18] | 1.12 [0.99, 1.35] | 1.03 [0.84, 1.15] |
| unif_mirrorW | 6,250 | 1.03 [0.98, 1.13] | 1.11 [1.02, 1.17] | 0.98 [0.89, 1.08] |
| bright_mirrorW | 5,941 | 0.96 [0.86, 1.01] | 0.99 [0.89, 1.07] | 0.89 [0.82, 1.00] |
| unif_Wonly | 6,590 | 1.15 [1.04, 1.20] | 1.12 [1.02, 1.31] | 1.11 [1.04, 1.19] |
| bright_Wonly | 6,213 | 1.07 [1.02, 1.19] | 1.03 [0.93, 1.23] | 1.11 [1.04, 1.20] |

## Totani's baseline row, quadrant halo counts (WN, WS, EN, ES)

| eROSITA model | WN | WS | EN | ES |
|---|---|---|---|---|
| none | 1,435 | 1,357 | 1,012 | 1,561 |
| unif | 1,049 | 1,367 | 864 | 1,547 |
| bright | 1,103 | 1,138 | 821 | 1,385 |
| unif_outer | 1,094 | 1,322 | 922 | 1,497 |
| bright_outer | 1,179 | 1,286 | 903 | 1,499 |
| unif_EW | 936 | 1,205 | 981 | 1,689 |
| bright_EW | 957 | 1,044 | 907 | 1,652 |
| unif_flex | 983 | 1,119 | 905 | 1,584 |
| bright_flex | 1,046 | 1,018 | 872 | 1,161 |
| unif_mirrorW | 1,203 | 1,310 | 777 | 1,496 |
| bright_mirrorW | 1,244 | 1,114 | 823 | 1,338 |
| unif_Wonly | 946 | 1,234 | 1,131 | 1,702 |
| bright_Wonly | 995 | 1,066 | 1,071 | 1,663 |
