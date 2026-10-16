# Does the halo's S/N at 20.8 GeV follow the eROSITA lobes? PSF-folded 2° cells, flat mask (2026-09-27)

Halo fitted as four quadrant pieces (W = l<0, native eROSITA; E = l>0, reconstructed). Medians over rows (data: 33-row subset + Totani's baseline; reference: 20 flat-bubble skies of the true row with a spherical halo, unif lobes true). Δ(−2lnL) H0 is the 9-bin background-only fit relative to the unif default (negative = better; 92 per extra template is the BIC cost). Quadrant LR: common amplitude vs 4 quadrants at 20.8 GeV, χ²₃. Lobe S/N and E/W: fitted lobe counts at 20.8 GeV.

## Data

| eROSITA model | lobe templates | Δ(−2lnL) H0 vs unif | halo total | S/N | S/N West | S/N East | quadrant LR p | lobe counts N / S | lobe E/W |
|---|---|---|---|---|---|---|---|---|---|
| none | 0 | +3054 | 6,085 | 1.12 | 0.94 | 1.35 | 0.0016 | 0 / 0 | nan |
| unif | 2 | +0 | 5,756 | 1.38 | 1.26 | 1.50 | 0.00045 | 1,254 / 280 | 0.75 |
| bright | 2 | +976 | 5,046 | 1.24 | 1.04 | 1.47 | 0.0018 | 1,005 / 608 | 0.78 |
| unif_outer | 2 | -77 | 5,787 | 1.27 | 1.17 | 1.37 | 0.0018 | 994 / 161 | 0.76 |
| bright_outer | 2 | +1050 | 5,663 | 1.24 | 1.08 | 1.43 | 0.0014 | 816 / 262 | 0.82 |
| unif_EW | 4 | -171 | 5,770 | 1.37 | 1.20 | 1.51 | 1.4e-05 | 1,236 / 401 | 0.30 |
| bright_EW | 4 | +756 | 5,229 | 1.35 | 1.08 | 1.61 | 0.013 | 1,077 / 482 | 0.27 |
| unif_flex | 12 | -1208 | 5,688 | 1.29 | 1.07 | 1.54 | 0.0061 | 1,313 / 586 | 0.50 |
| bright_flex | 12 | -1271 | 4,868 | 1.09 | 0.98 | 1.22 | 0.66 | 1,015 / 922 | 0.71 |
| unif_mirrorW | 2 | +1447 | 5,814 | 1.27 | 1.06 | 1.56 | 0.00079 | 976 / 615 | 1.00 |
| bright_mirrorW | 2 | +1771 | 5,182 | 1.11 | 0.88 | 1.40 | 0.0014 | 831 / 846 | 0.96 |
| unif_Wonly | 2 | +1138 | 5,984 | 1.28 | 1.22 | 1.32 | 2.3e-06 | 774 / 365 | 0.00 |
| bright_Wonly | 2 | +2370 | 5,682 | 1.23 | 1.06 | 1.37 | 1.4e-05 | 629 / 444 | 0.00 |

## Reference skies (spherical halo injected)

| eROSITA model | halo total | S/N [16–84%] | S/N West [16–84%] | S/N East [16–84%] |
|---|---|---|---|---|
| none | 7,098 | 1.13 [1.08, 1.19] | 1.20 [1.08, 1.30] | 1.08 [0.98, 1.18] |
| unif | 6,518 | 1.03 [0.93, 1.07] | 1.11 [1.00, 1.20] | 0.93 [0.86, 1.04] |
| bright | 5,763 | 0.94 [0.75, 0.99] | 0.92 [0.83, 1.08] | 0.87 [0.73, 1.02] |
| unif_outer | 6,532 | 1.06 [1.00, 1.14] | 1.16 [1.04, 1.26] | 0.98 [0.89, 1.04] |
| bright_outer | 6,478 | 1.07 [0.98, 1.13] | 1.13 [0.99, 1.20] | 1.03 [0.94, 1.11] |
| unif_EW | 6,553 | 1.04 [0.93, 1.08] | 1.06 [0.98, 1.15] | 0.95 [0.86, 1.06] |
| bright_EW | 5,625 | 0.93 [0.71, 1.00] | 0.98 [0.87, 1.15] | 0.78 [0.59, 0.95] |
| unif_flex | 6,457 | 0.96 [0.89, 1.12] | 1.09 [0.99, 1.20] | 0.92 [0.77, 1.06] |
| bright_flex | 5,981 | 0.96 [0.79, 1.15] | 1.10 [0.97, 1.23] | 0.80 [0.57, 1.02] |
| unif_mirrorW | 6,371 | 1.01 [0.91, 1.10] | 1.04 [0.95, 1.17] | 0.96 [0.83, 1.09] |
| bright_mirrorW | 5,974 | 0.89 [0.82, 0.96] | 0.91 [0.83, 1.02] | 0.85 [0.75, 0.99] |
| unif_Wonly | 6,725 | 1.09 [1.00, 1.15] | 1.06 [0.99, 1.16] | 1.09 [0.99, 1.19] |
| bright_Wonly | 6,388 | 1.08 [0.95, 1.12] | 0.95 [0.85, 1.12] | 1.10 [1.00, 1.22] |

## Totani's baseline row, quadrant halo counts (WN, WS, EN, ES)

| eROSITA model | WN | WS | EN | ES |
|---|---|---|---|---|
| none | 1,435 | 1,360 | 1,017 | 1,565 |
| unif | 1,045 | 1,371 | 867 | 1,551 |
| bright | 1,097 | 1,142 | 823 | 1,389 |
| unif_outer | 1,093 | 1,324 | 927 | 1,499 |
| bright_outer | 1,177 | 1,287 | 907 | 1,502 |
| unif_EW | 937 | 1,209 | 978 | 1,692 |
| bright_EW | 958 | 1,046 | 905 | 1,655 |
| unif_flex | 982 | 1,120 | 909 | 1,588 |
| bright_flex | 1,046 | 1,018 | 873 | 1,170 |
| unif_mirrorW | 1,201 | 1,315 | 780 | 1,500 |
| bright_mirrorW | 1,246 | 1,116 | 830 | 1,340 |
| unif_Wonly | 947 | 1,239 | 1,136 | 1,705 |
| bright_Wonly | 997 | 1,068 | 1,076 | 1,667 |
