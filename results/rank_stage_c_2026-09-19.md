# Stage C — 2026-09-19

Consensus selection (median Δ(−2lnL) across 73 model rows; BIC threshold p·ln 1080 = 63 per template): **bubble_totani → loopI_shells+nps → sun_moon_ic**

| step | chosen | median Δ | min | max | threshold | runner-up (median Δ) |
|---|---|---|---|---|---|---|
| 1 | bubble_totani | 29,647 | 21,329 | 37,961 | 126 | bubble_structured (19,926) |
| 2 | loopI_shells+nps | 7,197 | 4,887 | 8,491 | 189 | loopI_shells (7,195) |
| 3 | sun_moon_ic | 79 | 7 | 139 | 63 | erosita_outer_NS (75) |
| 4 | (stop) erosita_outer_NS | 55 | -33 | 158 | 126 | loopII (38) |

## Each model's own path

| model | variant | path |
|---|---|---|
| SL_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SL_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SL_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS → loopII |
| SL_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SS_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SS_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_outer_NS |
| SS_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps |
| SS_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells |
| SS_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| a2020 | - | bubble_totani → loopI_shells+nps |

## Ranking with the consensus set (bubble_totani, loopI_shells+nps, sun_moon_ic)

| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | a2020 | - | 22 | -1,587 | 8,549 [8,056, 9,054] | 43.2 | 0.71 | no | 0.98 | 3,283 / 2,784 |
| 2 | SS_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,480 | 6,074 [5,176, 7,048] | 17.9 | 0.18 | yes | 0.97 | -4 / 1,077 |
| 3 | SS_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,475 | 6,131 [5,227, 7,037] | 18.3 | 0.22 | yes | 0.99 | 71 / 1,248 |
| 4 | SS_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,468 | 6,187 [5,271, 7,123] | 20.3 | 0.32 | yes | 0.99 | 673 / 1,301 |
| 5 | SS_Z8_R20_T150_C2 | gas5_ics3 | 16 | -1,465 | 6,187 [5,326, 7,130] | 20.6 | 0.38 | yes | 0.98 | 763 / 1,609 |
| 6 | SS_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,453 | 6,356 [5,532, 7,417] | 23.1 | 0.49 | yes | 0.99 | 923 / 2,177 |
| 7 | SS_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,447 | 6,693 [5,569, 7,473] | 23.3 | 0.48 | yes | 0.97 | 998 / 2,246 |
| 8 | SL_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,443 | 6,299 [5,340, 7,237] | 19.6 | 0.44 | yes | 0.99 | 866 / 1,880 |
| 9 | SS_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,441 | 6,580 [5,653, 7,494] | 25.6 | 0.62 | yes | 0.99 | 1,660 / 2,424 |
| 10 | SS_Z6_R30_T150_C2 | gas5_ics3 | 16 | -1,439 | 6,412 [5,515, 7,299] | 26.0 | 0.60 | yes | 0.98 | 1,860 / 1,987 |
| 11 | SS_Z6_R20_T150_C2 | gas5_ics3 | 16 | -1,438 | 6,412 [5,523, 7,307] | 26.1 | 0.60 | yes | 0.97 | 1,875 / 1,949 |
| 12 | SS_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,434 | 6,693 [5,666, 7,512] | 25.8 | 0.61 | yes | 0.98 | 1,675 / 2,414 |
| 13 | SS_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,415 | 6,637 [5,774, 7,497] | 30.0 | 0.85 | no | 0.98 | 2,965 / 2,646 |
| 14 | SS_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,413 | 6,637 [5,780, 7,506] | 30.0 | 0.86 | no | 0.99 | 2,965 / 2,773 |
| 15 | SL_Z6_R30_T150_C2 | gas5_ics3 | 16 | -1,382 | 6,580 [5,661, 7,488] | 27.5 | 0.71 | yes | 0.96 | 2,413 / 2,288 |
| 16 | SL_Z6_R20_T150_C2 | gas5_ics3 | 16 | -1,379 | 6,580 [5,667, 7,493] | 27.5 | 0.71 | yes | 1.00 | 2,383 / 2,301 |
| 17 | SS_Z4_R30_T150_C2 | gas5_ics3 | 16 | -1,376 | 6,862 [6,009, 7,534] | 41.2 | 0.85 | no | 1.00 | 3,072 / 2,766 |
| 18 | SS_Z4_R20_T150_C2 | gas5_ics3 | 16 | -1,376 | 6,862 [6,012, 7,538] | 41.2 | 0.84 | no | 0.99 | 2,968 / 2,793 |
| 19 | SS_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -1,360 | 7,087 [5,911, 7,855] | 21.0 | 0.45 | yes | 0.98 | 1,103 / 2,093 |
| 20 | SS_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -1,355 | 6,862 [5,927, 7,799] | 22.8 | 0.57 | yes | 0.98 | 1,794 / 2,101 |
| 21 | SS_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -1,352 | 7,030 [6,172, 8,038] | 27.3 | 0.61 | yes | 1.00 | 1,796 / 2,480 |
| 22 | SS_Z4_R20_T150_C5 | gas5_ics3 | 16 | -1,345 | 7,030 [6,504, 7,595] | 43.0 | 0.91 | no | 0.98 | 3,522 / 2,874 |
| 23 | SL_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,345 | 7,199 [6,461, 7,716] | 28.0 | 0.72 | no | 0.66 | 2,935 / 2,269 |
| 24 | SS_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -1,345 | 6,974 [5,999, 7,891] | 21.6 | 0.49 | yes | 0.98 | 1,178 / 2,252 |
| 25 | SS_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -1,341 | 6,974 [6,094, 7,893] | 27.5 | 0.73 | yes | 0.99 | 2,784 / 2,302 |
| 26 | SS_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -1,339 | 6,974 [6,104, 7,836] | 27.6 | 0.75 | yes | 0.99 | 2,814 / 2,393 |
| 27 | SS_Z4_R30_T150_C5 | gas5_ics3 | 16 | -1,338 | 7,030 [6,500, 7,591] | 43.0 | 0.94 | no | 0.98 | 3,537 / 3,090 |
| 28 | SS_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -1,301 | 7,143 [6,215, 8,061] | 29.1 | 0.70 | yes | 0.98 | 2,260 / 2,762 |
| 29 | SS_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -1,291 | 7,255 [6,268, 8,111] | 29.6 | 0.70 | yes | 0.99 | 2,349 / 2,717 |
| 30 | SS_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -1,287 | 7,255 [6,277, 8,173] | 28.3 | 0.62 | yes | 0.99 | 1,763 / 2,742 |
| 31 | SL_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -1,276 | 7,199 [6,328, 8,114] | 27.9 | 0.75 | no | 1.00 | 3,234 / 2,160 |
| 32 | SS_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -1,275 | 7,143 [6,306, 7,867] | 40.0 | 0.84 | no | 0.98 | 3,505 / 2,508 |
| 33 | SS_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -1,274 | 7,143 [6,310, 7,873] | 40.0 | 0.84 | no | 0.99 | 3,399 / 2,636 |
| 34 | SS_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -1,273 | 7,199 [6,307, 8,069] | 33.0 | 0.79 | no | 1.00 | 2,995 / 2,685 |
| 35 | SS_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -1,267 | 7,199 [6,318, 8,078] | 33.1 | 0.81 | no | 0.97 | 3,010 / 2,799 |
| 36 | SS_Z4_R30_T100000_C5 | gas5_ics3 | 16 | -1,209 | 7,143 [6,343, 7,871] | 41.3 | 0.85 | no | 0.99 | 3,340 / 2,702 |
| 37 | SS_Z4_R20_T100000_C5 | gas5_ics3 | 16 | -1,208 | 7,143 [6,401, 7,923] | 41.2 | 0.84 | no | 0.99 | 3,370 / 2,602 |
| 38 | SS_Z4_R30_T100000_C5 | gas1_ics1 | 10 | -340 | 6,018 [5,529, 6,514] | 155.5 | 0.74 | no | 1.00 | 1,699 / 2,751 |
| 39 | SS_Z4_R20_T100000_C5 | gas1_ics1 | 10 | -331 | 6,018 [5,553, 6,534] | 158.1 | 0.75 | no | 0.96 | 1,789 / 2,752 |
| 40 | SS_Z4_R30_T150_C5 | gas1_ics1 | 10 | -328 | 5,962 [5,461, 6,443] | 152.9 | 0.72 | yes | 1.00 | 1,355 / 2,963 |
| 41 | SS_Z4_R20_T150_C5 | gas1_ics1 | 10 | -327 | 5,962 [5,460, 6,441] | 153.0 | 0.73 | yes | 0.98 | 1,414 / 2,963 |
| 42 | SS_Z6_R30_T100000_C5 | gas1_ics1 | 10 | -289 | 6,074 [5,588, 6,568] | 160.0 | 0.74 | no | 0.94 | 1,669 / 2,848 |
| 43 | SS_Z6_R20_T100000_C5 | gas1_ics1 | 10 | -267 | 6,074 [5,592, 6,569] | 160.8 | 0.75 | no | 0.88 | 1,683 / 2,848 |
| 44 | SS_Z6_R30_T150_C5 | gas1_ics1 | 10 | -252 | 5,962 [5,472, 6,455] | 152.8 | 0.76 | no | 1.00 | 1,609 / 2,936 |
| 45 | SS_Z8_R30_T100000_C5 | gas1_ics1 | 10 | -251 | 6,131 [5,619, 6,601] | 159.9 | 0.74 | no | 0.98 | 1,698 / 2,836 |
| 46 | SS_Z6_R20_T150_C5 | gas1_ics1 | 10 | -237 | 5,962 [5,476, 6,457] | 153.5 | 0.77 | no | 1.00 | 1,624 / 2,950 |
| 47 | SS_Z4_R30_T100000_C2 | gas1_ics1 | 10 | -233 | 6,131 [5,597, 6,582] | 160.1 | 0.72 | yes | 1.00 | 1,579 / 2,827 |
| 48 | SS_Z6_R20_T100000_C2 | gas1_ics1 | 10 | -227 | 6,131 [5,669, 6,644] | 165.9 | 0.73 | yes | 0.99 | 1,683 / 2,815 |
| 49 | SS_Z4_R20_T100000_C2 | gas1_ics1 | 10 | -225 | 5,906 [5,565, 6,611] | 160.1 | 0.75 | yes | 1.00 | 1,579 / 2,827 |
| 50 | SS_Z4_R30_T150_C2 | gas1_ics1 | 10 | -214 | 5,962 [5,492, 6,474] | 158.5 | 0.74 | yes | 0.97 | 1,400 / 3,026 |
| 51 | SS_Z8_R20_T100000_C5 | gas1_ics1 | 10 | -209 | 6,131 [5,645, 6,621] | 163.1 | 0.75 | no | 1.00 | 1,743 / 2,849 |
| 52 | SS_Z4_R20_T150_C2 | gas1_ics1 | 10 | -208 | 5,962 [5,491, 6,472] | 158.4 | 0.74 | yes | 0.97 | 1,400 / 3,026 |
| 53 | SS_Z6_R30_T100000_C2 | gas1_ics1 | 10 | -195 | 6,131 [5,663, 6,640] | 164.9 | 0.73 | yes | 1.00 | 1,623 / 2,842 |
| 54 | SS_Z8_R20_T150_C5 | gas1_ics1 | 10 | -173 | 6,074 [5,571, 6,546] | 159.8 | 0.78 | no | 1.00 | 1,848 / 2,884 |
| 55 | SL_Z6_R30_T100000_C2 | gas1_ics1 | 10 | -162 | 6,524 [6,053, 7,045] | 176.5 | 0.79 | no | 0.95 | 2,372 / 2,780 |
| 56 | SS_Z10_R30_T150_C5 | gas1_ics1 | 10 | -162 | 6,074 [5,598, 6,576] | 159.9 | 0.78 | no | 0.89 | 1,908 / 2,844 |
| 57 | SS_Z8_R30_T150_C5 | gas1_ics1 | 10 | -156 | 6,018 [5,555, 6,535] | 158.0 | 0.78 | no | 0.97 | 1,803 / 2,884 |
| 58 | SS_Z10_R20_T100000_C5 | gas1_ics1 | 10 | -155 | 6,187 [5,715, 6,687] | 167.8 | 0.75 | no | 0.98 | 1,833 / 2,836 |
| 59 | SS_Z8_R30_T100000_C2 | gas1_ics1 | 10 | -154 | 6,187 [5,700, 6,679] | 165.1 | 0.76 | no | 1.00 | 1,758 / 2,926 |
| 60 | SS_Z6_R30_T150_C2 | gas1_ics1 | 10 | -150 | 6,018 [5,503, 6,486] | 155.1 | 0.76 | no | 0.94 | 1,609 / 2,987 |
| 61 | SS_Z6_R20_T150_C2 | gas1_ics1 | 10 | -129 | 6,018 [5,509, 6,490] | 155.5 | 0.75 | no | 0.94 | 1,535 / 2,986 |
| 62 | SS_Z10_R30_T100000_C2 | gas1_ics1 | 10 | -120 | 6,243 [5,778, 6,755] | 169.9 | 0.75 | no | 0.96 | 1,833 / 2,846 |
| 63 | SS_Z8_R20_T100000_C2 | gas1_ics1 | 10 | -113 | 6,187 [5,717, 6,692] | 167.3 | 0.75 | no | 1.00 | 1,698 / 2,966 |
| 64 | SS_Z8_R30_T150_C2 | gas1_ics1 | 10 | -99 | 6,074 [5,568, 6,549] | 158.1 | 0.74 | no | 0.98 | 1,551 / 2,975 |
| 65 | SL_Z6_R30_T150_C5 | gas1_ics1 | 10 | -91 | 6,243 [5,774, 6,767] | 161.1 | 0.82 | no | 1.00 | 2,373 / 2,756 |
| 66 | SS_Z8_R20_T150_C2 | gas1_ics1 | 10 | -64 | 6,074 [5,601, 6,578] | 161.2 | 0.75 | no | 1.00 | 1,596 / 2,988 |
| 67 | SS_Z10_R20_T150_C5 | gas1_ics1 | 10 | -63 | 6,131 [5,630, 6,604] | 163.2 | 0.79 | no | 0.99 | 1,983 / 2,885 |
| 68 | SS_Z10_R20_T100000_C2 | gas1_ics1 | 10 | -57 | 6,299 [5,818, 6,786] | 174.3 | 0.77 | no | 1.00 | 1,937 / 2,926 |
| 69 | SS_Z10_R30_T150_C2 | gas1_ics1 | 10 | -15 | 6,131 [5,633, 6,612] | 161.6 | 0.77 | no | 0.99 | 1,775 / 2,949 |
| 70 | SS_Z10_R20_T150_C2 | gas1_ics1 | 10 | -7 | 6,299 [5,656, 6,668] | 165.4 | 0.77 | no | 0.99 | 1,894 / 2,975 |
| 71 | SL_Z6_R30_T150_C2 | gas1_ics1 | 10 | +0 | 6,412 [5,904, 6,898] | 167.9 | 0.82 | no | 0.99 | 2,299 / 2,938 |
| 72 | SL_Z6_R20_T150_C2 | gas1_ics1 | 10 | +11 | 6,412 [5,904, 6,897] | 168.0 | 0.82 | no | 1.00 | 2,314 / 2,938 |
| 73 | SL_Z10_R30_T150_C2 | gas1_ics1 | 10 | +141 | 6,693 [6,179, 7,187] | 178.1 | 0.86 | no | 1.00 | 2,749 / 3,021 |
| 74 | gll_iem_v07 | - | 7 | +787 | 3,150 [2,614, 3,663] | 36.4 | 0.50 | yes | 0.96 | 110 / 1,473 |
