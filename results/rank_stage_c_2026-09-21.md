# Stage C — 2026-09-21

Consensus selection (median Δ(−2lnL) across 257 model rows; BIC threshold p·ln 1080 = 63 per template): **bubble_totani → loopI_shells → erosita_unif_NS**

| step | chosen | median Δ | min | max | threshold | runner-up (median Δ) |
|---|---|---|---|---|---|---|
| 1 | bubble_totani | 28,819 | 20,615 | 39,682 | 126 | bubble_flat (19,194) |
| 2 | loopI_shells | 7,021 | 4,891 | 8,492 | 126 | loopI_shells+nps (7,019) |
| 3 | erosita_unif_NS | 127 | 70 | 275 | 126 | erosita_bright_NS (112) |
| 4 | (stop) loopII | 62 | 0 | 82 | 63 | sun_moon_ic (51) |

## Each model's own path

| model | variant | path |
|---|---|---|
| SL_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SL_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SL_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SL_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SO_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SO_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic → loopII |
| SO_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic → loopII |
| SO_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic → loopII |
| SO_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SO_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SO_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic → loopII |
| SO_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SS_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → erosita_unif_NS → sun_moon_ic |
| SS_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → sun_moon_ic |
| SS_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SS_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SS_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SS_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SS_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SY_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SY_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS |
| SY_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| a2020 | - | bubble_totani → loopI_shells+nps |

## Ranking with the consensus set (bubble_totani, loopI_shells, erosita_unif_NS)

| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | a2020 | - | 22 | -1,408 | 8,493 [7,941, 9,014] | 32.4 | 0.80 | no | 0.99 | 3,915 / 2,842 |
| 2 | SS_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,244 | 5,962 [5,056, 6,972] | 12.1 | 0.64 | yes | 1.00 | 1,608 / 2,197 |
| 3 | SS_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,240 | 6,018 [5,094, 7,004] | 12.2 | 0.67 | yes | 0.99 | 1,667 / 2,367 |
| 4 | SY_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,240 | 6,018 [5,033, 7,071] | 12.3 | 0.71 | yes | 0.99 | 1,772 / 2,473 |
| 5 | SS_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,239 | 6,356 [5,413, 7,297] | 15.7 | 0.87 | no | 0.99 | 2,716 / 2,793 |
| 6 | SS_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,233 | 6,074 [5,119, 7,013] | 13.9 | 0.66 | no | 1.00 | 1,803 / 2,206 |
| 7 | SY_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,232 | 6,018 [5,040, 7,069] | 12.4 | 0.69 | yes | 1.00 | 1,757 / 2,419 |
| 8 | SS_Z8_R20_T150_C2 | gas5_ics3 | 16 | -1,230 | 6,131 [5,165, 7,054] | 14.1 | 0.67 | no | 0.99 | 1,847 / 2,246 |
| 9 | SS_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,227 | 6,356 [5,437, 7,317] | 15.9 | 0.91 | no | 0.99 | 2,790 / 3,022 |
| 10 | SY_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,224 | 6,131 [5,114, 7,114] | 14.1 | 0.73 | no | 0.92 | 1,937 / 2,526 |
| 11 | SS_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,220 | 6,412 [5,486, 7,348] | 17.5 | 0.93 | no | 0.99 | 3,000 / 2,932 |
| 12 | SY_Z8_R20_T150_C2 | gas5_ics3 | 16 | -1,219 | 6,131 [5,149, 7,139] | 14.0 | 0.70 | no | 0.99 | 2,011 / 2,299 |
| 13 | SS_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,218 | 6,412 [5,499, 7,358] | 17.6 | 0.91 | no | 0.99 | 2,895 / 2,920 |
| 14 | SL_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,213 | 6,074 [5,088, 7,013] | 13.0 | 0.71 | no | 0.98 | 2,026 / 2,303 |
| 15 | SY_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,211 | 6,243 [5,258, 7,267] | 13.1 | 0.81 | no | 1.00 | 2,608 / 2,457 |
| 16 | SS_Z6_R30_T150_C2 | gas5_ics3 | 16 | -1,208 | 6,243 [5,290, 7,153] | 18.2 | 0.72 | no | 1.00 | 2,568 / 1,919 |
| 17 | SS_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,207 | 6,468 [5,534, 7,370] | 20.7 | 0.99 | no | 1.00 | 3,613 / 2,808 |
| 18 | SO_Z10_R30_T150_C2 | gas5_ics3 | 16 | -1,206 | 5,906 [5,000, 6,807] | 16.0 | 0.74 | no | 0.99 | 1,925 / 2,432 |
| 19 | SS_Z6_R20_T150_C2 | gas5_ics3 | 16 | -1,205 | 6,243 [5,297, 7,157] | 18.2 | 0.74 | no | 0.99 | 2,568 / 2,034 |
| 20 | SS_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,205 | 6,468 [5,537, 7,373] | 20.7 | 1.01 | no | 1.00 | 3,719 / 2,809 |
| 21 | SL_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,204 | 6,074 [5,104, 7,026] | 13.1 | 0.73 | no | 0.98 | 2,056 / 2,356 |
| 22 | SY_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,202 | 6,243 [5,276, 7,273] | 13.1 | 0.82 | no | 0.67 | 2,578 / 2,533 |
| 23 | SO_Z10_R20_T150_C2 | gas5_ics3 | 16 | -1,201 | 5,906 [5,014, 6,817] | 16.1 | 0.76 | no | 0.99 | 1,969 / 2,523 |
| 24 | SL_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,199 | 6,131 [5,352, 7,287] | 14.8 | 1.01 | no | 0.98 | 3,556 / 2,610 |
| 25 | SY_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,197 | 6,299 [5,337, 7,303] | 14.4 | 0.87 | no | 0.99 | 3,026 / 2,431 |
| 26 | SL_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,195 | 6,131 [5,189, 7,101] | 14.8 | 0.76 | no | 0.99 | 2,205 / 2,432 |
| 27 | SO_Z8_R30_T150_C2 | gas5_ics3 | 16 | -1,193 | 5,962 [5,099, 6,904] | 18.0 | 0.72 | no | 0.99 | 2,209 / 2,110 |
| 28 | SO_Z8_R20_T150_C2 | gas5_ics3 | 16 | -1,191 | 5,962 [5,105, 6,909] | 18.1 | 0.78 | no | 0.99 | 2,224 / 2,404 |
| 29 | SY_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,191 | 6,299 [5,346, 7,307] | 14.4 | 0.87 | no | 1.00 | 2,981 / 2,482 |
| 30 | SO_Z10_R30_T150_C5 | gas5_ics3 | 16 | -1,190 | 6,074 [5,206, 6,991] | 19.9 | 0.93 | no | 0.91 | 2,762 / 2,908 |
| 31 | SL_Z8_R20_T150_C2 | gas5_ics3 | 16 | -1,190 | 6,187 [5,223, 7,132] | 14.8 | 0.75 | no | 0.99 | 2,250 / 2,420 |
| 32 | SL_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,188 | 6,299 [5,389, 7,285] | 15.2 | 0.81 | no | 0.98 | 2,802 / 2,308 |
| 33 | SO_Z10_R20_T150_C5 | gas5_ics3 | 16 | -1,186 | 6,131 [5,216, 6,997] | 20.0 | 0.92 | no | 0.99 | 2,777 / 2,880 |
| 34 | SO_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,179 | 6,187 [5,267, 7,042] | 21.7 | 0.97 | no | 0.99 | 3,137 / 2,895 |
| 35 | SO_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,177 | 6,187 [5,265, 7,038] | 21.7 | 0.98 | no | 0.99 | 3,137 / 2,907 |
| 36 | SL_Z8_R30_T150_C5 | gas5_ics3 | 16 | -1,175 | 6,412 [5,461, 7,346] | 16.6 | 0.83 | no | 0.67 | 2,982 / 2,321 |
| 37 | SO_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,174 | 6,187 [5,319, 7,079] | 24.2 | 1.06 | no | 1.00 | 3,840 / 2,744 |
| 38 | SY_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -1,172 | 6,862 [5,856, 7,895] | 12.3 | 0.60 | yes | 0.99 | 2,037 / 2,089 |
| 39 | SO_Z6_R20_T150_C2 | gas5_ics3 | 16 | -1,172 | 6,074 [5,263, 7,080] | 24.7 | 0.82 | no | 1.00 | 2,914 / 2,055 |
| 40 | SL_Z8_R20_T150_C5 | gas5_ics3 | 16 | -1,171 | 6,412 [5,471, 7,355] | 16.7 | 0.83 | no | 0.99 | 3,012 / 2,333 |
| 41 | SY_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,170 | 6,356 [5,424, 7,332] | 17.8 | 0.92 | no | 0.99 | 3,477 / 2,397 |
| 42 | SY_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,169 | 6,356 [5,428, 7,266] | 17.8 | 0.93 | no | 0.99 | 3,492 / 2,447 |
| 43 | SO_Z6_R30_T150_C2 | gas5_ics3 | 16 | -1,168 | 6,187 [5,246, 7,054] | 22.2 | 0.78 | no | 0.99 | 2,899 / 1,948 |
| 44 | SY_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -1,164 | 6,805 [5,790, 7,791] | 13.4 | 0.62 | yes | 0.99 | 2,126 / 2,123 |
| 45 | SL_Z6_R30_T150_C2 | gas5_ics3 | 16 | -1,162 | 6,299 [5,351, 7,238] | 18.8 | 0.82 | no | 0.88 | 2,775 / 2,367 |
| 46 | SS_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -1,162 | 6,749 [5,799, 7,628] | 15.4 | 0.74 | yes | 0.93 | 2,458 / 2,509 |
| 47 | SS_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -1,161 | 6,974 [6,019, 7,864] | 18.5 | 0.82 | no | 1.00 | 2,924 / 2,808 |
| 48 | SO_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,161 | 6,243 [5,321, 7,066] | 24.4 | 1.05 | no | 0.99 | 3,825 / 2,704 |
| 49 | SS_Z4_R30_T150_C2 | gas5_ics3 | 16 | -1,161 | 6,580 [5,714, 7,355] | 29.8 | 0.91 | no | 0.99 | 3,451 / 2,552 |
| 50 | SS_Z4_R20_T150_C2 | gas5_ics3 | 16 | -1,161 | 6,580 [5,716, 7,359] | 29.8 | 0.92 | no | 0.99 | 3,225 / 2,805 |
| 51 | SS_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -1,161 | 6,805 [5,839, 7,712] | 14.1 | 0.68 | yes | 0.89 | 2,294 / 2,311 |
| 52 | SL_Z6_R20_T150_C2 | gas5_ics3 | 16 | -1,160 | 6,299 [5,357, 7,243] | 18.8 | 0.81 | no | 1.00 | 2,775 / 2,329 |
| 53 | SY_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -1,158 | 6,918 [5,901, 7,930] | 12.4 | 0.63 | yes | 0.99 | 2,185 / 2,141 |
| 54 | SY_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -1,157 | 6,805 [5,808, 7,802] | 13.5 | 0.64 | yes | 0.99 | 2,349 / 2,023 |
| 55 | SS_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -1,150 | 6,862 [5,900, 7,761] | 14.4 | 0.66 | yes | 0.99 | 2,427 / 2,071 |
| 56 | SY_Z6_R30_T150_C2 | gas5_ics3 | 16 | -1,148 | 6,356 [5,416, 7,290] | 19.0 | 0.81 | no | 0.99 | 2,655 / 2,522 |
| 57 | SY_Z6_R20_T150_C2 | gas5_ics3 | 16 | -1,147 | 6,356 [5,420, 7,293] | 19.0 | 0.84 | no | 1.00 | 2,670 / 2,699 |
| 58 | SO_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -1,144 | 6,637 [5,732, 7,507] | 21.0 | 0.85 | no | 0.97 | 2,877 / 2,795 |
| 59 | SS_Z4_R20_T150_C5 | gas5_ics3 | 16 | -1,143 | 6,637 [5,786, 7,392] | 30.2 | 1.11 | no | 0.99 | 4,365 / 2,987 |
| 60 | SS_Z4_R30_T150_C5 | gas5_ics3 | 16 | -1,143 | 6,637 [5,785, 7,389] | 30.2 | 1.06 | no | 0.95 | 3,943 / 3,076 |
| 61 | SL_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -1,140 | 6,918 [5,898, 7,953] | 12.7 | 0.61 | yes | 1.00 | 2,273 / 1,914 |
| 62 | SS_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -1,136 | 6,805 [5,901, 7,738] | 18.5 | 0.81 | no | 1.00 | 3,086 / 2,434 |
| 63 | SO_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -1,134 | 6,693 [5,783, 7,566] | 20.5 | 0.81 | no | 0.88 | 2,611 / 2,782 |
| 64 | SL_Z6_R20_T150_C5 | gas5_ics3 | 16 | -1,133 | 6,524 [5,572, 7,412] | 20.6 | 0.87 | no | 0.99 | 3,417 / 2,245 |
| 65 | SS_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -1,132 | 6,805 [5,912, 7,747] | 18.5 | 0.82 | no | 0.98 | 3,100 / 2,447 |
| 66 | SL_Z6_R30_T150_C5 | gas5_ics3 | 16 | -1,131 | 6,468 [5,564, 7,406] | 20.6 | 0.89 | no | 0.66 | 3,418 / 2,323 |
| 67 | SO_Z4_R20_T150_C5 | gas5_ics3 | 16 | -1,125 | 6,412 [5,861, 6,997] | 32.0 | 1.05 | no | 0.96 | 4,171 / 2,589 |
| 68 | SL_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -1,124 | 6,862 [5,873, 7,877] | 14.2 | 0.62 | yes | 0.66 | 2,199 / 2,039 |
| 69 | SY_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -1,123 | 6,805 [5,861, 7,792] | 15.9 | 0.74 | yes | 0.99 | 2,813 / 2,224 |
| 70 | SY_Z4_R20_T150_C5 | gas5_ics3 | 16 | -1,122 | 6,580 [5,718, 7,354] | 26.1 | 0.98 | no | 1.00 | 3,957 / 2,469 |
| 71 | SL_Z10_R20_T100000_C2 | gas5_ics3 | 16 | -1,121 | 6,918 [5,931, 7,972] | 12.7 | 0.60 | yes | 0.67 | 2,199 / 1,952 |
| 72 | SO_Z4_R30_T150_C5 | gas5_ics3 | 16 | -1,121 | 6,412 [5,859, 6,995] | 32.1 | 1.04 | no | 0.97 | 4,035 / 2,601 |
| 73 | SY_Z10_R30_T100000_C5 | gas5_ics3 | 16 | -1,121 | 7,199 [6,151, 8,125] | 17.2 | 0.70 | no | 0.99 | 2,711 / 2,319 |
| 74 | SL_Z8_R20_T100000_C2 | gas5_ics3 | 16 | -1,120 | 6,862 [5,886, 7,881] | 14.2 | 0.63 | yes | 0.86 | 2,244 / 2,052 |
| 75 | SY_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -1,119 | 6,805 [5,868, 7,797] | 15.9 | 0.74 | yes | 1.00 | 2,813 / 2,225 |
| 76 | SY_Z4_R30_T150_C5 | gas5_ics3 | 16 | -1,117 | 6,580 [5,717, 7,352] | 26.1 | 0.96 | no | 1.00 | 3,853 / 2,482 |
| 77 | SS_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -1,117 | 6,974 [6,042, 7,861] | 22.5 | 0.84 | no | 0.99 | 3,118 / 2,714 |
| 78 | SO_Z8_R30_T100000_C2 | gas5_ics3 | 16 | -1,116 | 6,524 [5,655, 7,437] | 19.7 | 0.79 | no | 1.00 | 2,476 / 2,705 |
| 79 | SS_Z10_R30_T100000_C5 | gas5_ics3 | 16 | -1,111 | 7,030 [6,097, 7,960] | 21.5 | 0.82 | no | 0.99 | 3,103 / 2,639 |
| 80 | SO_Z10_R30_T100000_C2 | gas5_ics3 | 16 | -1,110 | 6,580 [5,688, 7,478] | 18.5 | 0.76 | no | 0.90 | 2,371 / 2,617 |
| 81 | SY_Z4_R30_T150_C2 | gas5_ics3 | 16 | -1,106 | 6,637 [5,753, 7,403] | 28.5 | 0.94 | no | 1.00 | 3,661 / 2,604 |
| 82 | SO_Z4_R20_T150_C2 | gas5_ics3 | 16 | -1,104 | 6,468 [5,928, 7,063] | 34.5 | 1.00 | no | 0.99 | 3,871 / 2,613 |
| 83 | SO_Z4_R30_T150_C2 | gas5_ics3 | 16 | -1,104 | 6,524 [5,930, 7,066] | 34.5 | 0.99 | no | 0.99 | 3,856 / 2,613 |
| 84 | SS_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -1,102 | 7,087 [6,136, 7,991] | 19.1 | 0.80 | no | 0.99 | 2,969 / 2,705 |
| 85 | SS_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -1,101 | 7,030 [6,092, 7,911] | 20.1 | 0.84 | no | 0.92 | 3,162 / 2,768 |
| 86 | SY_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -1,101 | 7,030 [6,057, 8,007] | 17.8 | 0.90 | no | 0.66 | 3,553 / 2,744 |
| 87 | SY_Z4_R20_T150_C2 | gas5_ics3 | 16 | -1,097 | 6,637 [5,754, 7,405] | 28.5 | 0.97 | no | 0.89 | 3,887 / 2,579 |
| 88 | SL_Z10_R30_T100000_C5 | gas5_ics3 | 16 | -1,095 | 6,974 [5,995, 8,025] | 15.1 | 0.79 | no | 0.67 | 3,327 / 2,148 |
| 89 | SL_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -1,088 | 6,918 [5,939, 7,915] | 16.4 | 0.66 | yes | 0.99 | 2,544 / 1,999 |
| 90 | SS_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -1,087 | 6,974 [6,085, 7,900] | 22.3 | 0.91 | no | 1.00 | 3,462 / 2,855 |
| 91 | SY_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -1,087 | 6,918 [5,962, 7,923] | 15.3 | 0.67 | yes | 0.98 | 2,456 / 2,199 |
| 92 | SY_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -1,087 | 6,974 [5,972, 7,961] | 15.8 | 0.65 | yes | 0.67 | 2,411 / 2,151 |
| 93 | SS_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -1,084 | 6,862 [6,006, 7,682] | 27.7 | 0.96 | no | 1.00 | 4,032 / 2,587 |
| 94 | SO_Z8_R30_T100000_C5 | gas5_ics3 | 16 | -1,084 | 6,580 [5,714, 7,479] | 21.2 | 0.84 | no | 0.99 | 3,045 / 2,500 |
| 95 | SS_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -1,083 | 6,862 [6,010, 7,687] | 27.6 | 0.97 | no | 0.99 | 4,032 / 2,623 |
| 96 | SO_Z10_R30_T100000_C5 | gas5_ics3 | 16 | -1,083 | 6,580 [5,722, 7,498] | 20.8 | 0.87 | no | 0.99 | 2,999 / 2,718 |
| 97 | SL_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -1,083 | 6,918 [5,894, 8,026] | 15.2 | 0.63 | yes | 0.98 | 2,365 / 1,965 |
| 98 | SO_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -1,082 | 6,580 [5,687, 7,404] | 22.3 | 0.90 | no | 0.99 | 3,404 / 2,517 |
| 99 | SO_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -1,082 | 6,637 [5,745, 7,500] | 22.9 | 0.84 | no | 1.00 | 3,164 / 2,379 |
| 100 | SO_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -1,082 | 6,580 [5,684, 7,399] | 22.3 | 0.88 | no | 1.00 | 3,389 / 2,427 |
| 101 | SS_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -1,081 | 6,974 [6,077, 7,894] | 22.2 | 0.89 | no | 0.99 | 3,447 / 2,738 |
| 102 | SL_Z4_R30_T150_C2 | gas5_ics3 | 16 | -1,076 | 6,637 [5,751, 7,401] | 28.6 | 0.90 | no | 1.00 | 3,523 / 2,444 |
| 103 | SO_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -1,076 | 6,637 [5,746, 7,502] | 22.9 | 0.85 | no | 1.00 | 3,208 / 2,444 |
| 104 | SO_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -1,075 | 6,637 [5,744, 7,504] | 21.3 | 0.87 | no | 1.00 | 3,089 / 2,692 |
| 105 | SL_Z4_R20_T150_C2 | gas5_ics3 | 16 | -1,075 | 6,637 [5,756, 7,406] | 28.6 | 0.90 | no | 1.00 | 3,510 / 2,456 |
| 106 | SO_Z10_R20_T100000_C5 | gas5_ics3 | 16 | -1,075 | 6,637 [5,739, 7,510] | 21.0 | 0.83 | no | 0.89 | 2,820 / 2,656 |
| 107 | SL_Z6_R30_T100000_C2 | gas5_ics3 | 16 | -1,073 | 6,974 [6,036, 7,962] | 18.3 | 0.86 | no | 0.99 | 3,791 / 2,207 |
| 108 | SL_Z6_R20_T100000_C2 | gas5_ics3 | 16 | -1,071 | 6,918 [6,010, 7,862] | 17.2 | 0.69 | yes | 0.99 | 2,722 / 2,018 |
| 109 | SO_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -1,068 | 6,580 [5,830, 7,373] | 29.9 | 1.01 | no | 0.99 | 4,124 / 2,536 |
| 110 | SL_Z8_R20_T100000_C5 | gas5_ics3 | 16 | -1,066 | 7,030 [6,030, 8,001] | 16.5 | 0.62 | yes | 0.99 | 2,335 / 2,014 |
| 111 | SY_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -1,065 | 7,030 [6,103, 8,002] | 20.4 | 0.78 | no | 0.66 | 3,115 / 2,336 |
| 112 | SY_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -1,061 | 6,918 [5,978, 7,885] | 17.8 | 0.74 | no | 0.98 | 2,965 / 2,139 |
| 113 | SL_Z4_R20_T150_C5 | gas5_ics3 | 16 | -1,048 | 7,199 [6,180, 7,562] | 27.5 | 0.88 | no | 1.00 | 3,868 / 2,486 |
| 114 | SL_Z4_R30_T150_C5 | gas5_ics3 | 16 | -1,043 | 6,693 [5,828, 7,414] | 27.4 | 0.92 | no | 0.99 | 3,778 / 2,396 |
| 115 | SY_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -1,036 | 6,862 [5,987, 7,694] | 24.5 | 0.86 | no | 1.00 | 3,611 / 2,258 |
| 116 | SO_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -1,034 | 6,749 [5,871, 7,340] | 30.4 | 0.96 | no | 1.00 | 4,064 / 2,407 |
| 117 | SY_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -1,033 | 6,862 [5,989, 7,809] | 24.5 | 0.88 | no | 1.00 | 3,731 / 2,285 |
| 118 | SS_Z4_R30_T100000_C5 | gas5_ics3 | 16 | -1,032 | 6,862 [6,017, 7,669] | 28.0 | 0.99 | no | 0.99 | 4,138 / 2,639 |
| 119 | SL_Z6_R30_T100000_C5 | gas5_ics3 | 16 | -1,027 | 6,974 [6,054, 7,945] | 19.1 | 0.86 | no | 0.66 | 3,672 / 2,324 |
| 120 | SL_Z6_R20_T100000_C5 | gas5_ics3 | 16 | -1,026 | 7,143 [6,034, 7,956] | 18.3 | 0.69 | yes | 0.98 | 2,874 / 2,020 |
| 121 | SO_Z4_R20_T100000_C5 | gas5_ics3 | 16 | -1,014 | 6,749 [6,178, 7,293] | 28.4 | 0.97 | no | 1.00 | 4,066 / 2,500 |
| 122 | SS_Z4_R20_T100000_C5 | gas5_ics3 | 16 | -1,014 | 6,918 [6,081, 7,722] | 27.7 | 1.01 | no | 0.99 | 4,153 / 2,856 |
| 123 | SO_Z4_R30_T100000_C5 | gas5_ics3 | 16 | -1,013 | 6,749 [6,177, 7,292] | 28.4 | 0.95 | no | 0.99 | 3,930 / 2,486 |
| 124 | SL_Z4_R30_T100000_C2 | gas5_ics3 | 16 | -988 | 7,030 [6,153, 7,871] | 25.8 | 0.81 | no | 0.79 | 3,505 / 2,200 |
| 125 | SL_Z4_R20_T100000_C2 | gas5_ics3 | 16 | -987 | 7,030 [6,155, 7,872] | 25.8 | 0.80 | no | 1.00 | 3,400 / 2,226 |
| 126 | SY_Z4_R20_T100000_C5 | gas5_ics3 | 16 | -987 | 7,030 [6,140, 7,799] | 25.0 | 0.88 | no | 1.00 | 3,882 / 2,289 |
| 127 | SY_Z4_R30_T100000_C5 | gas5_ics3 | 16 | -987 | 7,030 [6,137, 7,796] | 25.0 | 0.88 | no | 0.96 | 3,882 / 2,289 |
| 128 | SL_Z4_R20_T100000_C5 | gas5_ics3 | 16 | -932 | 7,087 [6,215, 7,852] | 25.9 | 0.98 | no | 1.00 | 4,454 / 2,468 |
| 129 | SL_Z4_R30_T100000_C5 | gas5_ics3 | 16 | -930 | 7,087 [6,201, 7,841] | 24.8 | 0.83 | no | 0.66 | 3,731 / 2,152 |
| 130 | SY_Z4_R30_T100000_C5 | gas1_ics1 | 10 | -343 | 5,906 [5,370, 6,424] | 130.2 | 0.93 | no | 1.00 | 2,707 / 2,766 |
| 131 | SY_Z4_R20_T100000_C5 | gas1_ics1 | 10 | -338 | 5,906 [5,366, 6,420] | 130.0 | 0.93 | no | 1.00 | 2,721 / 2,766 |
| 132 | SS_Z4_R30_T100000_C5 | gas1_ics1 | 10 | -321 | 5,793 [5,277, 6,335] | 125.7 | 0.95 | no | 1.00 | 2,603 / 2,887 |
| 133 | SS_Z4_R20_T100000_C5 | gas1_ics1 | 10 | -318 | 5,849 [5,303, 6,360] | 127.4 | 0.92 | no | 0.99 | 2,484 / 2,915 |
| 134 | SY_Z4_R30_T150_C5 | gas1_ics1 | 10 | -306 | 5,737 [5,227, 6,287] | 122.3 | 0.95 | no | 1.00 | 2,527 / 2,938 |
| 135 | SY_Z4_R20_T150_C5 | gas1_ics1 | 10 | -303 | 5,737 [5,223, 6,284] | 122.2 | 0.95 | no | 1.00 | 2,527 / 2,938 |
| 136 | SL_Z4_R30_T100000_C5 | gas1_ics1 | 10 | -296 | 5,906 [5,370, 6,424] | 129.6 | 0.94 | no | 1.00 | 2,960 / 2,587 |
| 137 | SL_Z4_R20_T100000_C5 | gas1_ics1 | 10 | -295 | 5,906 [5,364, 6,418] | 129.3 | 0.94 | no | 1.00 | 2,945 / 2,587 |
| 138 | SY_Z4_R30_T100000_C2 | gas1_ics1 | 10 | -278 | 6,018 [5,471, 6,525] | 134.8 | 0.92 | no | 1.00 | 2,647 / 2,881 |
| 139 | SS_Z6_R30_T100000_C5 | gas1_ics1 | 10 | -278 | 5,849 [5,325, 6,378] | 128.8 | 0.92 | no | 1.00 | 2,483 / 2,916 |
| 140 | SY_Z4_R20_T100000_C2 | gas1_ics1 | 10 | -274 | 6,018 [5,468, 6,522] | 134.7 | 0.92 | no | 1.00 | 2,647 / 2,919 |
| 141 | SS_Z4_R30_T150_C5 | gas1_ics1 | 10 | -269 | 5,737 [5,217, 6,279] | 122.9 | 0.95 | no | 1.00 | 2,304 / 3,126 |
| 142 | SS_Z4_R20_T150_C5 | gas1_ics1 | 10 | -264 | 5,737 [5,215, 6,276] | 122.9 | 1.00 | no | 1.00 | 2,602 / 3,112 |
| 143 | SY_Z6_R30_T100000_C5 | gas1_ics1 | 10 | -257 | 6,018 [5,494, 6,551] | 134.4 | 0.95 | no | 1.00 | 2,810 / 2,886 |
| 144 | SS_Z6_R20_T100000_C5 | gas1_ics1 | 10 | -253 | 5,849 [5,329, 6,380] | 129.3 | 0.95 | no | 1.00 | 2,617 / 2,916 |
| 145 | SS_Z8_R30_T100000_C5 | gas1_ics1 | 10 | -244 | 5,849 [5,338, 6,392] | 129.2 | 0.94 | no | 1.00 | 2,632 / 2,850 |
| 146 | SL_Z4_R30_T150_C5 | gas1_ics1 | 10 | -241 | 5,793 [5,252, 6,311] | 123.1 | 0.97 | no | 1.00 | 2,856 / 2,745 |
| 147 | SL_Z4_R20_T150_C5 | gas1_ics1 | 10 | -241 | 5,793 [5,247, 6,307] | 122.9 | 0.97 | no | 1.00 | 2,856 / 2,745 |
| 148 | SY_Z6_R20_T100000_C5 | gas1_ics1 | 10 | -240 | 6,018 [5,493, 6,550] | 134.7 | 0.94 | no | 1.00 | 2,840 / 2,818 |
| 149 | SL_Z4_R30_T100000_C2 | gas1_ics1 | 10 | -240 | 6,074 [5,529, 6,582] | 137.2 | 0.94 | no | 1.00 | 2,945 / 2,753 |
| 150 | SL_Z4_R20_T100000_C2 | gas1_ics1 | 10 | -239 | 6,074 [5,522, 6,576] | 136.8 | 0.93 | no | 1.00 | 2,931 / 2,740 |
| 151 | SS_Z4_R30_T100000_C2 | gas1_ics1 | 10 | -228 | 5,849 [5,341, 6,396] | 129.4 | 0.92 | no | 1.00 | 2,349 / 3,055 |
| 152 | SL_Z6_R30_T100000_C5 | gas1_ics1 | 10 | -227 | 6,074 [5,550, 6,614] | 134.9 | 0.96 | no | 0.99 | 3,124 / 2,679 |
| 153 | SS_Z4_R20_T100000_C2 | gas1_ics1 | 10 | -221 | 5,849 [5,340, 6,395] | 129.5 | 0.92 | no | 1.00 | 2,453 / 2,949 |
| 154 | SS_Z10_R30_T100000_C5 | gas1_ics1 | 10 | -219 | 5,906 [5,397, 6,448] | 132.2 | 0.94 | no | 0.99 | 2,737 / 2,824 |
| 155 | SS_Z8_R20_T100000_C5 | gas1_ics1 | 10 | -216 | 5,906 [5,364, 6,414] | 131.2 | 0.94 | no | 0.99 | 2,677 / 2,876 |
| 156 | SL_Z6_R20_T100000_C5 | gas1_ics1 | 10 | -215 | 6,074 [5,546, 6,610] | 134.9 | 0.96 | no | 0.99 | 3,124 / 2,679 |
| 157 | SS_Z6_R30_T150_C5 | gas1_ics1 | 10 | -212 | 5,737 [5,205, 6,262] | 122.6 | 0.98 | no | 1.00 | 2,513 / 3,086 |
| 158 | SY_Z4_R30_T150_C2 | gas1_ics1 | 10 | -206 | 5,849 [5,313, 6,373] | 126.3 | 0.96 | no | 1.00 | 2,498 / 3,105 |
| 159 | SS_Z6_R20_T150_C5 | gas1_ics1 | 10 | -205 | 5,737 [5,209, 6,265] | 123.0 | 0.98 | no | 1.00 | 2,527 / 3,086 |
| 160 | SY_Z4_R20_T150_C2 | gas1_ics1 | 10 | -203 | 5,849 [5,310, 6,370] | 126.2 | 0.96 | no | 1.00 | 2,498 / 3,092 |
| 161 | SY_Z8_R30_T100000_C5 | gas1_ics1 | 10 | -199 | 6,074 [5,557, 6,620] | 136.0 | 0.95 | no | 1.00 | 2,930 / 2,819 |
| 162 | SO_Z10_R30_T100000_C5 | gas1_ics1 | 10 | -196 | 5,624 [5,091, 6,140] | 119.3 | 0.94 | no | 0.98 | 2,647 / 2,648 |
| 163 | SO_Z8_R30_T100000_C5 | gas1_ics1 | 10 | -192 | 5,624 [5,087, 6,136] | 119.4 | 0.95 | no | 1.00 | 2,602 / 2,751 |
| 164 | SS_Z6_R30_T100000_C2 | gas1_ics1 | 10 | -188 | 5,906 [5,398, 6,448] | 133.0 | 0.95 | no | 1.00 | 2,647 / 2,940 |
| 165 | SL_Z8_R30_T100000_C5 | gas1_ics1 | 10 | -181 | 6,187 [5,647, 6,719] | 137.0 | 0.96 | no | 0.99 | 3,244 / 2,719 |
| 166 | SS_Z8_R30_T150_C5 | gas1_ics1 | 10 | -179 | 5,793 [5,275, 6,329] | 126.1 | 0.99 | no | 0.99 | 2,692 / 3,033 |
| 167 | SO_Z4_R20_T100000_C5 | gas1_ics1 | 10 | -175 | 5,512 [4,984, 6,036] | 114.4 | 0.91 | no | 0.99 | 2,350 / 2,689 |
| 168 | SO_Z4_R30_T100000_C5 | gas1_ics1 | 10 | -174 | 5,512 [4,987, 6,040] | 114.6 | 0.96 | no | 0.99 | 2,603 / 2,689 |
| 169 | SO_Z6_R30_T100000_C5 | gas1_ics1 | 10 | -174 | 5,568 [5,054, 6,105] | 117.8 | 0.93 | no | 0.99 | 2,543 / 2,638 |
| 170 | SO_Z8_R20_T100000_C5 | gas1_ics1 | 10 | -167 | 5,624 [5,109, 6,156] | 120.8 | 0.96 | no | 0.99 | 2,632 / 2,764 |
| 171 | SY_Z6_R30_T150_C5 | gas1_ics1 | 10 | -167 | 5,849 [5,343, 6,404] | 127.4 | 0.99 | no | 0.99 | 2,766 / 3,042 |
| 172 | SS_Z10_R20_T100000_C5 | gas1_ics1 | 10 | -167 | 5,962 [5,419, 6,467] | 134.2 | 0.94 | no | 0.98 | 2,766 / 2,837 |
| 173 | SS_Z6_R20_T100000_C2 | gas1_ics1 | 10 | -166 | 5,906 [5,404, 6,453] | 133.6 | 0.95 | no | 1.00 | 2,647 / 2,967 |
| 174 | SS_Z4_R30_T150_C2 | gas1_ics1 | 10 | -165 | 5,793 [5,237, 6,298] | 124.3 | 0.96 | no | 0.98 | 2,349 / 3,189 |
| 175 | SY_Z6_R30_T100000_C2 | gas1_ics1 | 10 | -161 | 6,131 [5,628, 6,683] | 142.1 | 0.96 | no | 1.00 | 2,886 / 2,973 |
| 176 | SS_Z4_R20_T150_C2 | gas1_ics1 | 10 | -160 | 5,793 [5,235, 6,296] | 124.3 | 0.96 | no | 0.98 | 2,349 / 3,217 |
| 177 | SO_Z10_R20_T100000_C5 | gas1_ics1 | 10 | -155 | 5,624 [5,117, 6,164] | 121.0 | 0.96 | no | 0.99 | 2,662 / 2,754 |
| 178 | SS_Z8_R30_T100000_C2 | gas1_ics1 | 10 | -154 | 5,962 [5,419, 6,469] | 133.6 | 0.93 | no | 1.00 | 2,722 / 2,840 |
| 179 | SY_Z8_R20_T100000_C5 | gas1_ics1 | 10 | -152 | 6,131 [5,588, 6,647] | 138.4 | 0.95 | no | 1.00 | 2,975 / 2,832 |
| 180 | SY_Z6_R20_T150_C5 | gas1_ics1 | 10 | -152 | 5,849 [5,343, 6,403] | 127.5 | 0.99 | no | 0.99 | 2,766 / 3,015 |
| 181 | SY_Z6_R20_T100000_C2 | gas1_ics1 | 10 | -151 | 6,131 [5,630, 6,685] | 142.5 | 0.96 | no | 1.00 | 2,916 / 2,973 |
| 182 | SL_Z4_R30_T150_C2 | gas1_ics1 | 10 | -149 | 5,906 [5,394, 6,452] | 129.6 | 0.97 | no | 1.00 | 2,796 / 2,913 |
| 183 | SL_Z4_R20_T150_C2 | gas1_ics1 | 10 | -148 | 5,906 [5,390, 6,448] | 129.4 | 0.97 | no | 1.00 | 2,796 / 2,913 |
| 184 | SS_Z10_R30_T150_C5 | gas1_ics1 | 10 | -148 | 5,849 [5,302, 6,356] | 127.2 | 0.99 | no | 0.99 | 2,782 / 2,994 |
| 185 | SO_Z6_R20_T100000_C5 | gas1_ics1 | 10 | -148 | 5,568 [5,057, 6,107] | 118.0 | 0.91 | no | 0.99 | 2,543 / 2,549 |
| 186 | SL_Z6_R30_T100000_C2 | gas1_ics1 | 10 | -147 | 6,243 [5,740, 6,802] | 144.5 | 0.97 | no | 1.00 | 3,230 / 2,821 |
| 187 | SS_Z8_R20_T150_C5 | gas1_ics1 | 10 | -146 | 5,793 [5,290, 6,342] | 127.3 | 0.99 | no | 0.97 | 2,722 / 3,033 |
| 188 | SY_Z10_R30_T100000_C5 | gas1_ics1 | 10 | -142 | 6,187 [5,653, 6,715] | 140.1 | 0.96 | no | 1.00 | 3,065 / 2,862 |
| 189 | SL_Z8_R20_T100000_C5 | gas1_ics1 | 10 | -138 | 6,187 [5,676, 6,747] | 139.2 | 0.98 | no | 0.99 | 3,304 / 2,733 |
| 190 | SL_Z6_R20_T100000_C2 | gas1_ics1 | 10 | -135 | 6,243 [5,740, 6,801] | 144.7 | 0.97 | no | 1.00 | 3,245 / 2,821 |
| 191 | SS_Z10_R30_T100000_C2 | gas1_ics1 | 10 | -135 | 6,018 [5,483, 6,532] | 136.8 | 0.96 | no | 0.98 | 2,842 / 2,928 |
| 192 | SL_Z10_R30_T100000_C5 | gas1_ics1 | 10 | -122 | 6,299 [5,759, 6,834] | 141.3 | 0.97 | no | 0.99 | 3,379 / 2,735 |
| 193 | SY_Z8_R30_T100000_C2 | gas1_ics1 | 10 | -119 | 6,243 [5,688, 6,748] | 142.8 | 0.97 | no | 1.00 | 3,080 / 2,975 |
| 194 | SS_Z8_R20_T100000_C2 | gas1_ics1 | 10 | -115 | 5,962 [5,437, 6,484] | 135.1 | 0.96 | no | 0.97 | 2,752 / 2,952 |
| 195 | SS_Z6_R30_T150_C2 | gas1_ics1 | 10 | -108 | 5,737 [5,231, 6,287] | 124.1 | 0.99 | no | 0.98 | 2,528 / 3,163 |
| 196 | SL_Z8_R30_T100000_C2 | gas1_ics1 | 10 | -101 | 6,356 [5,831, 6,901] | 146.1 | 0.97 | no | 1.00 | 3,349 / 2,835 |
| 197 | SS_Z10_R20_T150_C5 | gas1_ics1 | 10 | -100 | 5,849 [5,334, 6,384] | 129.3 | 1.00 | no | 0.99 | 2,827 / 3,007 |
| 198 | SO_Z10_R30_T100000_C2 | gas1_ics1 | 10 | -100 | 5,681 [5,169, 6,216] | 123.4 | 0.95 | no | 1.00 | 2,737 / 2,649 |
| 199 | SY_Z10_R20_T100000_C5 | gas1_ics1 | 10 | -91 | 6,187 [5,664, 6,724] | 141.5 | 0.97 | no | 0.99 | 3,095 / 2,887 |
| 200 | SL_Z6_R30_T150_C5 | gas1_ics1 | 10 | -89 | 5,962 [5,438, 6,504] | 129.4 | 1.00 | no | 0.99 | 3,110 / 2,824 |
| 201 | SS_Z6_R20_T150_C2 | gas1_ics1 | 10 | -88 | 5,737 [5,237, 6,292] | 124.6 | 0.99 | no | 0.98 | 2,543 / 3,162 |
| 202 | SO_Z8_R30_T100000_C2 | gas1_ics1 | 10 | -85 | 5,681 [5,130, 6,179] | 121.5 | 0.98 | no | 0.99 | 2,662 / 2,880 |
| 203 | SY_Z8_R20_T100000_C2 | gas1_ics1 | 10 | -84 | 6,243 [5,697, 6,755] | 143.9 | 0.97 | no | 1.00 | 3,080 / 2,960 |
| 204 | SL_Z6_R20_T150_C5 | gas1_ics1 | 10 | -80 | 5,962 [5,439, 6,504] | 129.6 | 1.00 | no | 1.00 | 3,125 / 2,836 |
| 205 | SL_Z10_R20_T100000_C5 | gas1_ics1 | 10 | -78 | 6,299 [5,764, 6,837] | 142.2 | 0.98 | no | 0.99 | 3,409 / 2,747 |
| 206 | SL_Z8_R20_T100000_C2 | gas1_ics1 | 10 | -74 | 6,356 [5,835, 6,904] | 146.7 | 0.98 | no | 0.99 | 3,379 / 2,848 |
| 207 | SY_Z8_R30_T150_C5 | gas1_ics1 | 10 | -73 | 6,018 [5,484, 6,545] | 133.7 | 1.01 | no | 0.99 | 3,021 / 3,056 |
| 208 | SO_Z10_R30_T150_C5 | gas1_ics1 | 10 | -68 | 5,624 [5,102, 6,155] | 118.9 | 1.00 | no | 0.98 | 2,722 / 2,884 |
| 209 | SO_Z8_R30_T150_C5 | gas1_ics1 | 10 | -68 | 5,624 [5,098, 6,151] | 118.9 | 0.95 | no | 1.00 | 2,692 / 2,644 |
| 210 | SY_Z10_R30_T100000_C2 | gas1_ics1 | 10 | -68 | 6,299 [5,788, 6,849] | 147.0 | 0.99 | no | 1.00 | 3,215 / 2,991 |
| 211 | SS_Z10_R20_T100000_C2 | gas1_ics1 | 10 | -66 | 6,018 [5,521, 6,566] | 139.8 | 0.96 | no | 0.99 | 2,886 / 2,902 |
| 212 | SS_Z8_R30_T150_C2 | gas1_ics1 | 10 | -66 | 5,793 [5,283, 6,337] | 126.4 | 1.00 | no | 0.98 | 2,678 / 3,111 |
| 213 | SY_Z6_R30_T150_C2 | gas1_ics1 | 10 | -63 | 5,962 [5,413, 6,464] | 132.7 | 0.99 | no | 1.00 | 2,737 / 3,170 |
| 214 | SO_Z8_R20_T100000_C2 | gas1_ics1 | 10 | -62 | 5,681 [5,160, 6,205] | 123.3 | 0.96 | no | 0.99 | 2,692 / 2,776 |
| 215 | SO_Z10_R20_T100000_C2 | gas1_ics1 | 10 | -56 | 5,737 [5,199, 6,244] | 125.3 | 0.97 | no | 0.98 | 2,767 / 2,802 |
| 216 | SO_Z6_R30_T100000_C2 | gas1_ics1 | 10 | -52 | 5,624 [5,117, 6,165] | 121.1 | 0.97 | no | 0.99 | 2,602 / 2,879 |
| 217 | SY_Z6_R20_T150_C2 | gas1_ics1 | 10 | -47 | 5,962 [5,416, 6,466] | 133.0 | 0.99 | no | 1.00 | 2,752 / 3,157 |
| 218 | SO_Z8_R20_T150_C5 | gas1_ics1 | 10 | -45 | 5,624 [5,110, 6,162] | 119.6 | 0.94 | no | 1.00 | 2,692 / 2,618 |
| 219 | SL_Z10_R30_T100000_C2 | gas1_ics1 | 10 | -44 | 6,468 [5,949, 7,022] | 150.7 | 0.99 | no | 0.99 | 3,544 / 2,864 |
| 220 | SO_Z6_R20_T100000_C2 | gas1_ics1 | 10 | -42 | 5,624 [5,120, 6,168] | 121.4 | 0.92 | no | 0.99 | 2,602 / 2,550 |
| 221 | SY_Z8_R20_T150_C5 | gas1_ics1 | 10 | -41 | 6,018 [5,492, 6,551] | 134.4 | 1.01 | no | 0.99 | 3,036 / 3,056 |
| 222 | SS_Z10_R30_T150_C2 | gas1_ics1 | 10 | -36 | 5,849 [5,335, 6,388] | 128.8 | 1.00 | no | 0.97 | 2,782 / 3,058 |
| 223 | SS_Z8_R20_T150_C2 | gas1_ics1 | 10 | -36 | 5,849 [5,316, 6,368] | 128.5 | 1.00 | no | 1.00 | 2,722 / 3,111 |
| 224 | SO_Z10_R20_T150_C5 | gas1_ics1 | 10 | -35 | 5,624 [5,127, 6,179] | 120.4 | 1.01 | no | 0.99 | 2,767 / 2,897 |
| 225 | SO_Z4_R20_T100000_C2 | gas1_ics1 | 10 | -33 | 5,568 [5,031, 6,083] | 116.6 | 0.92 | no | 0.99 | 2,498 / 2,626 |
| 226 | SO_Z4_R30_T100000_C2 | gas1_ics1 | 10 | -31 | 5,568 [5,035, 6,087] | 116.8 | 0.92 | no | 0.99 | 2,498 / 2,626 |
| 227 | SO_Z6_R30_T150_C5 | gas1_ics1 | 10 | -27 | 5,568 [5,056, 6,111] | 116.9 | 0.94 | no | 0.99 | 2,543 / 2,708 |
| 228 | SO_Z4_R20_T150_C5 | gas1_ics1 | 10 | -24 | 5,512 [5,008, 6,067] | 114.4 | 0.93 | no | 1.00 | 2,468 / 2,684 |
| 229 | SO_Z6_R20_T150_C5 | gas1_ics1 | 10 | -23 | 5,568 [5,061, 6,116] | 117.2 | 0.95 | no | 0.99 | 2,557 / 2,708 |
| 230 | SO_Z4_R30_T150_C5 | gas1_ics1 | 10 | -21 | 5,568 [5,012, 6,071] | 114.6 | 0.93 | no | 1.00 | 2,483 / 2,684 |
| 231 | SY_Z10_R20_T100000_C2 | gas1_ics1 | 10 | -12 | 6,356 [5,816, 6,872] | 149.4 | 0.98 | no | 0.99 | 3,260 / 2,989 |
| 232 | SY_Z10_R30_T150_C5 | gas1_ics1 | 10 | -11 | 6,074 [5,547, 6,611] | 135.6 | 1.02 | no | 0.99 | 3,170 / 3,044 |
| 233 | SL_Z6_R30_T150_C2 | gas1_ics1 | 10 | +0 | 6,074 [5,567, 6,634] | 135.3 | 0.99 | no | 0.99 | 3,051 / 2,992 |
| 234 | SL_Z10_R20_T100000_C2 | gas1_ics1 | 10 | +1 | 6,524 [5,967, 7,038] | 152.4 | 0.99 | no | 0.99 | 3,589 / 2,850 |
| 235 | SL_Z8_R30_T150_C5 | gas1_ics1 | 10 | +4 | 6,131 [5,613, 6,684] | 135.9 | 1.01 | no | 0.99 | 3,379 / 2,840 |
| 236 | SS_Z10_R20_T150_C2 | gas1_ics1 | 10 | +12 | 5,906 [5,376, 6,425] | 131.4 | 1.00 | no | 1.00 | 2,842 / 3,084 |
| 237 | SL_Z6_R20_T150_C2 | gas1_ics1 | 10 | +12 | 6,074 [5,567, 6,633] | 135.4 | 1.00 | no | 0.99 | 3,051 / 3,005 |
| 238 | SY_Z8_R30_T150_C2 | gas1_ics1 | 10 | +14 | 6,074 [5,539, 6,600] | 136.3 | 1.01 | no | 1.00 | 3,006 / 3,146 |
| 239 | SL_Z8_R20_T150_C5 | gas1_ics1 | 10 | +29 | 6,131 [5,617, 6,687] | 136.3 | 1.00 | no | 0.99 | 3,304 / 2,852 |
| 240 | SY_Z10_R20_T150_C5 | gas1_ics1 | 10 | +34 | 6,074 [5,561, 6,623] | 136.9 | 1.03 | no | 0.99 | 3,200 / 3,057 |
| 241 | SY_Z8_R20_T150_C2 | gas1_ics1 | 10 | +41 | 6,074 [5,571, 6,631] | 138.6 | 1.02 | no | 0.99 | 3,051 / 3,145 |
| 242 | SO_Z10_R30_T150_C2 | gas1_ics1 | 10 | +42 | 5,681 [5,126, 6,181] | 120.0 | 1.00 | no | 0.99 | 2,752 / 2,922 |
| 243 | SO_Z8_R30_T150_C2 | gas1_ics1 | 10 | +59 | 5,624 [5,117, 6,170] | 119.8 | 0.99 | no | 0.98 | 2,707 / 2,872 |
| 244 | SY_Z10_R30_T150_C2 | gas1_ics1 | 10 | +64 | 6,131 [5,621, 6,686] | 139.0 | 1.03 | no | 0.99 | 3,186 / 3,147 |
| 245 | SL_Z10_R30_T150_C5 | gas1_ics1 | 10 | +65 | 6,243 [5,689, 6,764] | 137.9 | 1.01 | no | 0.99 | 3,439 / 2,853 |
| 246 | SO_Z8_R20_T150_C2 | gas1_ics1 | 10 | +67 | 5,681 [5,131, 6,184] | 120.7 | 0.94 | no | 0.98 | 2,707 / 2,632 |
| 247 | SO_Z10_R20_T150_C2 | gas1_ics1 | 10 | +83 | 5,681 [5,158, 6,208] | 121.9 | 0.96 | no | 0.99 | 2,782 / 2,669 |
| 248 | SL_Z8_R30_T150_C2 | gas1_ics1 | 10 | +86 | 6,243 [5,723, 6,793] | 141.1 | 1.00 | no | 0.99 | 3,246 / 3,007 |
| 249 | SO_Z6_R20_T150_C2 | gas1_ics1 | 10 | +98 | 5,624 [5,073, 6,128] | 117.9 | 0.94 | no | 1.00 | 2,558 / 2,734 |
| 250 | SL_Z10_R20_T150_C5 | gas1_ics1 | 10 | +101 | 6,243 [5,695, 6,769] | 138.6 | 1.02 | no | 0.99 | 3,515 / 2,866 |
| 251 | SO_Z6_R30_T150_C2 | gas1_ics1 | 10 | +106 | 5,737 [5,052, 6,143] | 117.6 | 0.91 | no | 1.00 | 2,558 / 2,670 |
| 252 | SO_Z4_R20_T150_C2 | gas1_ics1 | 10 | +107 | 5,568 [5,042, 6,099] | 116.0 | 0.97 | no | 1.00 | 2,438 / 2,963 |
| 253 | SY_Z10_R20_T150_C2 | gas1_ics1 | 10 | +109 | 6,187 [5,648, 6,709] | 140.9 | 1.03 | no | 0.99 | 3,201 / 3,147 |
| 254 | SO_Z4_R30_T150_C2 | gas1_ics1 | 10 | +110 | 5,568 [5,046, 6,103] | 116.2 | 0.97 | no | 1.00 | 2,438 / 2,937 |
| 255 | SL_Z8_R20_T150_C2 | gas1_ics1 | 10 | +110 | 6,299 [5,748, 6,817] | 142.5 | 1.00 | no | 0.99 | 3,291 / 3,020 |
| 256 | SL_Z10_R30_T150_C2 | gas1_ics1 | 10 | +146 | 6,356 [5,816, 6,892] | 144.5 | 1.00 | no | 0.99 | 3,351 / 3,035 |
| 257 | SL_Z10_R20_T150_C2 | gas1_ics1 | 10 | +185 | 6,356 [5,833, 6,906] | 145.6 | 1.01 | no | 0.99 | 3,381 / 3,045 |
| 258 | gll_iem_v07 | - | 7 | +1,009 | 2,531 [1,997, 3,146] | 20.4 | 0.58 | yes | 0.96 | 265 / 1,207 |
