# Stage C — 2026-09-23

Consensus selection across 257 model rows (both) by the **fraction-of-rows rule**: a candidate is eligible if it clears its own BIC threshold (p·ln 1080 = 63 per template) in ≥ 80% of rows; the eligible candidate that is the largest improvement in the most rows is selected. Median Δ is shown for reference only. Selected: **bubble_totani → loopI_shells+nps**

| step | chosen | clears BIC in | best among eligible in | median Δ (ref.) | runner-up (clears / best) |
|---|---|---|---|---|---|
| 1 | bubble_totani | 100% | 100% | 28,766 | bubble_flat (100% / 0%) |
| 2 | loopI_shells+nps | 100% | 93% | 7,059 | loopI_shells (100% / 7%) |
| 3 | (stop) sun_moon_ic | 77% | 0% | 88 | loopII (72% / 0%) |

## Each model's own path

| model | variant | path |
|---|---|---|
| SL_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells |
| SL_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SL_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SL_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SL_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells → erosita_bright_NS → loopII |
| SL_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SO_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → erosita_unif_NS → sun_moon_ic |
| SO_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SO_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → erosita_unif_NS → sun_moon_ic |
| SO_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → erosita_unif_NS → sun_moon_ic |
| SO_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → erosita_unif_NS → sun_moon_ic |
| SO_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SO_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SO_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SS_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SS_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SS_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SS_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SS_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SS_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → sun_moon_ic |
| SS_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z10_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z10_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z10_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z10_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z10_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z10_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z10_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z10_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z10_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z10_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z10_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z4_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z4_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z4_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z4_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z4_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z4_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z4_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z4_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z4_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z4_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z4_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z4_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z4_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z4_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z6_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z6_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z6_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z6_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z6_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z6_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells → sun_moon_ic |
| SY_Z6_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z6_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z6_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z6_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z6_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z8_R20_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R20_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z8_R20_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R20_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z8_R20_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z8_R20_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z8_R20_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z8_R20_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z8_R30_T100000_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R30_T100000_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z8_R30_T100000_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R30_T100000_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps |
| SY_Z8_R30_T150_C2 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS |
| SY_Z8_R30_T150_C2 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| SY_Z8_R30_T150_C5 | gas1_ics1 | bubble_totani → loopI_shells+nps → erosita_bright_NS → loopII |
| SY_Z8_R30_T150_C5 | gas5_ics3 | bubble_totani → loopI_shells+nps → sun_moon_ic |
| a2020 | - | bubble_totani → loopI_shells+nps |

## Ranking with the consensus set (bubble_totani, loopI_shells+nps)

| rank | model | variant | n | Δ(−2lnL) vs Totani baseline | halo @ 20.8 | TS @ 20.8 | quad Σ/joint | any neg | inj min | N / S |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | a2020 | - | 21 | -1,589 | 8,549 [7,765, 9,104] | 43.2 | 0.79 | no | 1.00 | 3,915 / 2,842 |
| 2 | SS_Z10_R30_T150_C2 | gas5_ics3 | 15 | -1,405 | 6,243 [5,292, 7,175] | 20.0 | 0.70 | no | 0.99 | 1,789 / 2,557 |
| 3 | SS_Z10_R20_T150_C2 | gas5_ics3 | 15 | -1,400 | 6,243 [5,312, 7,180] | 20.1 | 0.71 | no | 0.98 | 1,848 / 2,611 |
| 4 | SY_Z10_R30_T150_C2 | gas5_ics3 | 15 | -1,396 | 6,187 [5,252, 7,133] | 18.6 | 0.69 | yes | 1.00 | 1,787 / 2,455 |
| 5 | SS_Z8_R30_T150_C2 | gas5_ics3 | 15 | -1,393 | 6,243 [5,359, 7,215] | 22.0 | 0.68 | no | 0.99 | 1,923 / 2,336 |
| 6 | SY_Z10_R20_T150_C2 | gas5_ics3 | 15 | -1,389 | 6,187 [5,267, 7,146] | 18.6 | 0.69 | yes | 1.00 | 1,832 / 2,467 |
| 7 | SS_Z8_R20_T150_C2 | gas5_ics3 | 15 | -1,389 | 6,299 [5,401, 7,246] | 22.1 | 0.69 | no | 0.99 | 1,983 / 2,349 |
| 8 | SY_Z10_R30_T150_C5 | gas5_ics3 | 15 | -1,388 | 6,412 [5,467, 7,329] | 19.5 | 0.80 | no | 1.00 | 2,698 / 2,419 |
| 9 | SS_Z10_R30_T150_C5 | gas5_ics3 | 15 | -1,384 | 6,524 [5,613, 7,486] | 24.7 | 0.87 | no | 0.98 | 2,806 / 2,855 |
| 10 | SY_Z10_R20_T150_C5 | gas5_ics3 | 15 | -1,378 | 6,412 [5,496, 7,353] | 19.6 | 0.81 | no | 0.97 | 2,757 / 2,432 |
| 11 | SY_Z8_R30_T150_C2 | gas5_ics3 | 15 | -1,378 | 6,299 [5,358, 7,222] | 20.9 | 0.70 | no | 0.99 | 1,937 / 2,463 |
| 12 | SS_Z10_R20_T150_C5 | gas5_ics3 | 15 | -1,378 | 6,524 [5,615, 7,478] | 24.5 | 0.87 | no | 0.99 | 2,850 / 2,842 |
| 13 | SY_Z8_R20_T150_C2 | gas5_ics3 | 15 | -1,374 | 6,299 [5,396, 7,256] | 21.0 | 0.71 | no | 0.99 | 2,011 / 2,476 |
| 14 | SS_Z8_R30_T150_C5 | gas5_ics3 | 15 | -1,372 | 6,637 [5,701, 7,542] | 26.7 | 0.90 | no | 0.99 | 3,060 / 2,920 |
| 15 | SL_Z10_R30_T150_C2 | gas5_ics3 | 15 | -1,369 | 6,299 [5,340, 7,237] | 19.6 | 0.69 | no | 1.00 | 2,026 / 2,342 |
| 16 | SY_Z8_R30_T150_C5 | gas5_ics3 | 15 | -1,369 | 6,468 [5,573, 7,410] | 21.5 | 0.84 | no | 0.99 | 2,952 / 2,458 |
| 17 | SS_Z8_R20_T150_C5 | gas5_ics3 | 15 | -1,367 | 6,637 [5,700, 7,535] | 26.5 | 0.91 | no | 0.99 | 3,089 / 2,920 |
| 18 | SS_Z6_R30_T150_C2 | gas5_ics3 | 15 | -1,363 | 6,468 [5,556, 7,345] | 26.9 | 0.73 | no | 1.00 | 2,613 / 2,085 |
| 19 | SY_Z8_R20_T150_C5 | gas5_ics3 | 15 | -1,362 | 6,524 [5,588, 7,423] | 21.6 | 0.84 | no | 0.99 | 2,996 / 2,470 |
| 20 | SS_Z6_R20_T150_C2 | gas5_ics3 | 15 | -1,362 | 6,468 [5,562, 7,349] | 26.9 | 0.73 | no | 1.00 | 2,613 / 2,098 |
| 21 | SL_Z10_R20_T150_C2 | gas5_ics3 | 15 | -1,361 | 6,299 [5,357, 7,252] | 19.7 | 0.70 | no | 1.00 | 2,056 / 2,355 |
| 22 | SL_Z8_R30_T150_C2 | gas5_ics3 | 15 | -1,349 | 6,412 [5,466, 7,344] | 22.1 | 0.72 | no | 0.98 | 2,205 / 2,419 |
| 23 | SS_Z6_R30_T150_C5 | gas5_ics3 | 15 | -1,349 | 6,693 [5,788, 7,540] | 30.4 | 0.98 | no | 1.00 | 3,734 / 2,821 |
| 24 | SL_Z10_R30_T150_C5 | gas5_ics3 | 15 | -1,348 | 6,524 [5,614, 7,484] | 21.1 | 0.77 | no | 0.96 | 2,758 / 2,295 |
| 25 | SS_Z6_R20_T150_C5 | gas5_ics3 | 15 | -1,347 | 6,637 [5,784, 7,536] | 30.1 | 0.99 | no | 1.00 | 3,734 / 2,834 |
| 26 | SL_Z8_R20_T150_C2 | gas5_ics3 | 15 | -1,344 | 6,412 [5,499, 7,375] | 22.1 | 0.73 | no | 0.98 | 2,250 / 2,419 |
| 27 | SY_Z6_R30_T150_C2 | gas5_ics3 | 15 | -1,341 | 6,468 [5,568, 7,396] | 26.1 | 0.76 | no | 1.00 | 2,611 / 2,318 |
| 28 | SY_Z6_R20_T150_C2 | gas5_ics3 | 15 | -1,339 | 6,468 [5,569, 7,396] | 26.1 | 0.76 | no | 0.98 | 2,596 / 2,331 |
| 29 | SO_Z10_R30_T150_C2 | gas5_ics3 | 15 | -1,339 | 6,187 [5,367, 7,076] | 29.6 | 0.81 | no | 0.98 | 2,271 / 2,717 |
| 30 | SY_Z10_R30_T100000_C2 | gas5_ics3 | 15 | -1,338 | 6,918 [5,951, 7,891] | 18.9 | 0.62 | yes | 0.99 | 2,201 / 2,102 |
| 31 | SL_Z10_R20_T150_C5 | gas5_ics3 | 15 | -1,337 | 6,580 [5,624, 7,491] | 21.2 | 0.78 | no | 0.98 | 2,802 / 2,308 |
| 32 | SO_Z10_R20_T150_C2 | gas5_ics3 | 15 | -1,337 | 6,187 [5,349, 7,068] | 29.3 | 0.82 | no | 0.99 | 2,286 / 2,758 |
| 33 | SO_Z8_R30_T150_C2 | gas5_ics3 | 15 | -1,332 | 6,299 [5,447, 7,191] | 31.8 | 0.79 | no | 0.99 | 2,495 / 2,471 |
| 34 | SY_Z6_R30_T150_C5 | gas5_ics3 | 15 | -1,331 | 6,580 [5,695, 7,480] | 25.9 | 0.91 | no | 1.00 | 3,492 / 2,476 |
| 35 | SO_Z8_R20_T150_C2 | gas5_ics3 | 15 | -1,329 | 6,299 [5,437, 7,183] | 31.6 | 0.79 | no | 0.99 | 2,495 / 2,499 |
| 36 | SY_Z6_R20_T150_C5 | gas5_ics3 | 15 | -1,328 | 6,580 [5,700, 7,485] | 25.9 | 0.91 | no | 0.99 | 3,492 / 2,488 |
| 37 | SY_Z8_R30_T100000_C2 | gas5_ics3 | 15 | -1,327 | 6,862 [5,937, 7,833] | 20.6 | 0.64 | yes | 0.99 | 2,275 / 2,149 |
| 38 | SO_Z8_R30_T150_C5 | gas5_ics3 | 15 | -1,325 | 6,468 [5,611, 7,318] | 33.6 | 0.96 | no | 0.99 | 3,272 / 2,919 |
| 39 | SY_Z10_R20_T100000_C2 | gas5_ics3 | 15 | -1,324 | 6,974 [6,006, 7,934] | 19.2 | 0.62 | yes | 0.99 | 2,230 / 2,128 |
| 40 | SL_Z8_R30_T150_C5 | gas5_ics3 | 15 | -1,322 | 6,637 [5,739, 7,580] | 23.4 | 0.80 | no | 0.99 | 2,967 / 2,321 |
| 41 | SO_Z10_R30_T150_C5 | gas5_ics3 | 15 | -1,322 | 6,412 [5,633, 7,269] | 33.0 | 0.92 | no | 0.99 | 2,988 / 2,918 |
| 42 | SY_Z8_R20_T100000_C2 | gas5_ics3 | 15 | -1,320 | 6,918 [5,960, 7,850] | 20.7 | 0.64 | yes | 0.98 | 2,275 / 2,175 |
| 43 | SS_Z8_R30_T100000_C2 | gas5_ics3 | 15 | -1,319 | 7,087 [6,148, 8,015] | 27.0 | 0.77 | no | 0.98 | 2,669 / 2,782 |
| 44 | SO_Z8_R20_T150_C5 | gas5_ics3 | 15 | -1,318 | 6,468 [5,615, 7,322] | 33.7 | 0.96 | no | 0.99 | 3,302 / 2,932 |
| 45 | SO_Z6_R30_T150_C5 | gas5_ics3 | 15 | -1,316 | 6,524 [5,646, 7,319] | 35.4 | 1.01 | no | 0.99 | 3,825 / 2,744 |
| 46 | SL_Z8_R20_T150_C5 | gas5_ics3 | 15 | -1,316 | 6,693 [5,744, 7,584] | 23.5 | 0.80 | no | 0.99 | 2,997 / 2,333 |
| 47 | SS_Z10_R30_T100000_C2 | gas5_ics3 | 15 | -1,316 | 7,143 [6,183, 8,079] | 25.5 | 0.75 | no | 0.99 | 2,610 / 2,718 |
| 48 | SO_Z6_R20_T150_C5 | gas5_ics3 | 15 | -1,314 | 6,524 [5,653, 7,325] | 35.5 | 1.01 | no | 0.99 | 3,840 / 2,744 |
| 49 | SO_Z10_R20_T150_C5 | gas5_ics3 | 15 | -1,312 | 6,468 [5,647, 7,284] | 33.4 | 0.91 | no | 0.98 | 2,973 / 2,918 |
| 50 | SS_Z8_R20_T100000_C2 | gas5_ics3 | 15 | -1,311 | 7,087 [6,175, 8,035] | 27.3 | 0.78 | no | 0.99 | 2,699 / 2,795 |
| 51 | SO_Z6_R30_T150_C2 | gas5_ics3 | 15 | -1,309 | 6,468 [5,608, 7,353] | 35.9 | 0.80 | no | 1.00 | 3,050 / 2,118 |
| 52 | SO_Z6_R20_T150_C2 | gas5_ics3 | 15 | -1,308 | 6,468 [5,605, 7,351] | 35.9 | 0.80 | no | 1.00 | 3,035 / 2,131 |
| 53 | SL_Z6_R30_T150_C2 | gas5_ics3 | 15 | -1,308 | 6,580 [5,661, 7,488] | 27.5 | 0.78 | no | 1.00 | 2,775 / 2,329 |
| 54 | SL_Z6_R20_T150_C2 | gas5_ics3 | 15 | -1,306 | 6,580 [5,667, 7,493] | 27.5 | 0.78 | no | 1.00 | 2,775 / 2,341 |
| 55 | SS_Z10_R20_T100000_C2 | gas5_ics3 | 15 | -1,304 | 7,199 [6,239, 8,126] | 26.1 | 0.75 | no | 1.00 | 2,684 / 2,744 |
| 56 | SL_Z10_R30_T100000_C2 | gas5_ics3 | 15 | -1,303 | 7,087 [6,113, 8,001] | 19.5 | 0.57 | yes | 0.99 | 2,124 / 1,927 |
| 57 | SS_Z4_R30_T150_C2 | gas5_ics3 | 15 | -1,300 | 6,862 [6,009, 7,534] | 41.2 | 0.97 | no | 0.99 | 3,887 / 2,767 |
| 58 | SS_Z4_R20_T150_C2 | gas5_ics3 | 15 | -1,300 | 6,862 [6,012, 7,538] | 41.2 | 0.97 | no | 0.99 | 3,902 / 2,767 |
| 59 | SO_Z8_R30_T100000_C2 | gas5_ics3 | 15 | -1,293 | 6,805 [5,958, 7,683] | 31.4 | 0.81 | no | 1.00 | 2,776 / 2,769 |
| 60 | SS_Z6_R30_T100000_C2 | gas5_ics3 | 15 | -1,292 | 7,143 [6,258, 8,049] | 31.2 | 0.83 | no | 0.98 | 3,281 / 2,643 |
| 61 | SL_Z10_R20_T100000_C2 | gas5_ics3 | 15 | -1,290 | 7,087 [6,149, 8,035] | 19.6 | 0.58 | yes | 0.98 | 2,154 / 1,939 |
| 62 | SS_Z6_R20_T100000_C2 | gas5_ics3 | 15 | -1,288 | 7,143 [6,271, 8,059] | 31.3 | 0.84 | no | 1.00 | 3,296 / 2,669 |
| 63 | SL_Z8_R30_T100000_C2 | gas5_ics3 | 15 | -1,287 | 7,087 [6,130, 7,986] | 21.8 | 0.59 | yes | 0.98 | 2,154 / 2,039 |
| 64 | SO_Z10_R30_T100000_C2 | gas5_ics3 | 15 | -1,285 | 6,862 [5,991, 7,734] | 30.8 | 0.79 | no | 0.99 | 2,656 / 2,783 |
| 65 | SO_Z8_R20_T100000_C2 | gas5_ics3 | 15 | -1,284 | 6,862 [5,997, 7,717] | 31.6 | 0.82 | no | 0.99 | 2,820 / 2,783 |
| 66 | SY_Z6_R30_T100000_C2 | gas5_ics3 | 15 | -1,283 | 7,030 [6,096, 7,956] | 24.7 | 0.70 | yes | 1.00 | 2,843 / 2,060 |
| 67 | SY_Z6_R20_T100000_C2 | gas5_ics3 | 15 | -1,280 | 7,030 [6,104, 7,963] | 24.7 | 0.70 | yes | 1.00 | 2,828 / 2,060 |
| 68 | SL_Z8_R20_T100000_C2 | gas5_ics3 | 15 | -1,280 | 7,087 [6,145, 7,998] | 21.9 | 0.60 | yes | 0.99 | 2,169 / 2,052 |
| 69 | SO_Z6_R30_T100000_C2 | gas5_ics3 | 15 | -1,276 | 6,918 [6,051, 7,741] | 33.8 | 0.84 | no | 0.99 | 3,343 / 2,460 |
| 70 | SO_Z10_R20_T100000_C2 | gas5_ics3 | 15 | -1,276 | 6,862 [6,016, 7,749] | 31.2 | 0.79 | no | 0.99 | 2,656 / 2,783 |
| 71 | SS_Z4_R20_T150_C5 | gas5_ics3 | 15 | -1,276 | 6,974 [6,110, 7,588] | 43.0 | 1.03 | no | 0.98 | 4,260 / 2,911 |
| 72 | SS_Z4_R30_T150_C5 | gas5_ics3 | 15 | -1,276 | 6,974 [6,108, 7,585] | 43.0 | 1.03 | no | 0.98 | 4,260 / 2,911 |
| 73 | SL_Z6_R30_T150_C5 | gas5_ics3 | 15 | -1,275 | 6,749 [5,871, 7,637] | 28.3 | 0.87 | no | 1.00 | 3,447 / 2,391 |
| 74 | SY_Z10_R30_T100000_C5 | gas5_ics3 | 15 | -1,274 | 7,199 [6,270, 8,177] | 24.2 | 0.70 | no | 0.98 | 2,741 / 2,279 |
| 75 | SO_Z6_R20_T100000_C2 | gas5_ics3 | 15 | -1,274 | 6,918 [6,055, 7,745] | 33.9 | 0.85 | no | 1.00 | 3,358 / 2,500 |
| 76 | SL_Z6_R20_T150_C5 | gas5_ics3 | 15 | -1,272 | 6,805 [5,879, 7,643] | 28.3 | 0.86 | no | 0.98 | 3,447 / 2,391 |
| 77 | SY_Z4_R20_T150_C2 | gas5_ics3 | 15 | -1,267 | 6,862 [6,014, 7,596] | 39.1 | 0.92 | no | 1.00 | 3,766 / 2,567 |
| 78 | SO_Z10_R30_T100000_C5 | gas5_ics3 | 15 | -1,267 | 6,805 [5,978, 7,703] | 31.2 | 0.81 | no | 0.99 | 2,850 / 2,693 |
| 79 | SY_Z8_R30_T100000_C5 | gas5_ics3 | 15 | -1,267 | 7,143 [6,231, 8,097] | 25.5 | 0.72 | no | 0.98 | 2,831 / 2,327 |
| 80 | SY_Z4_R30_T150_C2 | gas5_ics3 | 15 | -1,267 | 6,862 [6,013, 7,595] | 39.2 | 0.92 | no | 1.00 | 3,766 / 2,567 |
| 81 | SO_Z8_R30_T100000_C5 | gas5_ics3 | 15 | -1,267 | 6,862 [5,995, 7,697] | 31.6 | 0.85 | no | 0.99 | 3,149 / 2,652 |
| 82 | SY_Z10_R20_T100000_C5 | gas5_ics3 | 15 | -1,261 | 7,199 [6,297, 8,184] | 24.3 | 0.70 | yes | 0.99 | 2,770 / 2,292 |
| 83 | SO_Z4_R30_T150_C5 | gas5_ics3 | 15 | -1,260 | 6,693 [5,934, 7,211] | 45.4 | 0.99 | no | 0.99 | 4,020 / 2,601 |
| 84 | SO_Z4_R20_T150_C5 | gas5_ics3 | 15 | -1,259 | 6,693 [5,936, 7,211] | 45.4 | 0.99 | no | 0.99 | 4,020 / 2,601 |
| 85 | SO_Z8_R20_T100000_C5 | gas5_ics3 | 15 | -1,257 | 6,862 [6,032, 7,728] | 31.9 | 0.85 | no | 0.99 | 3,164 / 2,640 |
| 86 | SY_Z4_R20_T150_C5 | gas5_ics3 | 15 | -1,257 | 6,862 [6,024, 7,560] | 37.2 | 0.98 | no | 0.98 | 4,260 / 2,469 |
| 87 | SY_Z4_R30_T150_C5 | gas5_ics3 | 15 | -1,256 | 6,862 [6,022, 7,559] | 37.2 | 0.98 | no | 0.98 | 4,260 / 2,469 |
| 88 | SO_Z10_R20_T100000_C5 | gas5_ics3 | 15 | -1,256 | 6,862 [6,006, 7,722] | 31.6 | 0.81 | no | 0.99 | 2,864 / 2,706 |
| 89 | SO_Z6_R30_T100000_C5 | gas5_ics3 | 15 | -1,255 | 6,862 [6,019, 7,641] | 33.2 | 0.87 | no | 1.00 | 3,434 / 2,518 |
| 90 | SL_Z10_R30_T100000_C5 | gas5_ics3 | 15 | -1,254 | 7,199 [6,247, 8,097] | 22.2 | 0.61 | yes | 0.98 | 2,514 / 1,897 |
| 91 | SO_Z6_R20_T100000_C5 | gas5_ics3 | 15 | -1,252 | 6,862 [6,024, 7,647] | 33.3 | 0.87 | no | 0.99 | 3,434 / 2,518 |
| 92 | SY_Z8_R20_T100000_C5 | gas5_ics3 | 15 | -1,252 | 7,199 [6,297, 8,157] | 25.8 | 0.71 | no | 0.99 | 2,815 / 2,314 |
| 93 | SO_Z4_R30_T150_C2 | gas5_ics3 | 15 | -1,252 | 6,749 [5,975, 7,295] | 48.1 | 0.96 | no | 0.98 | 3,856 / 2,627 |
| 94 | SO_Z4_R20_T150_C2 | gas5_ics3 | 15 | -1,252 | 6,749 [5,975, 7,295] | 48.1 | 0.96 | no | 0.98 | 3,856 / 2,627 |
| 95 | SS_Z8_R30_T100000_C5 | gas5_ics3 | 15 | -1,251 | 7,312 [6,383, 8,224] | 32.7 | 0.84 | no | 0.99 | 3,284 / 2,883 |
| 96 | SS_Z10_R30_T100000_C5 | gas5_ics3 | 15 | -1,249 | 7,368 [6,445, 8,310] | 31.7 | 0.81 | no | 0.99 | 3,135 / 2,806 |
| 97 | SL_Z10_R20_T100000_C5 | gas5_ics3 | 15 | -1,240 | 7,199 [6,265, 8,111] | 22.2 | 0.61 | yes | 0.98 | 2,514 / 1,897 |
| 98 | SS_Z8_R20_T100000_C5 | gas5_ics3 | 15 | -1,240 | 7,312 [6,420, 8,254] | 32.9 | 0.85 | no | 0.99 | 3,343 / 2,870 |
| 99 | SL_Z8_R30_T100000_C5 | gas5_ics3 | 15 | -1,240 | 7,143 [6,217, 8,042] | 23.8 | 0.63 | yes | 0.99 | 2,544 / 1,958 |
| 100 | SS_Z10_R20_T100000_C5 | gas5_ics3 | 15 | -1,238 | 7,368 [6,462, 8,322] | 31.9 | 0.82 | no | 0.98 | 3,194 / 2,833 |
| 101 | SS_Z6_R30_T100000_C5 | gas5_ics3 | 15 | -1,237 | 7,255 [6,364, 8,130] | 34.3 | 0.87 | no | 0.98 | 3,552 / 2,790 |
| 102 | SS_Z6_R20_T100000_C5 | gas5_ics3 | 15 | -1,232 | 7,255 [6,374, 8,138] | 34.4 | 0.88 | no | 0.99 | 3,537 / 2,816 |
| 103 | SL_Z6_R30_T100000_C2 | gas5_ics3 | 15 | -1,231 | 7,199 [6,297, 8,126] | 26.4 | 0.66 | yes | 1.00 | 2,767 / 2,018 |
| 104 | SL_Z6_R20_T100000_C2 | gas5_ics3 | 15 | -1,227 | 7,199 [6,303, 8,131] | 26.4 | 0.66 | yes | 0.99 | 2,752 / 2,018 |
| 105 | SY_Z6_R30_T100000_C5 | gas5_ics3 | 15 | -1,226 | 7,255 [6,351, 8,125] | 29.3 | 0.75 | no | 0.98 | 3,145 / 2,323 |
| 106 | SY_Z6_R20_T100000_C5 | gas5_ics3 | 15 | -1,223 | 7,255 [6,352, 8,124] | 29.2 | 0.75 | no | 1.00 | 3,145 / 2,323 |
| 107 | SL_Z8_R20_T100000_C5 | gas5_ics3 | 15 | -1,221 | 7,199 [6,305, 8,121] | 24.2 | 0.62 | yes | 0.99 | 2,544 / 1,946 |
| 108 | SL_Z4_R30_T150_C2 | gas5_ics3 | 15 | -1,217 | 6,974 [6,094, 7,626] | 41.0 | 0.87 | no | 1.00 | 3,600 / 2,456 |
| 109 | SL_Z4_R20_T150_C2 | gas5_ics3 | 15 | -1,217 | 6,974 [6,099, 7,631] | 41.0 | 0.87 | no | 1.00 | 3,600 / 2,456 |
| 110 | SS_Z4_R30_T100000_C2 | gas5_ics3 | 15 | -1,212 | 7,199 [6,391, 7,973] | 42.9 | 0.95 | no | 1.00 | 4,077 / 2,752 |
| 111 | SO_Z4_R30_T100000_C2 | gas5_ics3 | 15 | -1,211 | 7,030 [6,180, 7,573] | 43.3 | 0.95 | no | 1.00 | 4,154 / 2,522 |
| 112 | SO_Z4_R20_T100000_C2 | gas5_ics3 | 15 | -1,211 | 7,030 [6,181, 7,574] | 43.3 | 0.95 | no | 1.00 | 4,154 / 2,522 |
| 113 | SS_Z4_R20_T100000_C2 | gas5_ics3 | 15 | -1,210 | 7,255 [6,398, 7,981] | 43.0 | 0.94 | no | 1.00 | 4,077 / 2,765 |
| 114 | SY_Z4_R30_T100000_C2 | gas5_ics3 | 15 | -1,209 | 7,143 [6,258, 7,974] | 35.7 | 0.84 | no | 1.00 | 3,701 / 2,273 |
| 115 | SY_Z4_R20_T100000_C2 | gas5_ics3 | 15 | -1,209 | 7,143 [6,262, 7,978] | 35.8 | 0.84 | no | 0.99 | 3,701 / 2,273 |
| 116 | SL_Z6_R30_T100000_C5 | gas5_ics3 | 15 | -1,182 | 7,255 [6,367, 8,124] | 28.1 | 0.67 | yes | 1.00 | 2,874 / 1,993 |
| 117 | SL_Z4_R30_T150_C5 | gas5_ics3 | 15 | -1,180 | 7,030 [6,193, 7,681] | 40.1 | 0.96 | no | 0.98 | 4,334 / 2,398 |
| 118 | SL_Z4_R20_T150_C5 | gas5_ics3 | 15 | -1,180 | 7,030 [6,193, 7,681] | 40.1 | 0.96 | no | 0.98 | 4,334 / 2,410 |
| 119 | SL_Z6_R20_T100000_C5 | gas5_ics3 | 15 | -1,179 | 7,255 [6,367, 8,124] | 28.0 | 0.67 | yes | 0.98 | 2,859 / 1,993 |
| 120 | SS_Z4_R30_T100000_C5 | gas5_ics3 | 15 | -1,175 | 7,143 [6,342, 7,876] | 41.2 | 0.97 | no | 1.00 | 4,258 / 2,677 |
| 121 | SO_Z4_R20_T100000_C5 | gas5_ics3 | 15 | -1,175 | 6,974 [6,205, 7,506] | 41.5 | 0.94 | no | 1.00 | 4,081 / 2,474 |
| 122 | SO_Z4_R30_T100000_C5 | gas5_ics3 | 15 | -1,174 | 6,974 [6,202, 7,504] | 41.4 | 0.94 | no | 1.00 | 4,081 / 2,474 |
| 123 | SS_Z4_R20_T100000_C5 | gas5_ics3 | 15 | -1,166 | 7,255 [6,401, 7,920] | 41.2 | 0.94 | no | 1.00 | 4,183 / 2,666 |
| 124 | SY_Z4_R30_T100000_C5 | gas5_ics3 | 15 | -1,130 | 7,312 [6,467, 8,086] | 38.6 | 0.85 | no | 0.98 | 3,852 / 2,354 |
| 125 | SY_Z4_R20_T100000_C5 | gas5_ics3 | 15 | -1,130 | 7,312 [6,468, 8,088] | 38.6 | 0.85 | no | 0.98 | 3,852 / 2,354 |
| 126 | SL_Z4_R30_T100000_C2 | gas5_ics3 | 15 | -1,129 | 7,312 [6,469, 8,158] | 38.4 | 0.77 | no | 1.00 | 3,445 / 2,200 |
| 127 | SL_Z4_R20_T100000_C2 | gas5_ics3 | 15 | -1,129 | 7,312 [6,470, 8,159] | 38.4 | 0.77 | no | 0.98 | 3,445 / 2,200 |
| 128 | SL_Z4_R30_T100000_C5 | gas5_ics3 | 15 | -1,076 | 7,368 [6,460, 8,098] | 36.7 | 0.79 | no | 0.98 | 3,657 / 2,152 |
| 129 | SL_Z4_R20_T100000_C5 | gas5_ics3 | 15 | -1,075 | 7,368 [6,459, 8,098] | 36.7 | 0.79 | no | 0.98 | 3,657 / 2,152 |
| 130 | SS_Z4_R30_T100000_C5 | gas1_ics1 | 9 | -373 | 6,018 [5,529, 6,514] | 155.5 | 0.90 | no | 1.00 | 2,603 / 2,794 |
| 131 | SS_Z4_R20_T100000_C5 | gas1_ics1 | 9 | -366 | 6,018 [5,553, 6,534] | 158.1 | 0.90 | no | 1.00 | 2,633 / 2,807 |
| 132 | SY_Z4_R30_T100000_C5 | gas1_ics1 | 9 | -355 | 6,131 [5,657, 6,633] | 163.7 | 0.88 | no | 1.00 | 2,721 / 2,671 |
| 133 | SY_Z4_R20_T100000_C5 | gas1_ics1 | 9 | -351 | 6,131 [5,653, 6,629] | 163.5 | 0.88 | no | 1.00 | 2,721 / 2,671 |
| 134 | SL_Z4_R30_T100000_C5 | gas1_ics1 | 9 | -314 | 6,187 [5,695, 6,672] | 165.0 | 0.89 | no | 1.00 | 2,960 / 2,519 |
| 135 | SS_Z6_R30_T100000_C5 | gas1_ics1 | 9 | -313 | 6,074 [5,588, 6,568] | 159.9 | 0.89 | no | 1.00 | 2,602 / 2,822 |
| 136 | SL_Z4_R20_T100000_C5 | gas1_ics1 | 9 | -313 | 6,187 [5,689, 6,666] | 164.7 | 0.89 | no | 1.00 | 2,960 / 2,519 |
| 137 | SS_Z6_R20_T100000_C5 | gas1_ics1 | 9 | -291 | 6,074 [5,592, 6,569] | 160.7 | 0.90 | no | 0.98 | 2,617 / 2,822 |
| 138 | SS_Z4_R30_T150_C5 | gas1_ics1 | 9 | -288 | 5,962 [5,477, 6,463] | 153.1 | 0.90 | no | 1.00 | 2,363 / 3,008 |
| 139 | SY_Z4_R30_T150_C5 | gas1_ics1 | 9 | -284 | 6,018 [5,526, 6,508] | 154.6 | 0.89 | no | 1.00 | 2,527 / 2,859 |
| 140 | SS_Z4_R20_T150_C5 | gas1_ics1 | 9 | -283 | 5,962 [5,474, 6,461] | 153.1 | 0.90 | no | 0.98 | 2,363 / 3,008 |
| 141 | SY_Z4_R20_T150_C5 | gas1_ics1 | 9 | -282 | 6,018 [5,523, 6,505] | 154.4 | 0.89 | no | 1.00 | 2,527 / 2,859 |
| 142 | SY_Z4_R30_T100000_C2 | gas1_ics1 | 9 | -280 | 6,243 [5,754, 6,732] | 168.3 | 0.88 | no | 1.00 | 2,647 / 2,840 |
| 143 | SY_Z4_R20_T100000_C2 | gas1_ics1 | 9 | -276 | 6,243 [5,751, 6,729] | 168.2 | 0.88 | no | 1.00 | 2,647 / 2,840 |
| 144 | SS_Z8_R30_T100000_C5 | gas1_ics1 | 9 | -272 | 6,131 [5,619, 6,601] | 160.7 | 0.89 | no | 0.99 | 2,632 / 2,811 |
| 145 | SY_Z6_R30_T100000_C5 | gas1_ics1 | 9 | -271 | 6,299 [5,781, 6,766] | 167.5 | 0.89 | no | 1.00 | 2,810 / 2,791 |
| 146 | SS_Z4_R30_T100000_C2 | gas1_ics1 | 9 | -266 | 6,074 [5,598, 6,580] | 160.1 | 0.90 | no | 0.99 | 2,558 / 2,899 |
| 147 | SS_Z4_R20_T100000_C2 | gas1_ics1 | 9 | -259 | 6,074 [5,597, 6,579] | 160.2 | 0.90 | no | 1.00 | 2,558 / 2,899 |
| 148 | SY_Z6_R20_T100000_C5 | gas1_ics1 | 9 | -253 | 6,299 [5,781, 6,765] | 168.0 | 0.89 | no | 1.00 | 2,840 / 2,791 |
| 149 | SL_Z4_R30_T100000_C2 | gas1_ics1 | 9 | -249 | 6,356 [5,839, 6,818] | 172.6 | 0.89 | no | 1.00 | 2,945 / 2,726 |
| 150 | SL_Z4_R20_T100000_C2 | gas1_ics1 | 9 | -249 | 6,299 [5,833, 6,812] | 172.2 | 0.90 | no | 1.00 | 2,931 / 2,713 |
| 151 | SL_Z6_R30_T100000_C5 | gas1_ics1 | 9 | -247 | 6,356 [5,870, 6,863] | 168.9 | 0.91 | no | 0.99 | 3,124 / 2,651 |
| 152 | SL_Z4_R30_T150_C5 | gas1_ics1 | 9 | -238 | 6,074 [5,585, 6,565] | 157.4 | 0.92 | no | 0.99 | 2,856 / 2,706 |
| 153 | SL_Z4_R20_T150_C5 | gas1_ics1 | 9 | -237 | 6,074 [5,580, 6,561] | 157.2 | 0.92 | no | 0.99 | 2,856 / 2,706 |
| 154 | SL_Z6_R20_T100000_C5 | gas1_ics1 | 9 | -236 | 6,356 [5,868, 6,860] | 169.0 | 0.91 | no | 0.99 | 3,124 / 2,651 |
| 155 | SO_Z4_R20_T100000_C5 | gas1_ics1 | 9 | -234 | 5,737 [5,257, 6,230] | 145.1 | 0.89 | no | 0.99 | 2,543 / 2,543 |
| 156 | SO_Z4_R30_T100000_C5 | gas1_ics1 | 9 | -232 | 5,737 [5,261, 6,234] | 145.3 | 0.89 | no | 1.00 | 2,543 / 2,543 |
| 157 | SS_Z10_R30_T100000_C5 | gas1_ics1 | 9 | -231 | 6,187 [5,692, 6,669] | 165.5 | 0.89 | no | 0.99 | 2,737 / 2,798 |
| 158 | SS_Z8_R20_T100000_C5 | gas1_ics1 | 9 | -230 | 6,131 [5,645, 6,621] | 163.9 | 0.90 | no | 1.00 | 2,677 / 2,823 |
| 159 | SO_Z10_R30_T100000_C5 | gas1_ics1 | 9 | -226 | 5,906 [5,401, 6,373] | 152.8 | 0.89 | no | 0.99 | 2,647 / 2,596 |
| 160 | SO_Z8_R30_T100000_C5 | gas1_ics1 | 9 | -226 | 5,849 [5,382, 6,353] | 152.0 | 0.89 | no | 0.98 | 2,602 / 2,595 |
| 161 | SS_Z6_R30_T150_C5 | gas1_ics1 | 9 | -220 | 5,962 [5,472, 6,455] | 152.8 | 0.92 | no | 1.00 | 2,513 / 2,969 |
| 162 | SO_Z6_R30_T100000_C5 | gas1_ics1 | 9 | -219 | 5,793 [5,331, 6,304] | 149.0 | 0.88 | no | 0.99 | 2,543 / 2,581 |
| 163 | SS_Z6_R30_T100000_C2 | gas1_ics1 | 9 | -214 | 6,131 [5,663, 6,640] | 164.9 | 0.91 | no | 1.00 | 2,647 / 2,941 |
| 164 | SY_Z8_R30_T100000_C5 | gas1_ics1 | 9 | -212 | 6,356 [5,859, 6,851] | 169.0 | 0.91 | no | 1.00 | 2,930 / 2,833 |
| 165 | SO_Z6_R20_T100000_C5 | gas1_ics1 | 9 | -208 | 5,793 [5,334, 6,307] | 149.5 | 0.88 | no | 0.99 | 2,543 / 2,581 |
| 166 | SL_Z8_R30_T100000_C5 | gas1_ics1 | 9 | -206 | 6,468 [5,976, 6,982] | 170.2 | 0.92 | no | 0.99 | 3,244 / 2,706 |
| 167 | SS_Z6_R20_T150_C5 | gas1_ics1 | 9 | -201 | 5,962 [5,476, 6,457] | 153.5 | 0.92 | no | 1.00 | 2,527 / 2,983 |
| 168 | SO_Z8_R20_T100000_C5 | gas1_ics1 | 9 | -200 | 5,906 [5,403, 6,370] | 154.2 | 0.89 | no | 0.99 | 2,632 / 2,609 |
| 169 | SS_Z6_R20_T100000_C2 | gas1_ics1 | 9 | -189 | 6,131 [5,669, 6,644] | 165.9 | 0.91 | no | 1.00 | 2,647 / 2,954 |
| 170 | SY_Z4_R30_T150_C2 | gas1_ics1 | 9 | -182 | 6,131 [5,616, 6,600] | 159.5 | 0.90 | no | 1.00 | 2,498 / 3,001 |
| 171 | SO_Z10_R20_T100000_C5 | gas1_ics1 | 9 | -181 | 5,906 [5,428, 6,395] | 155.4 | 0.89 | no | 0.99 | 2,662 / 2,622 |
| 172 | SS_Z8_R30_T150_C5 | gas1_ics1 | 9 | -180 | 6,018 [5,555, 6,535] | 158.0 | 0.94 | no | 0.97 | 2,692 / 2,970 |
| 173 | SY_Z4_R20_T150_C2 | gas1_ics1 | 9 | -179 | 6,131 [5,613, 6,597] | 159.4 | 0.90 | no | 1.00 | 2,498 / 3,001 |
| 174 | SS_Z4_R30_T150_C2 | gas1_ics1 | 9 | -177 | 6,018 [5,506, 6,494] | 155.0 | 0.90 | no | 0.98 | 2,349 / 3,086 |
| 175 | SS_Z10_R20_T100000_C5 | gas1_ics1 | 9 | -173 | 6,187 [5,715, 6,687] | 168.9 | 0.90 | no | 0.98 | 2,766 / 2,823 |
| 176 | SS_Z8_R30_T100000_C2 | gas1_ics1 | 9 | -173 | 6,187 [5,700, 6,679] | 166.0 | 0.91 | no | 1.00 | 2,722 / 2,889 |
| 177 | SS_Z4_R20_T150_C2 | gas1_ics1 | 9 | -171 | 6,018 [5,504, 6,491] | 155.0 | 0.90 | no | 0.98 | 2,349 / 3,086 |
| 178 | SY_Z6_R30_T100000_C2 | gas1_ics1 | 9 | -166 | 6,412 [5,915, 6,896] | 176.7 | 0.91 | no | 1.00 | 2,901 / 2,934 |
| 179 | SY_Z8_R20_T100000_C5 | gas1_ics1 | 9 | -164 | 6,356 [5,890, 6,877] | 172.8 | 0.92 | no | 0.99 | 2,975 / 2,859 |
| 180 | SL_Z8_R20_T100000_C5 | gas1_ics1 | 9 | -162 | 6,524 [6,009, 7,010] | 173.8 | 0.92 | no | 0.99 | 3,304 / 2,719 |
| 181 | SL_Z6_R30_T100000_C2 | gas1_ics1 | 9 | -161 | 6,524 [6,053, 7,045] | 179.8 | 0.93 | no | 0.99 | 3,230 / 2,808 |
| 182 | SY_Z10_R30_T100000_C5 | gas1_ics1 | 9 | -156 | 6,468 [5,967, 6,959] | 174.6 | 0.92 | no | 0.99 | 3,065 / 2,875 |
| 183 | SL_Z10_R30_T100000_C5 | gas1_ics1 | 9 | -150 | 6,580 [6,098, 7,109] | 175.5 | 0.93 | no | 0.99 | 3,379 / 2,748 |
| 184 | SL_Z6_R20_T100000_C2 | gas1_ics1 | 9 | -149 | 6,524 [6,054, 7,044] | 180.1 | 0.93 | no | 0.99 | 3,245 / 2,808 |
| 185 | SY_Z6_R20_T100000_C2 | gas1_ics1 | 9 | -148 | 6,412 [5,917, 6,897] | 177.4 | 0.91 | no | 1.00 | 2,916 / 2,947 |
| 186 | SY_Z6_R30_T150_C5 | gas1_ics1 | 9 | -147 | 6,131 [5,647, 6,632] | 160.6 | 0.93 | no | 1.00 | 2,766 / 2,964 |
| 187 | SS_Z8_R20_T150_C5 | gas1_ics1 | 9 | -144 | 6,074 [5,571, 6,546] | 159.8 | 0.94 | no | 1.00 | 2,722 / 2,970 |
| 188 | SL_Z4_R30_T150_C2 | gas1_ics1 | 9 | -140 | 6,243 [5,724, 6,707] | 164.8 | 0.90 | no | 0.99 | 2,796 / 2,836 |
| 189 | SL_Z4_R20_T150_C2 | gas1_ics1 | 9 | -139 | 6,187 [5,721, 6,703] | 164.6 | 0.91 | no | 0.99 | 2,796 / 2,836 |
| 190 | SS_Z10_R30_T100000_C2 | gas1_ics1 | 9 | -137 | 6,243 [5,778, 6,755] | 170.9 | 0.92 | no | 0.97 | 2,842 / 2,889 |
| 191 | SS_Z10_R30_T150_C5 | gas1_ics1 | 9 | -136 | 6,074 [5,598, 6,576] | 159.9 | 0.94 | no | 0.98 | 2,782 / 2,918 |
| 192 | SY_Z6_R20_T150_C5 | gas1_ics1 | 9 | -132 | 6,131 [5,648, 6,631] | 161.0 | 0.93 | no | 0.99 | 2,766 / 2,964 |
| 193 | SS_Z8_R20_T100000_C2 | gas1_ics1 | 9 | -131 | 6,187 [5,717, 6,692] | 168.3 | 0.92 | no | 1.00 | 2,752 / 2,915 |
| 194 | SY_Z8_R30_T100000_C2 | gas1_ics1 | 9 | -126 | 6,468 [5,990, 6,981] | 177.1 | 0.93 | no | 1.00 | 3,080 / 2,949 |
| 195 | SL_Z8_R30_T100000_C2 | gas1_ics1 | 9 | -120 | 6,637 [6,155, 7,160] | 180.5 | 0.93 | no | 0.99 | 3,364 / 2,836 |
| 196 | SO_Z10_R30_T100000_C2 | gas1_ics1 | 9 | -117 | 5,962 [5,483, 6,452] | 157.7 | 0.91 | no | 0.98 | 2,737 / 2,675 |
| 197 | SS_Z6_R30_T150_C2 | gas1_ics1 | 9 | -114 | 6,018 [5,503, 6,486] | 154.6 | 0.93 | no | 1.00 | 2,528 / 3,061 |
| 198 | SO_Z8_R30_T100000_C2 | gas1_ics1 | 9 | -109 | 5,906 [5,429, 6,400] | 154.2 | 0.91 | no | 0.98 | 2,662 / 2,687 |
| 199 | SL_Z10_R20_T100000_C5 | gas1_ics1 | 9 | -106 | 6,637 [6,107, 7,114] | 177.3 | 0.93 | no | 0.99 | 3,409 / 2,761 |
| 200 | SO_Z6_R30_T100000_C2 | gas1_ics1 | 9 | -104 | 5,906 [5,397, 6,368] | 153.2 | 0.90 | no | 0.99 | 2,602 / 2,686 |
| 201 | SY_Z10_R20_T100000_C5 | gas1_ics1 | 9 | -103 | 6,468 [5,981, 6,968] | 177.1 | 0.93 | no | 0.99 | 3,095 / 2,901 |
| 202 | SL_Z8_R20_T100000_C2 | gas1_ics1 | 9 | -93 | 6,637 [6,160, 7,163] | 181.6 | 0.94 | no | 0.99 | 3,394 / 2,835 |
| 203 | SS_Z6_R20_T150_C2 | gas1_ics1 | 9 | -93 | 6,018 [5,509, 6,490] | 155.5 | 0.93 | no | 1.00 | 2,543 / 3,061 |
| 204 | SO_Z6_R20_T100000_C2 | gas1_ics1 | 9 | -92 | 5,906 [5,400, 6,370] | 153.6 | 0.90 | no | 0.99 | 2,602 / 2,686 |
| 205 | SL_Z6_R30_T150_C5 | gas1_ics1 | 9 | -92 | 6,243 [5,774, 6,767] | 163.6 | 0.94 | no | 0.99 | 3,110 / 2,786 |
| 206 | SY_Z8_R20_T100000_C2 | gas1_ics1 | 9 | -90 | 6,468 [6,000, 6,988] | 178.9 | 0.93 | no | 1.00 | 3,080 / 2,962 |
| 207 | SS_Z10_R20_T150_C5 | gas1_ics1 | 9 | -84 | 6,131 [5,630, 6,604] | 163.2 | 0.94 | no | 0.99 | 2,827 / 2,944 |
| 208 | SO_Z8_R20_T100000_C2 | gas1_ics1 | 9 | -83 | 5,962 [5,457, 6,425] | 156.9 | 0.90 | no | 1.00 | 2,692 / 2,700 |
| 209 | SL_Z6_R20_T150_C5 | gas1_ics1 | 9 | -83 | 6,243 [5,775, 6,767] | 163.9 | 0.95 | no | 0.99 | 3,125 / 2,786 |
| 210 | SO_Z4_R20_T100000_C2 | gas1_ics1 | 9 | -81 | 5,793 [5,306, 6,281] | 147.2 | 0.89 | no | 0.99 | 2,498 / 2,648 |
| 211 | SO_Z8_R30_T150_C5 | gas1_ics1 | 9 | -81 | 5,906 [5,399, 6,373] | 151.8 | 0.92 | no | 0.99 | 2,692 / 2,756 |
| 212 | SO_Z4_R30_T100000_C2 | gas1_ics1 | 9 | -79 | 5,793 [5,310, 6,285] | 147.4 | 0.89 | no | 1.00 | 2,498 / 2,648 |
| 213 | SO_Z10_R30_T150_C5 | gas1_ics1 | 9 | -76 | 5,906 [5,419, 6,393] | 152.4 | 0.93 | no | 1.00 | 2,722 / 2,743 |
| 214 | SY_Z10_R30_T100000_C2 | gas1_ics1 | 9 | -76 | 6,580 [6,104, 7,095] | 182.6 | 0.94 | no | 1.00 | 3,215 / 2,978 |
| 215 | SS_Z10_R20_T100000_C2 | gas1_ics1 | 9 | -73 | 6,299 [5,818, 6,786] | 175.5 | 0.92 | no | 0.99 | 2,886 / 2,915 |
| 216 | SO_Z10_R20_T100000_C2 | gas1_ics1 | 9 | -68 | 6,018 [5,512, 6,478] | 160.6 | 0.91 | no | 0.99 | 2,767 / 2,701 |
| 217 | SS_Z8_R30_T150_C2 | gas1_ics1 | 9 | -68 | 6,074 [5,568, 6,549] | 158.0 | 0.94 | no | 0.99 | 2,678 / 3,022 |
| 218 | SL_Z10_R30_T100000_C2 | gas1_ics1 | 9 | -64 | 6,805 [6,285, 7,293] | 186.1 | 0.94 | no | 0.99 | 3,544 / 2,864 |
| 219 | SY_Z8_R30_T150_C5 | gas1_ics1 | 9 | -59 | 6,299 [5,803, 6,788] | 169.1 | 0.96 | no | 1.00 | 3,021 / 3,031 |
| 220 | SO_Z8_R20_T150_C5 | gas1_ics1 | 9 | -56 | 5,906 [5,412, 6,383] | 153.0 | 0.92 | no | 1.00 | 2,692 / 2,769 |
| 221 | SO_Z4_R20_T150_C5 | gas1_ics1 | 9 | -56 | 5,793 [5,292, 6,271] | 145.0 | 0.89 | no | 0.99 | 2,468 / 2,680 |
| 222 | SO_Z4_R30_T150_C5 | gas1_ics1 | 9 | -52 | 5,793 [5,295, 6,274] | 145.2 | 0.89 | no | 0.99 | 2,468 / 2,680 |
| 223 | SO_Z6_R30_T150_C5 | gas1_ics1 | 9 | -49 | 5,849 [5,340, 6,317] | 148.2 | 0.91 | no | 0.98 | 2,543 / 2,755 |
| 224 | SO_Z6_R20_T150_C5 | gas1_ics1 | 9 | -44 | 5,849 [5,346, 6,322] | 148.8 | 0.91 | no | 0.99 | 2,557 / 2,755 |
| 225 | SY_Z6_R30_T150_C2 | gas1_ics1 | 9 | -40 | 6,243 [5,734, 6,720] | 167.3 | 0.94 | no | 1.00 | 2,737 / 3,119 |
| 226 | SS_Z10_R30_T150_C2 | gas1_ics1 | 9 | -39 | 6,131 [5,633, 6,612] | 161.6 | 0.94 | no | 0.97 | 2,782 / 2,982 |
| 227 | SO_Z10_R20_T150_C5 | gas1_ics1 | 9 | -37 | 5,906 [5,444, 6,416] | 154.7 | 0.94 | no | 0.98 | 2,767 / 2,770 |
| 228 | SS_Z8_R20_T150_C2 | gas1_ics1 | 9 | -35 | 6,074 [5,601, 6,578] | 161.1 | 0.95 | no | 0.98 | 2,722 / 3,034 |
| 229 | SY_Z8_R20_T150_C5 | gas1_ics1 | 9 | -27 | 6,299 [5,812, 6,795] | 170.3 | 0.97 | no | 0.99 | 3,051 / 3,031 |
| 230 | SY_Z6_R20_T150_C2 | gas1_ics1 | 9 | -23 | 6,243 [5,736, 6,721] | 167.8 | 0.94 | no | 1.00 | 2,752 / 3,131 |
| 231 | SL_Z10_R20_T100000_C2 | gas1_ics1 | 9 | -19 | 6,805 [6,306, 7,309] | 188.8 | 0.95 | no | 0.99 | 3,589 / 2,877 |
| 232 | SY_Z10_R20_T100000_C2 | gas1_ics1 | 9 | -17 | 6,637 [6,133, 7,118] | 186.5 | 0.94 | no | 1.00 | 3,260 / 2,990 |
| 233 | SL_Z8_R30_T150_C5 | gas1_ics1 | 9 | -6 | 6,468 [5,965, 6,963] | 171.4 | 0.96 | no | 0.99 | 3,379 / 2,814 |
| 234 | SY_Z10_R30_T150_C5 | gas1_ics1 | 9 | -1 | 6,356 [5,880, 6,869] | 171.5 | 0.98 | no | 0.99 | 3,170 / 3,032 |
| 235 | SL_Z6_R30_T150_C2 | gas1_ics1 | 9 | +0 | 6,412 [5,904, 6,898] | 170.2 | 0.93 | no | 0.99 | 3,051 / 2,941 |
| 236 | SL_Z6_R20_T150_C2 | gas1_ics1 | 9 | +11 | 6,412 [5,904, 6,897] | 170.4 | 0.93 | no | 0.99 | 3,051 / 2,941 |
| 237 | SS_Z10_R20_T150_C2 | gas1_ics1 | 9 | +18 | 6,187 [5,674, 6,650] | 165.5 | 0.95 | no | 1.00 | 2,842 / 3,008 |
| 238 | SL_Z8_R20_T150_C5 | gas1_ics1 | 9 | +18 | 6,468 [5,970, 6,966] | 172.1 | 0.96 | no | 0.99 | 3,379 / 2,827 |
| 239 | SY_Z8_R30_T150_C2 | gas1_ics1 | 9 | +26 | 6,356 [5,861, 6,850] | 171.7 | 0.96 | no | 0.99 | 3,006 / 3,121 |
| 240 | SO_Z8_R30_T150_C2 | gas1_ics1 | 9 | +33 | 5,906 [5,423, 6,397] | 152.8 | 0.94 | no | 0.99 | 2,707 / 2,821 |
| 241 | SO_Z10_R30_T150_C2 | gas1_ics1 | 9 | +38 | 5,962 [5,447, 6,423] | 153.7 | 0.93 | no | 0.98 | 2,752 / 2,795 |
| 242 | SY_Z10_R20_T150_C5 | gas1_ics1 | 9 | +46 | 6,412 [5,896, 6,882] | 173.7 | 0.97 | no | 0.99 | 3,200 / 3,044 |
| 243 | SL_Z10_R30_T150_C5 | gas1_ics1 | 9 | +49 | 6,580 [6,050, 7,057] | 173.3 | 0.95 | no | 0.99 | 3,439 / 2,828 |
| 244 | SY_Z8_R20_T150_C2 | gas1_ics1 | 9 | +59 | 6,412 [5,895, 6,880] | 175.1 | 0.96 | no | 0.99 | 3,051 / 3,133 |
| 245 | SO_Z8_R20_T150_C2 | gas1_ics1 | 9 | +60 | 5,906 [5,437, 6,409] | 154.1 | 0.94 | no | 0.99 | 2,707 / 2,834 |
| 246 | SO_Z6_R30_T150_C2 | gas1_ics1 | 9 | +71 | 5,849 [5,359, 6,337] | 149.2 | 0.92 | no | 0.99 | 2,558 / 2,820 |
| 247 | SY_Z10_R30_T150_C2 | gas1_ics1 | 9 | +72 | 6,468 [5,958, 6,950] | 175.3 | 0.98 | no | 0.99 | 3,186 / 3,122 |
| 248 | SO_Z4_R20_T150_C2 | gas1_ics1 | 9 | +79 | 5,793 [5,329, 6,310] | 146.9 | 0.90 | no | 1.00 | 2,438 / 2,783 |
| 249 | SL_Z8_R30_T150_C2 | gas1_ics1 | 9 | +79 | 6,580 [6,074, 7,074] | 177.2 | 0.95 | no | 0.99 | 3,246 / 2,982 |
| 250 | SO_Z6_R20_T150_C2 | gas1_ics1 | 9 | +81 | 5,849 [5,363, 6,339] | 149.6 | 0.92 | no | 0.99 | 2,558 / 2,833 |
| 251 | SO_Z10_R20_T150_C2 | gas1_ics1 | 9 | +81 | 5,962 [5,478, 6,449] | 156.3 | 0.94 | no | 0.99 | 2,782 / 2,821 |
| 252 | SO_Z4_R30_T150_C2 | gas1_ics1 | 9 | +83 | 5,849 [5,333, 6,314] | 147.1 | 0.89 | no | 1.00 | 2,438 / 2,783 |
| 253 | SL_Z10_R20_T150_C5 | gas1_ics1 | 9 | +85 | 6,580 [6,058, 7,062] | 174.7 | 0.96 | no | 0.99 | 3,469 / 2,841 |
| 254 | SL_Z8_R20_T150_C2 | gas1_ics1 | 9 | +103 | 6,580 [6,100, 7,098] | 179.3 | 0.95 | no | 0.99 | 3,291 / 2,982 |
| 255 | SY_Z10_R20_T150_C2 | gas1_ics1 | 9 | +123 | 6,468 [5,986, 6,973] | 178.2 | 0.98 | no | 0.99 | 3,216 / 3,134 |
| 256 | SL_Z10_R30_T150_C2 | gas1_ics1 | 9 | +134 | 6,693 [6,179, 7,187] | 181.1 | 0.95 | no | 0.99 | 3,351 / 2,996 |
| 257 | SL_Z10_R20_T150_C2 | gas1_ics1 | 9 | +174 | 6,693 [6,197, 7,202] | 182.8 | 0.95 | no | 0.99 | 3,381 / 3,008 |
| 258 | gll_iem_v07 | - | 6 | +749 | 3,150 [2,614, 3,663] | 36.4 | 0.66 | yes | 0.99 | 592 / 1,473 |
