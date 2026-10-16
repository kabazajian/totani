# Structure selection and scoring — cells1_psfuc_all (2026-09-29)

Bins used for selection AND inference: [0, 1, 2, 3, 4, 5, 6, 7, 8] (K = 9). Fraction-of-rows rule, cutoff 0.8; BIC p ln n with n = cells × K. Voters: 16 stratified 2012 models + A2020 per decomposition; Galp21 rows scored, not voting. 'favor' = TS > 4 at 20.8 GeV, positive, all four quadrants positive. Q⁺ over the K bins with the binomial χ² mixture; Z₂₀₁₈ = Φ⁻¹(1 − p/2).

## gas1_ics1, residual

Selected: **bubble_totani → erosita_unif_NS → loopI_shells+nps**

- 2012 subset: favor 17/17; halo20 median 6,829 [6,483–8,910]; Z₂₀₁₈ median 32.4 [27.7–44.3]; ⟨σv⟩ median +6.4, 2σ+ 17/17
- Galp21: favor 1/1; halo20 median 7,212 [7,212–7,212]; Z₂₀₁₈ median 35.0 [35.0–35.0]; ⟨σv⟩ median +6.3, 2σ+ 1/1
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}, 'G21': {'bubble_totani': 1}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median +1280 (flat larger in 16/17 rows)

## gas1_ics1, no residual

Selected: **bubble_flat → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic**

- 2012 subset: favor 17/17; halo20 median 5,821 [4,876–7,521]; Z₂₀₁₈ median 30.0 [17.2–45.7]; ⟨σv⟩ median +5.7, 2σ+ 17/17
- Galp21: favor 1/1; halo20 median 5,800 [5,800–5,800]; Z₂₀₁₈ median 28.3 [28.3–28.3]; ⟨σv⟩ median +5.6, 2σ+ 1/1
- own paths, bubble chosen: {'2012': {'bubble_flat': 16, 'bubble_structured': 1}, 'G21': {'bubble_flat': 1}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median +1280 (flat larger in 16/17 rows)

## gas5_ics3, residual

Selected: **bubble_totani → loopI_shells+nps → erosita_unif_NS**

- 2012 subset: favor 17/17; halo20 median 7,447 [7,088–8,910]; Z₂₀₁₈ median 15.9 [11.2–27.7]; ⟨σv⟩ median +7.7, 2σ+ 17/17
- Galp21: —
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median +2080 (flat larger in 17/17 rows)

## gas5_ics3, no residual

Selected: **bubble_flat → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic**

- 2012 subset: favor 10/17; halo20 median 6,335 [5,651–7,521]; Z₂₀₁₈ median 12.9 [6.2–34.0]; ⟨σv⟩ median +6.4, 2σ+ 17/17
- Galp21: —
- own paths, bubble chosen: {'2012': {'bubble_flat': 17}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median +2080 (flat larger in 17/17 rows)

