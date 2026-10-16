# Structure selection and scoring — cells1_psf_all (2026-09-29)

Bins used for selection AND inference: [0, 1, 2, 3, 4, 5, 6, 7, 8] (K = 9). Fraction-of-rows rule, cutoff 0.8; BIC p ln n with n = cells × K. Voters: 16 stratified 2012 models + A2020 per decomposition; Galp21 rows scored, not voting. 'favor' = TS > 4 at 20.8 GeV, positive, all four quadrants positive. Q⁺ over the K bins with the binomial χ² mixture; Z₂₀₁₈ = Φ⁻¹(1 − p/2).

## gas1_ics1, residual

Selected: **bubble_totani → loopI_shells+nps → erosita_bright_NS**

- 2012 subset: favor 17/17; halo20 median 6,900 [6,495–8,896]; Z₂₀₁₈ median 29.7 [24.8–37.1]; ⟨σv⟩ median +5.8, 2σ+ 17/17
- Galp21: favor 1/1; halo20 median 6,855 [6,855–6,855]; Z₂₀₁₈ median 32.1 [32.1–32.1]; ⟨σv⟩ median +6.1, 2σ+ 1/1
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}, 'G21': {'bubble_totani': 1}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median +1290 (flat larger in 16/17 rows)

## gas1_ics1, no residual

Selected: **bubble_flat → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic**

- 2012 subset: favor 17/17; halo20 median 5,802 [4,847–7,504]; Z₂₀₁₈ median 29.9 [17.2–45.6]; ⟨σv⟩ median +5.7, 2σ+ 17/17
- Galp21: favor 1/1; halo20 median 5,779 [5,779–5,779]; Z₂₀₁₈ median 28.2 [28.2–28.2]; ⟨σv⟩ median +5.5, 2σ+ 1/1
- own paths, bubble chosen: {'2012': {'bubble_flat': 16, 'bubble_structured': 1}, 'G21': {'bubble_flat': 1}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median +1290 (flat larger in 16/17 rows)

## gas5_ics3, residual

Selected: **bubble_totani → loopI_shells+nps → erosita_unif_NS**

- 2012 subset: favor 17/17; halo20 median 7,423 [7,071–8,896]; Z₂₀₁₈ median 15.7 [11.1–26.5]; ⟨σv⟩ median +7.6, 2σ+ 17/17
- Galp21: —
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median +2081 (flat larger in 17/17 rows)

## gas5_ics3, no residual

Selected: **bubble_flat → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic**

- 2012 subset: favor 8/17; halo20 median 6,301 [5,602–7,504]; Z₂₀₁₈ median 12.9 [6.2–34.0]; ⟨σv⟩ median +6.4, 2σ+ 17/17
- Galp21: —
- own paths, bubble chosen: {'2012': {'bubble_flat': 17}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median +2081 (flat larger in 17/17 rows)

