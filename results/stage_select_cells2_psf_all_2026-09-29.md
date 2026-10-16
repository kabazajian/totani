# Structure selection and scoring — cells2_psf_all (2026-09-29)

Bins used for selection AND inference: [0, 1, 2, 3, 4, 5, 6, 7, 8] (K = 9). Fraction-of-rows rule, cutoff 0.8; BIC p ln n with n = cells × K. Voters: 16 stratified 2012 models + A2020 per decomposition; Galp21 rows scored, not voting. 'favor' = TS > 4 at 20.8 GeV, positive, all four quadrants positive. Q⁺ over the K bins with the binomial χ² mixture; Z₂₀₁₈ = Φ⁻¹(1 − p/2).

## gas1_ics1, residual

Selected: **bubble_totani → loopI_shells+nps → erosita_bright_NS**

- 2012 subset: favor 17/17; halo20 median 6,642 [6,298–8,708]; Z₂₀₁₈ median 27.7 [22.2–34.5]; ⟨σv⟩ median +5.3, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 7,120 [6,573–7,545]; Z₂₀₁₈ median 36.0 [27.8–39.4]; ⟨σv⟩ median +6.9, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}, 'G21': {'bubble_totani': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median +1103 (flat larger in 15/17 rows)

## gas1_ics1, no residual

Selected: **bubble_flat → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic**

- 2012 subset: favor 17/17; halo20 median 5,745 [4,777–7,408]; Z₂₀₁₈ median 28.7 [16.2–44.5]; ⟨σv⟩ median +5.5, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 5,698 [5,248–6,109]; Z₂₀₁₈ median 28.7 [20.0–37.5]; ⟨σv⟩ median +5.7, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_flat': 15, 'bubble_structured': 2}, 'G21': {'bubble_flat': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median +1103 (flat larger in 15/17 rows)

## gas5_ics3, residual

Selected: **bubble_totani → loopI_shells+nps → erosita_unif_NS**

- 2012 subset: favor 17/17; halo20 median 7,463 [7,020–8,708]; Z₂₀₁₈ median 13.9 [9.6–23.7]; ⟨σv⟩ median +7.3, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 7,461 [6,999–7,918]; Z₂₀₁₈ median 14.7 [13.6–17.6]; ⟨σv⟩ median +7.4, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}, 'G21': {'bubble_totani': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median +1883 (flat larger in 17/17 rows)

## gas5_ics3, no residual

Selected: **bubble_flat → erosita_unif_NS → loopI_haslam_raw → sun_moon_ic**

- 2012 subset: favor 9/17; halo20 median 6,226 [5,411–7,408]; Z₂₀₁₈ median 11.5 [5.2–33.6]; ⟨σv⟩ median +6.0, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 6,719 [6,109–6,891]; Z₂₀₁₈ median 18.9 [13.9–23.2]; ⟨σv⟩ median +7.3, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_flat': 17}, 'G21': {'bubble_flat': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median +1883 (flat larger in 17/17 rows)

