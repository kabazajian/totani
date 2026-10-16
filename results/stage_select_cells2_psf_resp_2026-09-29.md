# Structure selection and scoring — cells2_psf_resp (2026-09-29)

Bins used for selection AND inference: [3, 4, 5, 6, 7, 8] (K = 6). Fraction-of-rows rule, cutoff 0.8; BIC p ln n with n = cells × K. Voters: 16 stratified 2012 models + A2020 per decomposition; Galp21 rows scored, not voting. 'favor' = TS > 4 at 20.8 GeV, positive, all four quadrants positive. Q⁺ over the K bins with the binomial χ² mixture; Z₂₀₁₈ = Φ⁻¹(1 − p/2).

## gas1_ics1, residual

Selected: **bubble_totani → loopI_shells+nps → erosita_bright_NS**

- 2012 subset: favor 17/17; halo20 median 6,642 [6,298–8,708]; Z₂₀₁₈ median 27.3 [19.6–28.1]; ⟨σv⟩ median +6.6, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 7,120 [6,573–7,545]; Z₂₀₁₈ median 31.8 [27.1–33.4]; ⟨σv⟩ median +7.3, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}, 'G21': {'bubble_totani': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median -185 (flat larger in 2/17 rows)

## gas1_ics1, no residual

Selected: **bubble_structured → erosita_unif_NS → loopI_haslam_raw**

- 2012 subset: favor 17/17; halo20 median 4,087 [3,175–5,796]; Z₂₀₁₈ median 12.6 [7.5–17.4]; ⟨σv⟩ median +3.0, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 4,045 [3,629–4,411]; Z₂₀₁₈ median 12.0 [9.1–14.9]; ⟨σv⟩ median +3.0, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_structured': 15, 'bubble_flat': 2}, 'G21': {'bubble_structured': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median -185 (flat larger in 2/17 rows)

## gas5_ics3, residual

Selected: **bubble_totani → loopI_shells+nps**

- 2012 subset: favor 17/17; halo20 median 7,464 [7,091–8,708]; Z₂₀₁₈ median 15.5 [10.9–21.8]; ⟨σv⟩ median +7.6, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 7,521 [6,999–8,023]; Z₂₀₁₈ median 15.0 [13.3–16.1]; ⟨σv⟩ median +7.7, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_totani': 17}, 'G21': {'bubble_totani': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0, 'bubble_totani': 1.0}; flat − structured improvement per row median +234 (flat larger in 16/17 rows)

## gas5_ics3, no residual

Selected: **bubble_flat → erosita_outer_NS**

- 2012 subset: favor 17/17; halo20 median 6,175 [5,291–7,594]; Z₂₀₁₈ median 11.2 [6.2–21.3]; ⟨σv⟩ median +6.1, 2σ+ 17/17
- Galp21: favor 12/12; halo20 median 6,757 [6,169–7,290]; Z₂₀₁₈ median 15.2 [12.4–18.1]; ⟨σv⟩ median +6.9, 2σ+ 12/12
- own paths, bubble chosen: {'2012': {'bubble_flat': 16, 'bubble_structured': 1}, 'G21': {'bubble_flat': 12}}

- step 1: clears {'bubble_flat': 1.0, 'bubble_structured': 1.0}; flat − structured improvement per row median +234 (flat larger in 16/17 rows)

