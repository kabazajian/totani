#!/usr/bin/env python
"""Totani's five GALPROP parameter variants (his sec 3.4.4, Figure 15 bottom).

    python scripts/totani_galprop_variants.py     # ~6 min once the archives are down

Each is one step away from his baseline SLZ6R30T150C2 on the Ackermann+ 2012
grid: SNR source distribution (SS_Z6_R30_T150_C2), z_h = 10 kpc
(SL_Z10_R30_T150_C2), R_h = 20 kpc (SL_Z6_R20_T150_C2), T_S = 1e5 K
(SL_Z6_R30_T100000_C2), E(B-V) cut 5 mag (SL_Z6_R30_T150_C5). The archives
are downloaded from galprop.stanford.edu/PaperIISuppMaterial into
`external/galprop_paperII/` (gitignored, ~1.1-1.4 GB each) and only the three
total skymaps each variant needs (pi0_decay, bremss, ics_isotropic) are
extracted.

Each variant is fitted as Totani's baseline configuration (gas + ics, two
free templates per bin) with BOTH bubble templates, and scored the standard
way. He reports (sec 3.4.4) that most variants leave the result almost
unchanged, but SNR sources make the halo positive in every bin and near zero
in the lowest -- which is a specific prediction to check.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from reduce_totani import TOTANI_CENTERS_GEV, TOTANI_EDGES_GEV  # noqa: E402
from systematics import bubble_cell_template, replace_component  # noqa: E402
from totani.galprop import TOTANI_SET  # noqa: E402
from totani.stack import build_background_stack, halo_cell_template  # noqa: E402
from totani_diffuse import score  # noqa: E402

RUN = ROOT / "runs" / "fff2de5bc7e0"
EDGES = TOTANI_EDGES_GEV[:10]
ARCH = ROOT / "external" / "galprop_paperII"

#: archive stem -> (galdef stem, what changed from the baseline)
VARIANTS = {
    "SS_Z6_R30_T150_C2": ("54_SNR_z6kpc_R30kpc_Ts150K_EBV2mag", "SNR source distribution"),
    "SL_Z10_R30_T150_C2": ("54_Lorimer_z10kpc_R30kpc_Ts150K_EBV2mag", "z_h = 10 kpc"),
    "SL_Z6_R20_T150_C2": ("54_Lorimer_z6kpc_R20kpc_Ts150K_EBV2mag", "R_h = 20 kpc"),
    "SL_Z6_R30_T100000_C2": ("54_Lorimer_z6kpc_R30kpc_Ts100000K_EBV2mag", "T_S = 1e5 K"),
    "SL_Z6_R30_T150_C5": ("54_Lorimer_z6kpc_R30kpc_Ts150K_EBV5mag", "E(B-V) cut 5 mag"),
}
NEEDED = ("pi0_decay", "bremss", "ics_isotropic")


def ensure_extracted(archive: str, model_set: str) -> None:
    tar = ARCH / f"{archive}.tar"
    if not tar.exists():
        raise FileNotFoundError(f"{tar} not downloaded")
    members = [f"OUTPUT/{c}_healpix_{model_set}.gz" for c in NEEDED]
    missing = [m for m in members if not (ARCH / m).exists()]
    if missing:
        subprocess.run(["tar", "xf", str(tar), *missing], cwd=ARCH, check=True)


def main() -> int:
    import resource
    t0 = time.time()
    J = halo_cell_template(RUN, EDGES)
    flat = bubble_cell_template(RUN, EDGES, "ackermann")
    results = {"date": "2026-09-18", "run": "fff2de5bc7e0", "baseline_set": TOTANI_SET,
               "E_centers_gev": list(TOTANI_CENTERS_GEV[:9]), "variants": {}}
    sets = {"baseline_SLZ6R30T150C2": (TOTANI_SET, "Totani baseline")}
    sets.update({k: v for k, v in VARIANTS.items()})
    for tag, (model_set, what) in sets.items():
        if tag != "baseline_SLZ6R30T150C2":
            ensure_extracted(tag, model_set)
        for bubble in ("macias", "ackermann"):
            # `model_set` is threaded through galprop_on_target via the stack's
            # `galprop_model_set`; see build_background_stack.
            st = build_background_stack(RUN, EDGES, bubble="macias", diffuse="galprop_totani",
                                        galprop_model_set=model_set)
            if bubble == "ackermann":
                st = replace_component(st, "bubble", flat)
            name = f"{tag}_{'structured' if bubble == 'macias' else 'flat'}"
            d = score(st, J); d["what"] = what; d["model_set"] = model_set
            results["variants"][name] = d
            h5 = d["halo_intervals"][5]; q5 = d["summary"]["bin5"]
            print(f"  {name:38s} halo@20.8 {h5['counts']:6.0f} [{h5['lo']:.0f}, {h5['hi']:.0f}] "
                  f"TS {d['ts'][5]:5.1f} quad {q5['sum_over_baseline']:+.2f} neg={q5['any_negative']} "
                  f"({time.time()-t0:.0f}s)", flush=True)
            print("      halo/bin   " + " ".join(f"{h['counts']:8.0f}" for h in d["halo_intervals"]))
    out = ROOT / "results" / "fits" / "totani_galprop_variants_2026-09-18.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"\nsaved {out}\ntotal {time.time()-t0:.0f}s; peak RSS: %.2f GB"
          % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
