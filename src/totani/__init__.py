"""Independent check of the Totani (arXiv:2507.07209) 20 GeV halo excess.

Importing this package puts the VENDORED gce-fisher code (`totani/_vendor/`)
on `sys.path` -- see `totani._paths` for why the code is copied and pinned
rather than imported from the sibling checkout, and what is still read from
that checkout (data only).
"""
from __future__ import annotations

from ._paths import (GCE_FISHER_ROOT, VENDOR_ROOT, ensure_on_path,
                     ensure_scripts_on_path, gce_fisher_data)

ensure_on_path()

__all__ = ["GCE_FISHER_ROOT", "VENDOR_ROOT", "ensure_on_path",
           "ensure_scripts_on_path", "gce_fisher_data"]
