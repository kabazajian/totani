# Vendored from gce-fisher

Copied **byte for byte** on 2026-09-18 from `~/work/gce-fisher` at commit
`62f05c79d841e14a881f8c51455753d470f29d2a` (2026-09-16). None of these files
was modified in that working tree at the time of the copy.

| here | upstream |
|---|---|
| `src/gcereduce/__init__.py` | `src/gcereduce/__init__.py` |
| `src/gcereduce/reduce.py` | `src/gcereduce/reduce.py` |
| `src/gcereduce/regrid.py` | `src/gcereduce/regrid.py` |
| `src/gcereduce/dm_templates.py` | `src/gcereduce/dm_templates.py` |
| `src/gcereduce/forward.py` | `src/gcereduce/forward.py` (lazy import of `dm_templates`) |
| `src/gcereduce/galp21.py` | `src/gcereduce/galp21.py` |
| `scripts/a2020_templates_official.py` | `scripts/a2020_templates_official.py` |
| `scripts/a2020_stage2_diffuse.py` | `scripts/a2020_stage2_diffuse.py` |
| `scripts/a2020_stage2_reduce.py` | `scripts/a2020_stage2_reduce.py` (lazy import of the other two) |

The layout mirrors upstream so the intra-imports (`from gcereduce.regrid
import ...`, `from a2020_stage2_reduce import ...`) work unchanged and
`diff -r` against a checkout is meaningful.

## Pinned

`MANIFEST.sha256` holds the sha256 of every `.py` here and
`tests/test_imports.py::test_vendored_files_match_manifest` asserts it. To
change a vendored file: edit it, regenerate the manifest
(`cd src/totani/_vendor && find . -name '*.py' | sort | xargs shasum -a 256 >
MANIFEST.sha256`), and record what changed and why below.

## Known harmless artefacts of the byte-identical copy

The two A2020 scripts compute `ROOT = Path(__file__).resolve().parents[1]`
and derive default data paths (`SRC`, `SRC_NORMALIZED`, `POHL`) from it, and
insert `ROOT/src` and `ROOT/scripts` on `sys.path`. Here `ROOT` is this
`_vendor/` directory, so those defaults point at directories that do not
exist. Nothing in `totani` uses them: `templates.delivered_paths` builds the
delivery paths from `GCE_FISHER_ROOT` and `ics_on_target` passes `src_dir`
explicitly. `ROOT/src` on `sys.path` is exactly the vendored `gcereduce`,
which is what we want.

## Local modifications

None.

## Why not a git submodule

gce-fisher is a private checkout without a remote this repo can pin to, and a
submodule would still pull the whole 2,786-line CONVENTIONS and every script.
Nine files is the actual dependency.
