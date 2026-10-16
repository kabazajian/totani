"""Portable manuscript rendering from validated numeric figure bundles.

No analysis loaders, likelihood code, Astropy or network access is needed
here. Text is typeset by LaTeX in Computer Modern to match the manuscript, so
a TeX installation (latex, dvipng) must be on PATH. Dense map layers are rasterized; annotations,
mask edges, markers and error bars remain vectors in the PDF.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


TEXT_WIDTH_INCHES = 510.0 / 72.27
SPATIAL_HEIGHT = 6.6
STYLE = {
    # LaTeX with Computer Modern Roman, matching the REVTeX manuscript text.
    "font.family": "serif", "font.serif": ["Computer Modern Roman"],
    "text.latex.preamble": r"\usepackage{amsmath}", "font.size": 10,
    "axes.labelsize": 10, "axes.titlesize": 10,
    "xtick.labelsize": 10, "ytick.labelsize": 10,
    "legend.fontsize": 10, "axes.linewidth": 0.8,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "pdf.fonttype": 42,
    "ps.fonttype": 42, "text.usetex": True,
}
EXCLUDED = "#dedede"
MISSING = "#dcb8a3"
BOUNDARY = "#242424"
COMMON = "#0072b2"
OUTSIDE = "#d55e00"


def _check_provenance(metadata, test_data):
    """Direct callers cannot turn a fixture into an unmarked production PDF."""
    provenance = metadata.get("provenance", "scientific")
    if provenance not in ("scientific", "test"):
        raise ValueError("Unknown scientific provenance")
    if provenance == "test" and not test_data:
        raise ValueError("Production mode rejects TEST DATA")


def _save(fig, output_dir, basename, test_data):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if test_data:
        fig.text(0.5, 0.49, "TEST DATA - NOT SCIENTIFIC RESULTS", ha="center",
                 va="center", rotation=24, fontsize=20, color="#9b1616",
                 alpha=0.32, zorder=100)
        fig.text(0.5, 0.018, "TEST DATA - NOT SCIENTIFIC RESULTS", ha="center",
                 fontsize=10, color="#9b1616", zorder=100)
    pdf, png = (output_dir / f"{basename}.{suffix}" for suffix in ("pdf", "png"))
    # Preserve the exact placed width: tight bounding boxes change page size.
    fig.savefig(pdf, dpi=600)
    fig.savefig(png, dpi=300)
    rgba = plt.imread(png)
    luminance = rgba[..., :3] @ np.array([0.2126, 0.7152, 0.0722])
    gray = output_dir / f"{basename}_grayscale.png"
    plt.imsave(gray, luminance, cmap="gray", vmin=0, vmax=1, dpi=300)
    plt.close(fig)
    return [pdf, png, gray]


def _boundary(ax, mask, longitude_edges, latitude_edges, *, color=BOUNDARY,
              linestyle="solid", linewidth=0.8, zorder=6):
    """Draw exact pixel edges of a Boolean mask, never interpolated contours."""
    mask = np.asarray(mask, dtype=bool)
    xe, ye = np.asarray(longitude_edges), np.asarray(latitude_edges)
    padded = np.pad(mask, 1)
    # Boundaries parallel to latitude: one edge for every horizontal flip.
    vertical = padded[1:-1, 1:] != padded[1:-1, :-1]
    segments = []
    for ix in np.flatnonzero(vertical.any(axis=0)):
        indices = np.flatnonzero(vertical[:, ix])
        breaks = np.flatnonzero(np.diff(indices) != 1)
        starts, ends = indices[np.r_[0, breaks + 1]], indices[np.r_[breaks, len(indices) - 1]]
        segments.extend([[(xe[ix], ye[start]), (xe[ix], ye[end + 1])]
                         for start, end in zip(starts, ends)])
    horizontal = padded[1:, 1:-1] != padded[:-1, 1:-1]
    for iy in np.flatnonzero(horizontal.any(axis=1)):
        indices = np.flatnonzero(horizontal[iy])
        breaks = np.flatnonzero(np.diff(indices) != 1)
        starts, ends = indices[np.r_[0, breaks + 1]], indices[np.r_[breaks, len(indices) - 1]]
        segments.extend([[(xe[start], ye[iy]), (xe[end + 1], ye[iy])]
                         for start, end in zip(starts, ends)])
    collection = LineCollection(segments, colors=color, linewidths=linewidth,
                                linestyles=linestyle, zorder=zorder)
    ax.add_collection(collection)
    return collection


def _overlay(ax, mask, extent, color, *, alpha=1, zorder=4):
    cmap = colors.ListedColormap([color])
    layer = np.ma.masked_where(~np.asarray(mask, bool), np.zeros(np.shape(mask)))
    ax.imshow(layer, origin="lower", extent=extent, interpolation="nearest",
              cmap=cmap, vmin=0, vmax=1, alpha=alpha, rasterized=True, zorder=zorder)


def _sky_axes(ax, longitude_edges, latitude_edges, *, labels=True):
    ax.set_xlim(longitude_edges[0], longitude_edges[-1])
    ax.set_ylim(latitude_edges[0], latitude_edges[-1])
    ax.set_aspect("equal", adjustable="box")
    # Edges are inherited from the WCS; no flip or array-handedness guess.
    lo, hi = min(longitude_edges), max(longitude_edges)
    blo, bhi = min(latitude_edges), max(latitude_edges)
    ticks = np.array([-60, -30, 0, 30, 60])
    ax.set_xticks(ticks[(ticks >= lo) & (ticks <= hi)])
    ax.set_yticks(ticks[(ticks >= blo) & (ticks <= bhi)])
    if labels:
        ax.set_xlabel(r"Galactic longitude $\ell$ [deg]", labelpad=3)
        ax.set_ylabel(r"Galactic latitude $b$ [deg]", labelpad=3)


def render_spatial(metadata, arrays, output_dir, test_data=False,
                   width_inches=TEXT_WIDTH_INCHES):
    """Render the four actual input-intensity maps and their exact footprint."""
    _check_provenance(metadata, test_data)
    xe = np.asarray(arrays["longitude_edges"])
    ye = np.asarray(arrays["latitude_edges"])
    extent = [xe[0], xe[-1], ye[0], ye[-1]]
    observed = np.asarray(arrays["observed"], bool)
    titles = [("nfw", "(a) NFW annihilation template"),
              ("flat", r"(b) Flat \textit{Fermi} bubbles template"),
              ("structured", r"(c) Structured \textit{Fermi} bubbles template"),
              ("erosita", "(d) eROSITA bubbles template weights")]
    with plt.rc_context(STYLE):
        # One-line titles, latitude labels only on the left column and
        # longitude labels only on the bottom row give the width to the maps.
        fig = plt.figure(figsize=(width_inches, SPATIAL_HEIGHT))
        gs = fig.add_gridspec(2, 2, left=0.085, right=0.94, bottom=0.135,
                              top=0.955, wspace=0.2, hspace=0.17)
        any_missing = False
        for k, (key, title) in enumerate(titles):
            ax = fig.add_subplot(gs[k // 2, k % 2])
            descriptor = metadata["maps"][key]
            display = descriptor["display_normalization"]
            values = np.asarray(arrays[key], float) / float(display["factor"])
            finite = np.isfinite(values)
            visible = values[finite & observed]
            if not visible.size:
                raise ValueError(f"{key}: no finite observed pixels to display")
            vmin = min(float(np.min(visible)), 0.0)
            vmax = max(float(np.max(visible)), np.finfo(float).eps)
            cmap = plt.get_cmap("RdBu_r" if vmin < 0 else "cividis").copy()
            if vmin < 0:
                vmax = max(vmax, -vmin)
                vmin = -vmax
            if display["stretch"] == "asinh":
                norm = colors.AsinhNorm(linear_width=0.03 * vmax,
                                        vmin=vmin, vmax=vmax)
            else:
                norm = colors.Normalize(vmin=vmin, vmax=vmax)
            image = ax.imshow(np.ma.masked_where(~finite | ~observed, values),
                              origin="lower", extent=extent, cmap=cmap,
                              norm=norm, interpolation="nearest", rasterized=True)
            _overlay(ax, ~observed, extent, EXCLUDED)
            _overlay(ax, ~finite, extent, MISSING, zorder=5)
            any_missing |= bool(np.any(~finite))
            _boundary(ax, arrays["source_excluded"], xe, ye,
                      color="#777777", linewidth=0.35)
            _boundary(ax, arrays["flat_support"], xe, ye)
            if key == "nfw":
                _boundary(ax, arrays["outside_common"], xe, ye,
                          color=OUTSIDE, linestyle="dashed", linewidth=0.85)
            if key == "erosita":
                _boundary(ax, arrays["erosita_support"], xe, ye,
                          color=BOUNDARY, linewidth=0.8)
                # Provenance masks come from metadata, never screen position.
                native = np.asarray(arrays["erosita_native"], bool)
                reconstructed = np.asarray(arrays["erosita_reconstructed"], bool)
                _boundary(ax, native, xe, ye, color="#555555",
                          linestyle="dashed", linewidth=0.65)
            ax.set_title(title, loc="left", pad=6)
            _sky_axes(ax, xe, ye)
            if key == "erosita":
                centers_x = (xe[:-1] + xe[1:]) / 2
                x_left, x_right = ax.get_xlim()
                for mask, name in [(native, "Native"), (reconstructed, "Reconstructed")]:
                    row, col = np.nonzero(mask & observed)
                    if row.size:
                        # Each name sits in the top corner on its own support's
                        # longitude side, inside the axes and clear of tick labels.
                        on_left = (np.median(centers_x[col]) > 0) == (x_left > x_right)
                        ax.text(0.02 if on_left else 0.98, 0.975, name,
                                transform=ax.transAxes, ha="left" if on_left else "right",
                                va="top", fontsize=10, color="black",
                                bbox=dict(facecolor="white", alpha=0.82,
                                          linewidth=0, pad=1.2), zorder=9)
            if k % 2:
                ax.set_ylabel("")
                ax.tick_params(labelleft=False)
            if k < 2:
                ax.set_xlabel("")
            cb = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.028, shrink=0.92)
            cb.ax.tick_params(labelsize=10)
            if display["stretch"] == "asinh":
                ticks = [-vmax, -.2 * vmax, 0, .2 * vmax, vmax] if vmin < 0 else [0, .03 * vmax, .2 * vmax, vmax]
            else:
                ticks = [-vmax, -.5 * vmax, 0, .5 * vmax, vmax] if vmin < 0 else [0, .5 * vmax, vmax]
            # Small maxima get a multiplier so tick labels stay narrow enough
            # not to run into the neighboring panel's latitude labels.
            power = int(np.floor(np.log10(vmax))) if vmax < 0.1 else 0
            quantity = r"w/w_0" if descriptor["quantity"] == "template_weight" else r"I/I_0"
            cb.ax.set_title(rf"${quantity}$" if not power else rf"${quantity}$" + "\n" + rf"$[10^{{{power}}}]$",
                            fontsize=10, pad=6, loc="left")
            cb.set_ticks(ticks, labels=[f"{value / 10.0**power:.2g}" for value in ticks])
            cb.ax.yaxis.get_offset_text().set_fontsize(10)
        handles = [Line2D([], [], color=BOUNDARY, lw=0.8, label=r"Flat \textit{Fermi} bubbles boundary"),
                   Line2D([], [], color=OUTSIDE, lw=0.85, ls="--", label="Common outside region"),
                   Patch(facecolor=EXCLUDED, label="Excluded")]
        if any_missing:
            handles.append(Patch(facecolor=MISSING, label="Missing"))
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
                   ncol=len(handles), frameon=False, handlelength=2.4, columnspacing=1.8)
        return _save(fig, output_dir, "template_overview", test_data)


def _schematic(ax, schematic, partition, labels):
    xe = np.asarray(schematic["longitude_edges"])
    ye = np.asarray(schematic["latitude_edges"])
    regions = np.asarray(schematic[partition], dtype=int)
    observed = np.asarray(schematic["observed"], bool)
    extent = [xe[0], xe[-1], ye[0], ye[-1]]
    cmap = colors.ListedColormap(["#f2f2f2", "#dbdbdb", "#bcbcbc", "#969696"])
    ax.imshow(np.ma.masked_where(~observed, regions), extent=extent, origin="lower",
              interpolation="nearest", cmap=cmap, vmin=-0.5, vmax=3.5,
              rasterized=True)
    _overlay(ax, ~observed, extent, "white")
    cx, cy = (xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2
    for index, label in enumerate(labels):
        mask = (regions == index) & observed
        _boundary(ax, mask, xe, ye, linewidth=0.5)
        row, col = np.nonzero(mask)
        if row.size:
            target_x, target_y = np.median(cx[col]), np.median(cy[row])
            if partition == "rot" and label in ("E", "W"):
                # The E/W wedges straddle the excluded plane: their centroid
                # is not itself in the observed region. Select an actual
                # interior label pixel in the northern part of that mask.
                northern = cy[row] > 0
                target_x, target_y = np.median(cx[col[northern]]), np.median(cy[row[northern]])
            closest = np.argmin((cx[col] - target_x)**2 + (cy[row] - target_y)**2)
            ax.text(cx[col[closest]], cy[row[closest]], label, ha="center", va="center",
                    fontsize=10, bbox=dict(facecolor="white", alpha=0.88, linewidth=0, pad=1))
    _sky_axes(ax, xe, ye, labels=False)
    ax.set_xticks([-60, 0, 60])
    ax.set_yticks([-60, 0, 60])
    ax.set_xlabel(r"$\ell$ [deg]; E: $\ell>0$, W: $\ell<0$", fontsize=10, labelpad=2)
    ax.set_ylabel(r"$b$ [deg]", fontsize=10, labelpad=1)


def _number(value):
    return f"{float(value):.1f}"


def render_regional(data, output_dir, test_data=False,
                    width_inches=TEXT_WIDTH_INCHES):
    """Render signed common-basis amplitudes and saved profile intervals."""
    _check_provenance(data, test_data)
    illustrated = int(data["illustrated_bin"])
    panels = {(p["config"], p["partition"]): p for p in data["panels"]}
    configurations = ["structured", "flat", "flat_flexlobes"]
    row_titles = [r"Structured \textit{Fermi} bubbles; uniform eROSITA N/S",
                  r"Flat \textit{Fermi} bubbles; uniform eROSITA N/S",
                  r"Flat \textit{Fermi} bubbles; twelve X-ray-brightness-weighted pieces"]
    bounds = []
    for panel in panels.values():
        for name in ["coefficients", "profile_lo", "profile_hi"]:
            bounds.extend(np.asarray(panel[name][illustrated], float).ravel())
        bounds.append(float(panel["common"][illustrated]))
    bounds = np.asarray(bounds)
    if not np.all(np.isfinite(bounds)):
        raise ValueError("regional coefficients and illustrated profile bounds must be finite")
    maximum = max(float(np.max(np.abs(bounds))), np.finfo(float).tiny)
    exponent = int(np.floor(np.log10(maximum)))
    scale = 10.0 ** exponent
    low, high = min(float(np.min(bounds / scale)), 0), max(float(np.max(bounds / scale)), 0)
    padding = max((high - low) * 0.18, 0.05)
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(width_inches, 6.5))
        gs = fig.add_gridspec(4, 2, height_ratios=[1.14, 1, 1, 1],
                              left=0.105, right=0.975, bottom=0.06, top=0.955,
                              hspace=0.77, wspace=0.19)
        schematic_axes = []
        for col, partition in enumerate(["std", "rot"]):
            # The left schematic aligns with the panels below and the right one
            # stays centered, which leaves room for the legend between them.
            ax = fig.add_subplot(gs[0, col], anchor="W" if col == 0 else "C")
            schematic_axes.append(ax)
            _schematic(ax, data["schematics"], partition,
                       data["partitions"][partition]["labels"])
            ax.set_title("Standard quadrants" if partition == "std" else r"$45^\circ$-rotated partition",
                         pad=5)
        for row, config in enumerate(configurations):
            for col, partition in enumerate(["std", "rot"]):
                panel = panels[(config, partition)]
                ax = fig.add_subplot(gs[row + 1, col])
                coefficients = np.asarray(panel["coefficients"][illustrated], float) / scale
                lo = np.asarray(panel["profile_lo"][illustrated], float) / scale
                hi = np.asarray(panel["profile_hi"][illustrated], float) / scale
                errors = np.array([coefficients - lo, hi - coefficients])
                if np.any(errors < 0):
                    raise ValueError("profile bounds do not bracket regional estimates")
                ax.axhline(0, color="#b7b7b7", lw=0.8, ls=":", zorder=1)
                ax.axhline(float(panel["common"][illustrated]) / scale,
                           color=COMMON, lw=1.0, zorder=2)
                ax.errorbar(np.arange(4), coefficients, yerr=errors,
                            fmt="o", markersize=4.5, color="#242424", elinewidth=1,
                            capsize=2.7, capthick=0.85, zorder=3)
                ax.set_xticks(np.arange(4), panel["labels"])
                ax.set_xlim(-0.45, 3.45)
                # Headroom keeps the one-line annotation above every point.
                ax.set_ylim(low - padding, high + padding * 2.2)
                # No top ticks: they would cut through the annotation line.
                ax.tick_params(axis="x", pad=5, top=False)
                label = chr(ord("a") + row * 2 + col)
                ax.text(0.025, 0.965, f"({label})", transform=ax.transAxes,
                        ha="left", va="top", fontsize=10)
                calibration = panel.get("calibration")
                annotation = rf"$Q_6={_number(panel['q_restricted'])}$"
                if calibration is not None and calibration.get("exceedances") is not None:
                    annotation += f"; {calibration['exceedances']}/{calibration['n_skies']}"
                    if not calibration.get("matched_to_data_row", False):
                        annotation += " reference null"
                else:
                    annotation += "; calibration unavailable"
                annotation = annotation.replace("; ", ", ")
                ax.text(0.975, 0.965, annotation, transform=ax.transAxes,
                        ha="right", va="top", fontsize=10)
                if col:
                    ax.tick_params(labelleft=False)
            # One foreground description spans each pair without crowding axes.
            pair_y = gs[row + 1, 0].get_position(fig).y1
            fig.text(0.54, pair_y + 0.012, row_titles[row], ha="center", va="bottom", fontsize=10)
        energy = np.sqrt(np.prod(np.asarray(data["energy_edges_gev"], float)[illustrated:illustrated + 2]))
        lower_panels = gs[1:, 0].get_position(fig)
        fig.text(0.026, (lower_panels.y0 + lower_panels.y1) / 2,
                 rf"NFW coefficient at {energy:.1f} GeV [$10^{{{exponent}}}$]",
                 ha="center", va="center", rotation=90, fontsize=10)
        handles = [Line2D([], [], color="#242424", marker="o", linestyle="none",
                          markersize=4.5, label="Regional fits,\n" r"68\% profile intervals"),
                   Line2D([], [], color=COMMON, lw=1, label="Common fit to\nthe same data")]
        for ax in schematic_axes:
            ax.apply_aspect()
        left_box, right_box = (ax.get_position() for ax in schematic_axes)
        # The right schematic's latitude label and ticks extend about 0.07 left.
        fig.legend(handles=handles, loc="center",
                   bbox_to_anchor=((left_box.x1 + right_box.x0 - 0.07) / 2, (left_box.y0 + left_box.y1) / 2),
                   ncol=1, frameon=False, handlelength=2.3, labelspacing=0.6)
        return _save(fig, output_dir, "regional_nfw_normalizations", test_data)
