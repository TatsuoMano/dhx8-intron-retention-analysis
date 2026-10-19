#!/usr/bin/env python3
"""Re-plot the graphs that were pasted from GraphPad Prism into MainFigures/SuppFigures .ai.

The .ai files contain these graphs only as pasted vector art (no linked file), so a
figcomp spec has nothing to link. This script reads the Prism .pzfx sources in
Manuscript/Figures/figures_for_final and writes one PDF+PNG per graph (figprov sidecar
included) so the composite can be laid out. Statistics annotations (p-values, brackets)
that Prism computed are NOT reproduced here; only the data (mean ± SD, points/lines).

Usage: python3 prism_panels.py [--out panels]
"""
import argparse, os, sys, math
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.expanduser("~/bin"))
import figprov

FFF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "Manuscript", "Figures", "figures_for_final")
FFF = os.path.normpath(FFF)
N_CELLS_1C = [520, 338, 876, 534]   # cells per condition (WT siCont, WT siDhx8, KO siCont, KO siDhx8), from the Fig 1C legend
FS = 7  # pt (labels); tick labels TS
TS = 6
STATS = {}
# v2 cell sizes (174 mm NAR layout, figcomp cells): overrides the per-function defaults when --v2 is given
V2 = {"fig1c_sg_fraction": (1.66, 1.9), "fig1e_peif2a_quant": (2.53, 2.5), "fig2e_lineprofile_p62": (2.18, 2.3),
      "fig8b_tdp43_quant": (1.57, 1.5), "fig8b_tdp35_quant": (1.57, 1.5), "s2b_sg_timecourse": (2.46, 1.8),
      "s5b_lineprofile_edc4": (2.46, 1.8), "s7a_tdp43_qpcr": (1.18, 1.19), "s7b_tdp43_quant": (1.18, 1.19), "s7c_chx_chase": (1.73, 1.19)}
# cell sizes of the v2 specs: Figure1 width 161 mm (u 18.4 mm), Figure8 166 mm (u 19), Figure2/S2/S5 170 mm (u 19.5), S7 174 mm / 11 cols
USE_V2 = False


def _ns(root):
    return root.tag[:root.tag.index("}") + 1] if root.tag.startswith("{") else ""


def read_pzfx(path):
    """Return {table_title: {"x": [...], "rowtitles": [...], "cols": [(name, [[subcol values]...])]}}."""
    root = ET.parse(path).getroot()
    ns = _ns(root)
    out = {}
    for tbl in root.iter(ns + "Table"):
        title = "".join(tbl.find(ns + "Title").itertext())
        def col_vals(col):
            subs = []
            for s in col.findall(ns + "Subcolumn"):
                vals = []
                for d in s.findall(ns + "d"):
                    t = "".join(d.itertext()).strip()
                    vals.append(float(t) if t not in ("", "*") else np.nan)
                subs.append(vals)
            return subs
        xc = tbl.find(ns + "XColumn")
        x = [v for s in col_vals(xc) for v in s] if xc is not None else []
        rt = tbl.find(ns + "RowTitlesColumn")
        rowtitles = ["".join(d.itertext()) for s in rt.findall(ns + "Subcolumn") for d in s.findall(ns + "d")] if rt is not None else []
        cols = []
        for yc in tbl.findall(ns + "YColumn"):
            name = "".join(yc.find(ns + "Title").itertext()) if yc.find(ns + "Title") is not None else ""
            cols.append((name, col_vals(yc)))
        key = title
        i = 2
        while key in out:            # duplicate table titles (e.g. "SGs Aria" x2)
            key = "%s#%d" % (title, i); i += 1
        out[key] = {"x": x, "rowtitles": rowtitles, "cols": cols}
    return out


def oneway(table):
    """OneWay table -> [(name, values)]"""
    return [(n, np.array([v for v in subs[0] if not np.isnan(v)])) for n, subs in table["cols"]]


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=TS, width=0.6, length=2.5)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.6)


def err(v, kind):
    s = v.std(ddof=1) if len(v) > 1 else 0.0
    return s / np.sqrt(len(v)) if kind == "sem" else s


def bars(ax, groups, colors, ylabel, xticklabels=None, rot=90, width=0.65, points=True, ylim=None, errkind="sem"):
    xs = np.arange(len(groups))
    for i, ((name, v), c) in enumerate(zip(groups, colors)):
        ax.bar(i, v.mean(), width=width, color=c, edgecolor="black", linewidth=0.6,
               yerr=err(v, errkind), capsize=2, error_kw={"lw": 0.6, "capthick": 0.6})
        if points:
            jit = (np.arange(len(v)) - (len(v) - 1) / 2) * 0.12
            ax.plot(i + jit, v, "o", ms=1.8, color="black", mec="black", mfc="black", zorder=3)
    ax.set_xticks(xs)
    ax.set_xticklabels(xticklabels or [g[0] for g in groups], rotation=rot, fontsize=FS,
                       ha="right" if rot not in (0, 90) else "center")
    ax.set_ylabel(ylabel, fontsize=FS)
    if ylim:
        ax.set_ylim(*ylim)
    style(ax)


def fmtp(p):
    return "<0.0001" if p < 0.0001 else (">0.9999" if p > 0.9999 else "%.4f" % p)   # Prism's printed-P convention


def below_ticks(fig, ax, pad_pt=3.0):
    """y (axes fraction) just below the lowest x tick label; needs a draw for the label extents."""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = ax.get_xaxis_transform().inverted()
    bbs = [t.get_window_extent(r) for t in ax.get_xticklabels() if t.get_text()]
    lo = min(inv.transform((bb.x0, bb.y0))[1] for bb in bbs) if bbs else -0.05
    return lo - pad_pt / (ax.get_window_extent(r).height / fig.dpi * 72.0)


def brackets(ax, pairs, ymax, step=0.09, fs=None):
    """pairs: [(i, j, p)] drawn as brackets above ymax (axis units). Narrow spans first; brackets that do not overlap
    horizontally share a row (greedy packing), so 14 comparisons need ~5 rows instead of 14."""
    fs = fs or FS
    lo, hi = ax.get_ylim()
    rows, placed = [], []                        # each row: list of (i, j) occupied spans
    for i, j, p in sorted(pairs, key=lambda t: (abs(t[1] - t[0]), t[0])):
        a, b = min(i, j), max(i, j)
        for r, spans in enumerate(rows):
            if all(b + 1.3 < sa or a - 1.3 > sb for sa, sb in spans):       # 1.3 x-unit clearance so that the P texts do not touch
                spans.append((a, b)); break
        else:
            rows.append([(a, b)]); r = len(rows) - 1
        placed.append((i, j, p, r))
    # row step: at least `step` of the data range, and at least enough printed height for a P text to clear the next row's
    # line (va="bottom" puts the descent at y + 0.2 h, so the digits reach y + 0.2 h + 0.95 fs; + 1 pt clearance; the brackets of
    # the next row start at y + h with 0.15 h end ticks). The axes height is measured after a
    # draw; the new y range (which the step itself enlarges) is solved in closed form, then checked once more after the redraw.
    need_pt = (0.95 * fs + 1.0) / 0.8          # the text must also clear the end ticks of the brackets one row up (y + h)
    R = len(rows)
    h = (hi - lo) * step
    for _ in range(2):
        ax.figure.canvas.draw()
        ax_h = ax.get_window_extent().height * 72.0 / ax.figure.dpi
        base = (ymax - lo)
        if ax_h > need_pt * (R + 1.0):
            h = max(h, need_pt * base / (ax_h - need_pt * (R + 1.0)))
        ax.set_ylim(lo, max(hi, ymax + h * (R + 1.0)))
    y0 = ymax + h * 0.6
    for i, j, p, r in placed:
        y = y0 + r * h
        ax.plot([i, i, j, j], [y, y + h * 0.15, y + h * 0.15, y], color="black", lw=0.6, clip_on=False)
        ax.text((i + j) / 2, y + h * 0.2, fmtp(p), ha="center", va="bottom", fontsize=fs)
    ax.set_ylim(lo, max(hi, y0 + R * h + h * 0.4))
    if R >= 3:   # a tall bracket stack (Fig 1C, 1E): the axis line and ticks stop at the first tick above the data, the brackets
        #          sit in the space above it (the tick step is the one of the enlarged range, so the labels keep their step)
        tk = [t for t in ax.get_yticks() if lo - 1e-9 <= t <= ax.get_ylim()[1] + 1e-9]
        top = min([t for t in tk if t >= ymax - 1e-9] or [max(tk)])
        ax.set_yticks([t for t in tk if t <= top + 1e-9])
        ax.spines["left"].set_bounds(lo, top)


def _stats():
    try:
        import prism_stats
        return prism_stats.compute()
    except Exception as exc:  # noqa: BLE001
        print("stats unavailable:", exc); return {}


BUMP = {}   # Figure 1 cells now match the panel sizes (cols 7); no bump needed   # Figure 1 rows are solved ~0.93 x cell; keep >= 6/7 pt after placement


def save(fig, stem, inputs, note, figsize):
    if USE_V2 and stem in V2:
        figsize = V2[stem]; fig.set_size_inches(*figsize)
        if stem in BUMP:
            import matplotlib.text
            for t in fig.findobj(matplotlib.text.Text):
                if t.get_text(): t.set_fontsize(t.get_fontsize() * BUMP[stem])
            for ax in fig.axes:   # tick labels are re-created at draw time from tick_params, so bump them there too
                ax.tick_params(labelsize=TS * BUMP[stem])
    figprov.save_figure(fig, os.path.join(OUT, stem), inputs=inputs, script=__file__,
                        note=note + " | figsize=%s" % (figsize,), savefig_kwargs={"bbox_inches": None})
    plt.close(fig)


# ---------------------------------------------------------------------------
def fs_of(stem, default):
    """figsize at creation (the v2 cell size when --v2), so that layout-dependent placements (below_ticks) see the final size."""
    return V2[stem] if USE_V2 and stem in V2 else default


def fig1c():
    src = os.path.join(FFF, "Fig01", "Fig1c.pzfx")
    t = read_pzfx(src)["Data 3"]
    conds = t["rowtitles"]                        # siCont./siDhx8 (WT), siCont./siDhx8 (Sidt2 KO)
    sg = np.array([[s[r] for s in t["cols"][0][1]] for r in range(4)])      # 4 rows x 10 fields
    df = np.array([[s[r] for s in t["cols"][1][1]] for r in range(4)])
    figsize = fs_of("fig1c_sg_fraction", (1.7, 2.1))
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    x = np.arange(4)
    m_sg = np.nanmean(sg, 1); s_sg = np.nanstd(sg, 1, ddof=1) / np.sqrt(np.sum(~np.isnan(sg), 1)); m_df = np.nanmean(df, 1)   # SEM (legend)
    ax.bar(x, m_sg, 0.65, color="#FFE600", edgecolor="black", linewidth=0.6, label="SGs",
           yerr=s_sg, capsize=2, error_kw={"lw": 0.6, "capthick": 0.6})
    ax.bar(x, m_df, 0.65, bottom=m_sg, color="white", edgecolor="black", linewidth=0.6, label="diffuse",
           yerr=s_sg, capsize=2, error_kw={"lw": 0.6, "capthick": 0.6})
    n = [int(np.sum(~np.isnan(sg[r]))) for r in range(4)]
    for i, ncell in enumerate(N_CELLS_1C):        # cell counts (legend of Fig 1C), written inside the bars as in the current figure
        ax.text(i, 50, "n=%d" % ncell, rotation=90, ha="center", va="center", fontsize=FS)
    ax.set_xticks(x); ax.set_xticklabels(conds, rotation=90, fontsize=FS)
    ax.set_ylabel("% of cells with G3BP1 SGs\nor diffuse signal", fontsize=FS)
    ax.set_ylim(0, 105)
    xt = ax.get_xaxis_transform(); YG = below_ticks(fig, ax, 5.0)               # just below the rotated tick labels
    for xc, lab, (a, b) in ((0.5, "WT", (-0.35, 1.35)), (2.5, "Sidt2 KO", (1.65, 3.35))):
        ax.text(xc, YG, lab, ha="center", va="top", transform=xt, fontsize=FS)
        ax.plot([a, b], [YG + 0.02, YG + 0.02], color="black", lw=0.6, transform=xt, clip_on=False)
    ax.legend(fontsize=FS, frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, handlelength=1, columnspacing=0.8)
    st = STATS.get("1C", {})
    if st:
        brackets(ax, [(0, 1, st["WT: siCont vs siDhx8"]), (2, 3, st["KO: siCont vs siDhx8"]), (0, 2, st["siCont: WT vs KO"]), (1, 3, st["siDhx8: WT vs KO"])], 105)
    style(ax)
    save(fig, "fig1c_sg_fraction", [src], "Fig 1C: %% cells with G3BP1 SGs vs diffuse; fields n=%s; P = two-way ANOVA + Šídák (prism_stats.py)" % (n,), figsize)


def fig1e():
    src = os.path.join(FFF, "Fig01", "Fig1f. P-EIF2a定量_240112 .pzfx")
    t = read_pzfx(src)["P-EIF2α"]
    g = oneway(t)      # 8 columns: (AS-: WT -,+ ; KO -,+ ; AS+: WT -,+ ; KO -,+)
    figsize = fs_of("fig1e_peif2a_quant", (2.3, 2.1))
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    cols = ["white", "#FFBFDE", "#D5D5D5", "#FF89D7", "#D5D5D5", "#FD5BA3", "#D5D5D5", "#FD5BA3"]   # master (.ai) bar fills, left to right
    bars(ax, g, cols, "Relative P-EIF2α levels\n(% of siCont.)", xticklabels=["–", "+"] * 4, rot=0)
    xt = ax.get_xaxis_transform(); Y1, Y2, Y3 = -0.06, -0.19, -0.34       # rows below the axis (axes fraction)
    ax.text(-0.75, Y1, "siDhx8:", ha="right", va="top", transform=xt, fontsize=FS)
    ax.text(-0.75, Y2, "Sidt2:", ha="right", va="top", transform=xt, fontsize=FS)
    for i, lab in enumerate(["WT", "KO", "WT", "KO"]):
        ax.text(i * 2 + 0.5, Y2, lab, ha="center", va="top", transform=xt, fontsize=FS)
        ax.plot([i * 2 - 0.35, i * 2 + 1.35], [Y2 - 0.13, Y2 - 0.13], color="black", lw=0.6, transform=xt, clip_on=False)
    ax.text(1.5, Y3, "AS (–)", ha="center", va="top", transform=xt, fontsize=FS)
    ax.text(5.5, Y3, "AS (+)", ha="center", va="top", transform=xt, fontsize=FS)
    ax.tick_params(axis="x", pad=1.5)
    st = STATS.get("1E", {})
    if st:
        N = ["AS-/WT/-", "AS-/WT/+", "AS-/KO/-", "AS-/KO/+", "AS+/WT/-", "AS+/WT/+", "AS+/KO/-", "AS+/KO/+"]
        PAIRS = [(0, 1), (0, 2), (0, 3), (2, 3), (1, 3), (4, 5), (4, 6), (4, 7), (5, 7), (6, 7), (0, 4), (1, 5), (2, 6), (0, 7)]
        ymax = max(v.max() + v.std(ddof=1) for _, v in g)
        brackets(ax, [(i, j, st["%s vs %s" % (N[i], N[j])]["tukey"]) for i, j in PAIRS], ymax, step=0.095)
    save(fig, "fig1e_peif2a_quant", [src], "Fig 1E: relative P-EIF2α (% of siCont.), n=3 each; P = one-way ANOVA + Tukey HSD (matches the printed values)", figsize)


def lineprofile(src, key, stem, note, labels, colors):
    t = read_pzfx(src)[key]
    x = np.array(t["x"])
    figsize = (2.3, 1.7)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    for (name, subs), lab, c in zip(t["cols"], labels, colors):
        y = np.array(subs[0])[: len(x)]
        ax.plot(x, y, color=c, lw=0.9, label=lab)
    ax.set_xlabel("Distance (µm)", fontsize=FS); ax.set_ylabel("Intensity (%)", fontsize=FS)
    ax.set_xlim(0, 4); ax.set_ylim(0, 100); ax.set_xticks([0, 1, 2, 3, 4]); ax.set_yticks([0, 25, 50, 75, 100])
    ax.legend(fontsize=FS, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, handlelength=1.5)
    style(ax)
    save(fig, stem, [src], note, figsize)


def fig8b():
    src = os.path.join(FFF, "Fig08", "Fig8b4. Fujikake-_定量データ_250716.pzfx")
    T = read_pzfx(src)
    for key, stem, yl in [("TDP-43", "fig8b_tdp43_quant", "TDP-43/Tubulin"), ("TDP-35", "fig8b_tdp35_quant", "TDP-35/Tubulin")]:
        g = oneway(T[key])
        figsize = (1.15, 1.6)
        fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
        bars(ax, g, ["white", "#F5B27A", "#F5B27A"], yl, ylim=(0, 3), errkind="sd")   # legend: mean ± SD
        ax.set_xlabel("RNAi:", fontsize=FS)
        st = STATS.get("8B", {}).get(key)
        if st:
            brackets(ax, [(0, 1, st["dDHX8 #1"]), (0, 2, st["dDHX8 #2"])], max(v.max() + v.std(ddof=1) for _, v in g), step=0.1)
        save(fig, stem, [src], "Fig 8B: %s in Drosophila, n=3; P = one-way ANOVA + Dunnett vs Cont. (prism_stats.py)" % key, figsize)


def s2b():
    src = os.path.join(FFF, "SupFig02", "SupFig02a. 2SGs定量_230609IC.pzfx")
    t = read_pzfx(src)["SGs number"]
    x = np.array(t["x"])
    figsize = (2.4, 1.8)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    for (name, subs), c in zip(t["cols"], ["#9A9A9A", "#F28C28"]):
        a = np.array(subs)            # replicates x time
        m = np.nanmean(a, 0); s = np.nanstd(a, 0, ddof=1) / np.sqrt(np.sum(~np.isnan(a), 0))   # SEM (matches the master bars)
        ax.errorbar(x, m, yerr=s, color=c, lw=0.9, marker="o", ms=2, capsize=1.5, elinewidth=0.6, label=name)
    ax.set_xlabel("Time (min)", fontsize=FS); ax.set_ylabel("SGs (number/cell)", fontsize=FS)
    ax.set_xlim(0, 21); ax.set_ylim(0, 80)
    ax.legend(fontsize=FS, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, handlelength=1.5)
    style(ax)
    save(fig, "s2b_sg_timecourse", [src], "Fig S2B: SG number per cell after arsenite, mean ± SEM; the printed P (0.9154 / 0.0003) could not be reproduced by t-tests T=1 vs 20 min, test unknown, not annotated", figsize)


def s7():
    src_a = os.path.join(FFF, "SupFig12", "SupFig12a._DHX8KD_tdp43_actb_gapdh まとめ_250306.pzfx")
    g = oneway(read_pzfx(src_a)["TDP-43/ActB qPCR (DHX8 #1のみ)"])
    figsize = (1.1, 1.8)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    bars(ax, g, ["white", "#F28C28"], "Tdp-43/Actb\nexpression", ylim=(0, 10))
    if "S7A" in STATS: brackets(ax, [(0, 1, STATS["S7A"])], max(v.max() + v.std(ddof=1) for _, v in g), step=0.1)
    save(fig, "s7a_tdp43_qpcr", [src_a], "Fig S7A: Tdp-43/Actb qPCR, n=3; P = unpaired t-test", figsize)

    src_bc = os.path.join(FFF, "SupFig12", "SupFig12b&c_231222 定量.pzfx")
    T = read_pzfx(src_bc)
    g = oneway(T["T=0 TDP-43"])
    figsize = (1.1, 1.8)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    bars(ax, g, ["white", "#F28C28"], "TDP-43 level\n(% of siCont.)", ylim=(0, 200))
    if "S7B" in STATS: brackets(ax, [(0, 1, STATS["S7B"])], max(v.max() + v.std(ddof=1) for _, v in g), step=0.1)
    save(fig, "s7b_tdp43_quant", [src_bc], "Fig S7B: TDP-43 protein (% of siCont.), n=3; P = unpaired t-test", figsize)

    t = T["TDP-43"]; x = np.array(t["x"])
    figsize = (2.2, 1.8)
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    for (name, subs), c in zip(t["cols"], ["black", "#F28C28"]):
        a = np.array(subs)   # error bars: SD (the siDhx8 whiskers of the current figure are SD-sized; convention not stated in the legend)
        ax.errorbar(x, np.nanmean(a, 0), yerr=np.nanstd(a, 0, ddof=1), color=c, lw=0.9, marker="o", ms=2.5, capsize=1.5, elinewidth=0.6, label=name)
    ax.set_xlabel("Time (h)", fontsize=FS); ax.set_ylabel("TDP-43 level\n(% of 0 h)", fontsize=FS)
    ax.set_xticks([0, 8, 16, 24]); ax.set_ylim(0, 150)
    ax.legend(fontsize=FS, frameon=False, loc="lower left", handlelength=1.5)
    style(ax)
    save(fig, "s7c_chx_chase", [src_bc], "Fig S7C: CHX chase, TDP-43 % of 0 h, mean ± SD n=3 (the siDhx8 whiskers of the current figure are SD-sized; siCont whiskers look SEM-sized, convention unknown)", figsize)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="panels")
    ap.add_argument("--v2", action="store_true", help="v2 cell sizes (174 mm layout), tick labels 6 pt")
    a = ap.parse_args()
    USE_V2 = a.v2
    if a.v2:
        TS = 6.2   # tick labels stay >= 6 pt even after the small solve-mode scaling
    OUT = os.path.abspath(a.out); os.makedirs(OUT, exist_ok=True)
    figprov.apply_arial()
    STATS = _stats()
    fig1c(); fig1e()
    lineprofile(os.path.join(FFF, "Fig02", "Fig2e. Line Plofile_231225 IC .pzfx"),
                "3. siDHX8, AS(-), G3BP1-GFP, TDP-43-delNLS-mChe, p62", "fig2e_lineprofile_p62",
                "Fig 2E: line profile G3BP1-GFP / TDP-43-ΔNLS-mChe / p62", ["G3BP1-GFP", "TDP-43-ΔNLS-mChe", "p62"], ["#2ECC40", "#FF00FF", "#00CFFF"])
    lineprofile(os.path.join(FFF, "SupFig14", "SupFig14B. Line Plofile_231225 IC .pzfx"),
                "6. siDHX8, G3BP1-GFP, TDP-43-delNLS-mChe, EDC4", "s5b_lineprofile_edc4",
                "Fig S5B: line profile G3BP1-GFP / TDP-43-ΔNLS-mChe / EDC4", ["G3BP1-GFP", "TDP-43-ΔNLS-mChe", "EDC4"], ["#2ECC40", "#FF00FF", "#00CFFF"])
    fig8b(); s2b(); s7()
    print("done ->", OUT)
