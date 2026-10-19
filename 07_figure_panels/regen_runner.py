#!/usr/bin/env python3
"""Re-render an existing figure script at a given figsize with a font floor, WITHOUT editing the script.

The original script is executed as-is (its own data loading, styling and text), but matplotlib is patched so that
  * every savefig / figprov.save_figure call is intercepted (the original output paths are never written),
  * the figure selected by --select (substring of the original output name; default: all) is resized to --figsize,
    optionally reduced to one of its axes (--panel i, the other axes are removed), text sizes are clamped
    (tick labels >= --floor-tick, other text >= --floor-label, everything <= --cap), the font family is forced to Arial,
    the layout is recomputed (constrained layout, or tight_layout on failure),
  * the result is written with figprov.save_figure to --out (or --out-dir/<original basename>) with a sidecar that
    records the script, the files it read (open/pandas hooks), the figsize and the clamp settings.

    python3 regen_runner.py SCRIPT --figsize W H --out STEM [--select SUBSTR] [--panel I] [--cwd DIR]
                             [--floor-label 7] [--floor-tick 6] [--cap 8] [--title-off] [--out-dir DIR] [--ylim LO HI] [--xlim LO HI] [--text-move PATTERN DX DY] [--text-offset-pt PATTERN DX DY] [--clear-xaxis [PAD_PT]] [--keep-look] [--keep-ticks] [--axes-rect L B R T] [--leader-to-edge] [--legend-ncol N] [--tick-replace OLD NEW] [--sup-mathtext]

Runs on the workstation or on the server (needs ~/bin/figprov.py). The panel-regeneration contract of the
figure-composite skill (references/panel-regeneration.md) is satisfied: <stem>.pdf/.png/.pdf.prov.json, bbox_inches=None.
"""
import argparse, builtins, os, runpy, sys, io, contextlib, traceback

sys.path.insert(0, os.path.expanduser("~/bin"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.text
import matplotlib.legend
from matplotlib.figure import Figure
import figprov

READ_FILES = []
_real_open = builtins.open
_real_save_figure = figprov.save_figure
_real_savefig = Figure.savefig
_real_plt_savefig = plt.savefig


def _rec(path, mode):
    try:
        if isinstance(path, (str, bytes, os.PathLike)) and ("r" in str(mode)) and ("w" not in str(mode)):
            p = os.path.abspath(os.fspath(path))
            if os.path.isfile(p) and p not in READ_FILES and not p.endswith((".py", ".pyc")):
                READ_FILES.append(p)
    except Exception:  # noqa: BLE001
        pass


def _open(path, mode="r", *a, **k):
    _rec(path, mode)
    return _real_open(path, mode, *a, **k)


def all_texts(fig):
    """Every Text artist, including matplotlib table cells (their Text objects are not children of the Table)."""
    import matplotlib.table
    ts = list(fig.findobj(matplotlib.text.Text))
    for tb in fig.findobj(matplotlib.table.Table):
        ts.extend(c.get_text() for c in tb.get_celld().values())
    return ts


def clamp_text(fig, floor_label, floor_tick, cap, family="Arial"):
    ticks = set()
    for ax in fig.axes:
        for lab in ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"):
            ticks.add(id(lab))
        for t in getattr(ax, "texts", []):
            pass
    import re as _re
    n = 0
    for t in all_texts(fig):
        if not t.get_text():
            continue
        s = t.get_fontsize()
        # tick labels that are names (pathways, genes, conditions) are read as labels: only numeric tick labels get the tick floor
        numeric = bool(_re.fullmatch(r"[\s\d.,+\-\u2212eE%$^{}\\a-z]*", t.get_text())) and not _re.search(r"[A-Za-df-z]{2,}", t.get_text())
        numeric = numeric or "mathdefault" in t.get_text()      # 10^x tick labels of log axes
        floor = floor_tick if (id(t) in ticks and numeric) else floor_label
        new = min(max(s, floor), cap)
        if abs(new - s) > 0.01:
            t.set_fontsize(new); n += 1
        try:
            t.set_fontfamily(family)
        except Exception:  # noqa: BLE001
            pass
    return n


def scale_artists(fig, r, min_lw=0.4):
    """Scale everything that is measured in points by r, as shrinking the finished figure would: line widths, marker sizes,
    scatter sizes, patch edges, tick marks (they are Line2D markers) and text. Lines keep at least min_lw."""
    import matplotlib.lines as mlines, matplotlib.collections as mcoll, matplotlib.patches as mpatches
    import numpy as np
    n = 0; seen = set()
    arts = list(fig.findobj())
    for lg in fig.findobj(matplotlib.legend.Legend):
        arts.extend(lg.findobj())
    for art in arts:
        if id(art) in seen: continue
        seen.add(id(art))
        try:
            if isinstance(art, mlines.Line2D):
                lw = art.get_linewidth()
                if lw: art.set_linewidth(max(min_lw, lw * r))
                art.set_markersize(max(1.2, art.get_markersize() * r) if art.get_marker() not in (None, "None", "", " ") and art.get_markersize() > 3 else art.get_markersize() * r)
                art.set_markeredgewidth(art.get_markeredgewidth() * r)
                n += 1
            elif isinstance(art, mcoll.Collection):
                sz = np.asarray(art.get_sizes(), dtype=float) if hasattr(art, "get_sizes") else np.array([])
                if sz.size: art.set_sizes(np.maximum(sz * r * r, 1.2))
                lws = np.asarray(art.get_linewidths(), dtype=float)
                if lws.size: art.set_linewidths([0 if w == 0 else max(min_lw * 0.6, w * r) for w in lws])
                n += 1
            elif isinstance(art, mpatches.Patch):
                lw = art.get_linewidth()
                if lw: art.set_linewidth(max(min_lw, lw * r))
                n += 1
        except Exception:  # noqa: BLE001
            pass
    for t in all_texts(fig):
        try:
            t.set_fontsize(t.get_fontsize() * r); n += 1
        except Exception:  # noqa: BLE001
            pass
    return n


def _legend_out(fig, a):
    """Legends below the axes as figure legends, each centred under its own axes (several axes with legends would
    otherwise be superimposed at the same anchor); bottom margin reserved explicitly (constrained layout does not
    reserve space for a legend anchored outside the axes, so it would fall off the page)."""
    pend = []; ncol_max = 1
    for ax in fig.axes:
        lg = ax.get_legend()
        if lg is None: continue
        h = list(getattr(lg, "legend_handles", None) or getattr(lg, "legendHandles", []))   # keep custom handles (Patch lists)
        l = [t.get_text() for t in lg.get_texts()]
        if not h or len(h) != len(l):
            h, l = ax.get_legend_handles_labels()
        fs = lg.get_texts()[0].get_fontsize() if lg.get_texts() else a.floor_label
        lg.remove()
        ncol = a.legend_ncol or (2 if len(l) >= 3 else max(1, len(l))); ncol_max = max(ncol_max, ncol)     # >= 3 entries: two columns (fits 2.2 in panels)
        pend.append((ax, h, l, fs, ncol))
    if pend:
        try:
            import math as _m
            fig.set_layout_engine("none")
            # room for the legend: its rows at the legend font size plus a gap, as a fraction of the page height
            rows = max(_m.ceil(len(l) / float(nc)) for _, _, l, _, nc in pend); fsz = max(fs for _, _, _, fs, _ in pend)
            Hpt = fig.get_size_inches()[1] * 72.0
            need = (rows * fsz * 1.45 + 5.0) / Hpt
            st = getattr(fig, "_suptitle", None); top = 1.0
            if st is not None and st.get_text():      # tight_layout ignores the suptitle: keep its lines free
                top = 1.0 - ((st.get_text().count("\n") + 1) * st.get_fontsize() * 1.3 + 4.0) / Hpt
                st.set_y(1.0 - 2.0 / Hpt); st.set_va("top")
                try: st.set_in_layout(False)      # newer matplotlib counts the suptitle in tight_layout as well: reserve it once
                except Exception: pass  # noqa: E701
            fig.tight_layout(rect=(0, min(0.35, need), 1, top))
        except Exception:  # noqa: BLE001
            pass
        for ax, h, l, fs, ncol in pend:
            pos = ax.get_position(); xc = (pos.x0 + pos.x1) / 2 if len(fig.axes) > 1 else 0.5
            lg2 = fig.legend(h, l, loc="lower center", bbox_to_anchor=(xc, 0.0), ncol=ncol, fontsize=fs, frameon=False,
                             handlelength=1.2, columnspacing=1.0, borderaxespad=0.2)
            try:      # a legend wider than the space beside its axes is shifted back inside the page
                fig.canvas.draw()
                bb = lg2.get_window_extent(fig.canvas.get_renderer()); Wpx = fig.get_size_inches()[0] * fig.dpi
                d = 0.0
                if bb.x1 > Wpx - 2: d = (Wpx - 2 - bb.x1) / Wpx
                if bb.x0 + d * Wpx < 2: d = (2 - bb.x0) / Wpx
                if d: lg2.set_bbox_to_anchor((xc + d, 0.0), transform=fig.transFigure)
            except Exception:  # noqa: BLE001
                pass


def _clear_xaxis(fig, pad_pt):
    """Artists placed with ax.get_xaxis_transform() at y < 0 (group rules and group labels under the tick labels) keep their
    axes-fraction y when the panel shrinks, so they strike through the enlarged tick labels. Shift them all together (relative
    offsets preserved) until the topmost one lies pad_pt below the lowest tick label. Needs a renderer for the tick-label extents."""
    import numpy as np
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    for ax in fig.axes:
        xt = ax.get_xaxis_transform()
        arts = [t for t in ax.texts if t.get_transform() is xt and t.get_position()[1] < 0]
        lines = [l for l in ax.lines if l.get_transform() is xt and len(l.get_ydata()) and max(l.get_ydata()) < 0]
        if not arts and not lines:
            continue
        bbs = [t.get_window_extent(r) for t in ax.get_xticklabels() if t.get_text()]
        if not bbs:
            continue
        inv = xt.inverted()
        tick_bottom = min(inv.transform((bb.x0, bb.y0))[1] for bb in bbs)
        ax_h_pt = ax.get_window_extent(r).height / fig.dpi * 72.0
        limit = tick_bottom - pad_pt / ax_h_pt
        top = max([t.get_position()[1] for t in arts] + [max(l.get_ydata()) for l in lines])
        if top <= limit:
            continue
        d = limit - top
        for t in arts:
            x0, y0 = t.get_position(); t.set_position((x0, y0 + d))
        for l in lines:
            l.set_ydata(np.asarray(l.get_ydata(), dtype=float) + d)
        print("regen_runner: below-axis artists shifted by %.3f (axes fraction) to clear the tick labels" % d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("script")
    ap.add_argument("--figsize", nargs=2, type=float, required=True, metavar=("W_IN", "H_IN"))
    ap.add_argument("--out", help="output stem (single selected figure)")
    ap.add_argument("--out-dir", help="output directory (several figures; original basenames kept)")
    ap.add_argument("--select", default="", help="substring of the original output filename to keep (default: all)")
    ap.add_argument("--panel", type=int, help="keep only axes index I of the selected figure")
    ap.add_argument("--cwd", help="directory to run the script in (default: the script's directory)")
    ap.add_argument("--floor-label", type=float, default=7.0)
    ap.add_argument("--floor-tick", type=float, default=6.0)
    ap.add_argument("--cap", type=float, default=8.0)
    ap.add_argument("--title-off", action="store_true", help="remove axes titles and suptitle")
    ap.add_argument("--xrot", type=float, default=0.0, help="rotate x tick labels (degrees)")
    ap.add_argument("--suptitle-off", action="store_true", help="remove the figure suptitle only")
    ap.add_argument("--text-replace", nargs=2, action="append", metavar=("OLD", "NEW"), help="replace OLD by NEW in every text artist (repeatable)")
    ap.add_argument("--label-wrap", type=int, default=0, help="wrap axis labels longer than N characters onto two lines")
    ap.add_argument("--title-regex", nargs=2, metavar=("PATTERN", "REPL"), help="re.sub applied to every axes title")
    ap.add_argument("--legend-out", action="store_true", help="move axes legends below the axes")
    ap.add_argument("--ylim", nargs=2, type=float, metavar=("LO", "HI"), help="set the y range of the kept axes (e.g. to keep in-axes annotations inside the frame after enlarging them)")
    ap.add_argument("--text-move", nargs=3, action="append", metavar=("PATTERN", "DX", "DY"), help="shift every text whose content matches PATTERN by (DX, DY) in the text's own coordinates (repeatable)")
    ap.add_argument("--text-offset-pt", nargs=3, action="append", metavar=("PATTERN", "DX", "DY"), help="shift every text whose content matches PATTERN by (DX, DY) printed points, independent of its coordinate system (repeatable)")
    ap.add_argument("--clear-xaxis", nargs="?", const=3.0, type=float, metavar="PAD_PT", help="push artists drawn below the x axis in axes-fraction coordinates (group rules / group labels via get_xaxis_transform) down until they clear the tick labels by PAD_PT (default 3)")
    ap.add_argument("--xlim", nargs=2, type=float, metavar=("LO", "HI"), help="set the x range of the kept axes")
    ap.add_argument("--axes-ylim", nargs=3, action="append", type=float, metavar=("I", "LO", "HI"), help="set the y range of axes I only (repeatable), e.g. room for a legend above the data")
    ap.add_argument("--legend-loc", help="matplotlib loc for the axes legends (e.g. 'upper left'), when the enlarged legend covers data at its original place")
    ap.add_argument("--text-size", nargs=2, action="append", metavar=("PATTERN", "PT"), help="font size for the texts matching PATTERN, applied after the floors (repeatable)")
    ap.add_argument("--text-rotate", nargs=2, action="append", metavar=("PATTERN", "DEG"), help="rotation for the texts matching PATTERN (repeatable), e.g. labels inside narrow boxes")
    ap.add_argument("--text-remove", action="append", metavar="PATTERN", help="remove the texts matching PATTERN (titles drawn with fig.text that --title-off does not reach)")
    ap.add_argument("--axes-rect", nargs=4, type=float, metavar=("L", "B", "R", "T"), help="position of the kept axes as fractions of the page (single-axes panels): same-format panels get the same rectangle, so their axes line up and have the same size")
    ap.add_argument("--keep-ticks", action="store_true", help="keep the tick positions of the original drawing (a smaller panel would otherwise get fewer ticks)")
    ap.add_argument("--leader-to-edge", action="store_true", help="leader lines of annotations end at the near edge of their label (top or bottom centre) instead of aiming at its centre, so that they do not cross neighbouring labels")
    ap.add_argument("--legend-ncol", type=int, help="number of columns of a legend moved out with --legend-out")
    ap.add_argument("--tick-replace", nargs=2, action="append", metavar=("OLD", "NEW"), help="replace OLD by NEW in fixed tick labels (category names); \\n = line break")
    ap.add_argument("--dpi", type=float, default=600.0, help="resolution of rasterised artists in the PDF and of the PNG twin")
    ap.add_argument("--sup-mathtext", nargs="?", const="regular", choices=["regular", "bold", "italic"], help="write Unicode superscripts (10⁻⁵) as mathtext superscripts (10^{-5}) set in Arial (Arial has no superscript glyphs); 'bold' / 'italic' when the texts that carry them are bold / italic (a string with mathtext is set entirely in the mathtext font, which drops the weight and style of the text)")
    ap.add_argument("--sci-notation", action="store_true", help="write e-notation (3.3e-15, p = 7.3e-05) as 3.3×10 with the exponent as a mathtext superscript (leading zeros of the exponent dropped); sets the mathtext font like --sup-mathtext (regular unless --sup-mathtext says otherwise)")
    ap.add_argument("--sup-raise", type=float, help="with --sup-mathtext / --sci-notation: raise of a superscript's baseline as a fraction of what mathtext takes for the x-height (default 0.7). For Arial, mathtext reads the x-height from the font's PCLT table with a formula that evaluates to about the font size, so the default raises a superscript by 0.7 of the text size (4.9 pt on 7-pt text, far above the cap height); 0.33 gives Illustrator's default superscript shift (33 %% of the size)")
    ap.add_argument("--sup-min-pt", type=float, help="with --sup-mathtext: smallest printed size of a superscript (pt). mathtext shrinks a superscript to 70%% of its text (7 pt -> 4.9 pt); the shrink factor is raised so that the smallest text carrying a superscript prints it at this size (5.5 for a 7-pt text -> 0.786)")
    ap.add_argument("--keep-look", action="store_true", help="keep the proportions of the original drawing: line widths, marker sizes, tick marks and text are scaled by the change of size (as if the original figure had been shrunk to --figsize), then the text floors are applied")
    ap.add_argument("--look-scale", type=float, help="scale used by --keep-look (default: new size / original size of the figure, or of the kept axes with --panel)")
    ap.add_argument("--min-lw", type=float, default=0.4, help="smallest line width after --keep-look (pt)")
    ap.add_argument("--layout", default="constrained", choices=["constrained", "tight", "none"])
    ap.add_argument("--note", default="")
    a = ap.parse_args()
    if not a.out and not a.out_dir:
        ap.error("--out or --out-dir is required")
    script = os.path.abspath(a.script)
    cwd = a.cwd or os.path.dirname(script)
    if a.out: a.out = os.path.abspath(a.out)          # resolve before chdir into the script's directory
    if a.out_dir: a.out_dir = os.path.abspath(a.out_dir)
    figsize = tuple(a.figsize)
    done = []

    def handle(fig, requested):
        base = os.path.splitext(os.path.basename(str(requested)))[0]
        if a.select and a.select not in base:
            return
        if a.out and done:
            return   # single-output mode: first match only
        w0, h0 = fig.get_size_inches()
        if a.panel is not None:
            try:
                pos0 = fig.axes[a.panel].get_position(); w0, h0 = w0 * pos0.width / 0.78, h0 * pos0.height / 0.78   # axes + its share of the margins
            except Exception:  # noqa: BLE001
                pass
        look = (a.look_scale or min(figsize[0] / w0, figsize[1] / h0)) if a.keep_look else None
        if a.keep_ticks or a.tick_replace:
            import matplotlib.ticker as _mt
            try: fig.canvas.draw()
            except Exception: pass  # noqa: E701
            for ax in fig.axes:
                for axis in (ax.xaxis, ax.yaxis):
                    fmt = axis.get_major_formatter(); loc = axis.get_major_locator()
                    if a.tick_replace:      # works for fixed, category and function formatters: the drawn labels are frozen
                        locs0 = list(axis.get_majorticklocs()); labs0 = [t.get_text() for t in axis.get_majorticklabels()]
                        labs1 = list(labs0)
                        for old_, new_ in a.tick_replace:
                            labs1 = [x.replace(old_.replace("\\n", "\n"), new_.replace("\\n", "\n")) for x in labs1]
                        if labs1 != labs0 and len(locs0) == len(labs1):
                            axis.set_major_locator(_mt.FixedLocator(locs0)); axis.set_major_formatter(_mt.FixedFormatter(labs1))
                            loc = axis.get_major_locator()
                    if a.keep_ticks and not isinstance(loc, (_mt.FixedLocator, _mt.NullLocator)) and axis.get_scale() == "linear":
                        lo, hi = sorted(axis.get_view_interval())
                        locs = [float(v) for v in axis.get_majorticklocs() if lo - 1e-9 <= v <= hi + 1e-9]
                        if 2 <= len(locs) <= 12: axis.set_major_locator(_mt.FixedLocator(locs))
        if a.panel is not None:
            keep = fig.axes[a.panel]
            for ax in list(fig.axes):
                if ax is not keep:
                    ax.remove()
            try:   # let the surviving axes use the whole figure instead of its old gridspec slot
                import matplotlib.gridspec as _gs
                keep.set_subplotspec(_gs.GridSpec(1, 1, figure=fig)[0, 0])
            except Exception:  # noqa: BLE001
                pass
        for i_, lo_, hi_ in (a.axes_ylim or []):
            fig.axes[int(i_)].set_ylim(lo_, hi_)
        if a.ylim or a.xlim:
            for ax in fig.axes:
                if a.ylim: ax.set_ylim(*a.ylim)
                if a.xlim: ax.set_xlim(*a.xlim)
        if a.xrot:
            for ax in fig.axes:
                for lab in ax.get_xticklabels():
                    lab.set_rotation(a.xrot); lab.set_ha("right"); lab.set_rotation_mode("anchor")
        pass  # legend-out handled after text clamping (see _legend_out)
        if a.title_off:
            for ax in fig.axes:
                ax.set_title("")
            fig.suptitle("")
        if a.suptitle_off:
            fig.suptitle("")
        if a.title_regex:
            import re as _re
            pat, rep = a.title_regex
            for ax in fig.axes:
                ax.title.set_text(_re.sub(pat, rep, ax.get_title()))     # set_text keeps the weight; set_title would reset it
        if a.sup_mathtext or a.sci_notation:
            SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-", "⁺": "+"}
            import re as _re3
            SCI = _re3.compile(r"(?<![\w.])(\d+(?:\.\d+)?)[eE]([-+]?)0*(\d+)(?![\w.])")      # 3.3e-15, 7.3e-05 (not gene names or 1e5 inside words)
            for t in all_texts(fig):
                txt = t.get_text()
                if txt and "$" not in txt and any(c in SUP for c in txt):
                    txt = _re3.sub("[" + "".join(SUP) + "]+", lambda m: "$^{" + "".join(SUP[c] for c in m.group(0)) + "}$", txt); t.set_text(txt)
                if a.sci_notation and txt and "$" not in txt and SCI.search(txt):
                    t.set_text(SCI.sub(lambda m: m.group(1) + "×10$^{" + ("-" if m.group(2) == "-" else "") + m.group(3) + "}$", txt))
        if a.text_replace:
            for t in all_texts(fig):
                for old, new in a.text_replace:
                    old, new = old.replace("\\n", "\n"), new.replace("\\n", "\n")     # a literal \n on the command line = line break
                    if old in t.get_text(): t.set_text(t.get_text().replace(old, new))
        if a.text_move:
            import re as _re
            for pat, dx, dy in a.text_move:
                for t in all_texts(fig):
                    if t.get_text() and _re.search(pat, t.get_text()):
                        x0, y0 = t.get_position(); t.set_position((x0 + float(dx), y0 + float(dy)))
        if a.text_offset_pt:
            import re as _re
            from matplotlib.transforms import offset_copy as _oc
            for pat, dx, dy in a.text_offset_pt:
                for t in all_texts(fig):
                    if t.get_text() and _re.search(pat, t.get_text()):
                        t.set_transform(_oc(t.get_transform(), fig=fig, x=float(dx), y=float(dy), units="points"))
        if a.label_wrap:
            import textwrap as _tw
            for ax in fig.axes:
                for get, set_ in ((ax.get_xlabel, ax.set_xlabel), (ax.get_ylabel, ax.set_ylabel)):
                    lab = get()
                    if len(lab) > a.label_wrap and "\n" not in lab:
                        set_("\n".join(_tw.wrap(lab, a.label_wrap)))
        if look:
            n_scaled = scale_artists(fig, look, a.min_lw)
            print("regen_runner: keep-look scale %.3f (original %.2f x %.2f in), %d artists scaled" % (look, w0, h0, n_scaled))
        fig.set_size_inches(*figsize)
        if a.layout == "constrained":
            try:
                st = getattr(fig, "_suptitle", None)
                if st is not None and st.get_text():       # a suptitle placed by hand is ignored by the layout and runs into the axes titles
                    fig.suptitle(st.get_text(), fontsize=st.get_fontsize(), fontweight=st.get_fontweight())
                fig.set_layout_engine("constrained")
            except Exception:  # noqa: BLE001
                pass
        n = clamp_text(fig, a.floor_label, a.floor_tick, a.cap)
        if a.sup_min_pt and (a.sup_mathtext or a.sci_notation):     # after the floors: the shrink is set from the final size of the carrying texts
            carriers = [t.get_fontsize() for t in all_texts(fig) if "$" in (t.get_text() or "")]
            if carriers:
                import matplotlib._mathtext as _mt          # SHRINK_FACTOR is read at layout time by Char/Box.shrink()
                _mt.SHRINK_FACTOR = max(0.7, min(1.0, a.sup_min_pt / min(carriers)))
                print("regen_runner: mathtext superscripts at %.3f of their text (smallest carrier %.1f pt -> %.1f pt)" % (_mt.SHRINK_FACTOR, min(carriers), min(carriers) * _mt.SHRINK_FACTOR))
        import re as _re2
        for pat in (a.text_remove or []):
            for t in all_texts(fig):
                if t.get_text() and _re2.search(pat, t.get_text()): t.set_text("")
        for pat, pt in (a.text_size or []):
            for t in all_texts(fig):
                if t.get_text() and _re2.search(pat, t.get_text()): t.set_fontsize(float(pt))
        for pat, deg in (a.text_rotate or []):
            for t in all_texts(fig):
                if t.get_text() and _re2.search(pat, t.get_text()): t.set_rotation(float(deg))
        if a.legend_loc:
            for ax in fig.axes:
                lg = ax.get_legend()
                if lg is not None:
                    try: lg.set_loc(a.legend_loc)
                    except Exception: lg._loc = matplotlib.legend.Legend.codes.get(a.legend_loc, 0)  # noqa: E701,SLF001
        if a.legend_out:
            _legend_out(fig, a)
        if a.layout == "tight":
            try:
                fig.tight_layout()
            except Exception:  # noqa: BLE001
                pass
        if a.clear_xaxis is not None:
            _clear_xaxis(fig, a.clear_xaxis)
        if a.leader_to_edge:
            import matplotlib.text as _mtx
            for an in fig.findobj(_mtx.Annotation):
                if getattr(an, "arrowprops", None) is not None:
                    an._arrow_relpos = (0.5, 1.0) if an.get_va() == "top" else (0.5, 0.0) if an.get_va() == "bottom" else (0.5, 0.5)  # noqa: SLF001
        if a.axes_rect:
            try: fig.set_layout_engine("none")
            except Exception: pass  # noqa: E701
            L_, B_, R_, T_ = a.axes_rect
            for ax in fig.axes[:1]:
                ax.set_position([L_, B_, R_ - L_, T_ - B_])
        stem = a.out if a.out else os.path.join(a.out_dir, base)
        os.makedirs(os.path.dirname(os.path.abspath(stem)) or ".", exist_ok=True)
        Figure.savefig = _real_savefig; plt.savefig = _real_plt_savefig      # real writers while saving
        matplotlib.rcParams["savefig.dpi"] = a.dpi                           # rasterised artists (dense scatters) inside the PDF
        _real_save_figure(fig, stem, inputs=[f for f in READ_FILES if os.path.getsize(f) < 200_000_000],
                            script=script, note=("regen_runner: figsize=%s panel=%s floors label>=%g tick>=%g cap<=%g; %d texts clamped; original output %s. %s"
                                                 % (figsize, a.panel, a.floor_label, a.floor_tick, a.cap, n, base, a.note)).strip(),
                            dpi=int(a.dpi), savefig_kwargs={"bbox_inches": None})
        Figure.savefig = fake_savefig; plt.savefig = fake_plt_savefig
        try:
            import json as _json
            fig.canvas.draw()
            boxes = [[round(v, 4) for v in (ax.get_position().x0, 1 - ax.get_position().y1, ax.get_position().x1, 1 - ax.get_position().y0)] for ax in fig.axes]
            with _real_open(stem + ".axes.json", "w") as fh:     # [x0, y0, x1, y1] as fractions of the page, y down
                _json.dump({"figsize_in": list(figsize), "axes": boxes, "look_scale": look}, fh)
        except Exception:  # noqa: BLE001
            pass
        done.append(stem)

    # --- patches -------------------------------------------------------------
    def fake_savefig(self, fname, *args, **kwargs):
        handle(self, fname)
    def fake_plt_savefig(fname, *args, **kwargs):
        handle(plt.gcf(), fname)
    Figure.savefig = fake_savefig
    plt.savefig = fake_plt_savefig
    figprov.save_figure = lambda fig, out_path, **kw: handle(fig, out_path)
    figprov.save_fig = lambda fig, stem, outdir=".", **kw: handle(fig, stem)
    plt.show = lambda *args, **kwargs: None
    builtins.open = _open
    try:
        import pandas as pd
        _rc = pd.read_csv
        def _read_csv(p, *args, **kw):
            _rec(p, "r"); return _rc(p, *args, **kw)
        pd.read_csv = _read_csv
        _rt = pd.read_table
        def _read_table(p, *args, **kw):
            _rec(p, "r"); return _rt(p, *args, **kw)
        pd.read_table = _read_table
    except Exception:  # noqa: BLE001
        pass
    if a.sup_mathtext or a.sci_notation:      # mathtext in Arial (the default mathtext font is DejaVu, which would mix with Arial)
        matplotlib.rcParams["mathtext.fontset"] = "custom"
        w = {"bold": ":bold", "italic": ":italic"}.get(a.sup_mathtext or "regular", "")
        matplotlib.rcParams["mathtext.rm"] = "Arial" + w; matplotlib.rcParams["mathtext.it"] = "Arial:italic" + (w if w == ":bold" else ""); matplotlib.rcParams["mathtext.bf"] = "Arial:bold" + (w if w == ":italic" else "")
        matplotlib.rcParams["mathtext.default"] = "regular"
        if a.sup_raise is not None:
            import matplotlib._mathtext as _mt0
            _mt0.FontConstantsBase.sup1 = a.sup_raise           # Arial has no entry in _font_constant_mapping, so the base constants apply
    matplotlib.rcParams["font.family"] = "Arial"
    matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]
    matplotlib.rcParams["pdf.fonttype"] = 42

    os.chdir(cwd)
    sys.argv = [script]
    sys.path.insert(0, cwd)
    try:
        runpy.run_path(script, run_name="__main__")
    except SystemExit:
        pass
    except Exception:
        traceback.print_exc()
    if not done:
        print("regen_runner: no figure matched --select %r (script wrote nothing?)" % a.select); sys.exit(2)
    for d in done:
        print("regen_runner ->", d + ".pdf")


if __name__ == "__main__":
    main()
