#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
figprov.py — figure-provenance save helper for DHX8-Sakai (and any project).

Drop-in replacement for ad-hoc ``fig.savefig(...)`` calls. On every call it:

  (a) writes the figure as a PDF + PNG twin (project convention, Arial-safe);
  (b) embeds provenance in the PDF/PNG file metadata
      (matplotlib ``savefig(metadata=...)``);
  (c) writes a sidecar ``<out>.prov.json`` recording output/script/input
      sha256s, git commit (if any), python + key package versions, a copied
      script-snapshot path, an ISO timestamp (JST + UTC), hostname, cwd, argv,
      and a free-text note;
  (d) snapshots the generating script into ``<figdir>/_provenance/`` so the
      exact code survives even if the working copy is deleted;
  (e) FAILS SAFE: provenance capture never raises into your plot. Any error is
      turned into a warning and the figure is still written.

Design notes
------------
* Pure standard library + matplotlib. ``pandas``/``numpy``/... versions are
  probed opportunistically (only if importable) and never required.
* Headless-safe: uses ``datetime.now().astimezone()`` (no wall-clock ambiguity),
  no display, no network. Behaves identically on macOS (local) and ncnp-linux.
* Deterministic layout: sidecars + script snapshots live *next to the figure*
  so ``mscript-fig-update --src <fig.pdf>`` naturally carries them along and a
  reviewer can trace any panel back to its code + inputs.
* Not a git repo? ``git.commit`` is null with a reason; nothing breaks. If the
  tree is later converted to git, the same call starts recording the real
  commit automatically.

Primary API (matches project spec)
-----------------------------------

    import sys, os
    sys.path.insert(0, os.path.expanduser("~/bin"))  # if ~/bin not on PYTHONPATH
    import matplotlib
    matplotlib.use("Agg")                            # headless (server + CI)
    import matplotlib.pyplot as plt
    from figprov import save_figure

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])

    save_figure(
        fig,
        out_path="figures/fig6e_sg_overlap_summary.pdf",
        inputs=["data/limma_ir_results_dhx8.tsv",
                "data/delta_psi_results.tsv"],
        script=__file__,
        note="7-row SG overlap, unified p/IR/DEG format",
    )

Writes ``fig6e_sg_overlap_summary.pdf`` + ``.png`` (both with embedded
provenance), ``fig6e_sg_overlap_summary.pdf.prov.json``, and a script snapshot
under ``figures/_provenance/``.

``~/bin`` is not on ``sys.path`` by default. Either insert it as above, or
export ``PYTHONPATH=$HOME/bin`` in your shell / job submission script.

A legacy name-based helper ``save_fig(fig, name, outdir=..., note=...)`` is kept
for backward compatibility with earlier callers; new code should use
``save_figure`` with an explicit ``out_path`` and ``inputs=``.
"""

# NOTE: no `from __future__ import annotations` — figprov must import on the
# older Python that server analysis envs use (e.g. conda clip3 = Python 3.6),
# and that future import requires 3.7+. Keep this file free of 3.7+-only syntax.

import datetime as _dt
import hashlib as _hashlib
import inspect as _inspect
import json as _json
import os as _os
import platform as _platform
import shutil as _shutil
import socket as _socket
import subprocess as _subprocess
import sys as _sys
import warnings as _warnings

__all__ = ["save_figure", "provenance_dict", "apply_arial", "save_fig",
           "check_layout", "label_points"]

# Chunk size for streaming sha256 of large matrices (e.g. 42M intron_counts.txt).
_HASH_CHUNK = 1 << 20  # 1 MiB

# Packages worth pinning in provenance for a plotting/analysis pipeline.
# Probed opportunistically; missing ones are simply omitted.
_PKGS_OF_INTEREST = (
    "matplotlib",
    "numpy",
    "pandas",
    "scipy",
    "seaborn",
    "statsmodels",
    "sklearn",
    "matplotlib_venn",
    "upsetplot",
)

# Arial fallback chain per global CLAUDE.md "Figures & Plots".
_ARIAL_STACK = ["Arial", "Helvetica", "Liberation Sans",
                "Nimbus Sans", "DejaVu Sans"]
_ARIAL_APPLIED = False


def _warn(msg):
    """Emit a non-fatal provenance warning (never raises)."""
    _warnings.warn("[figprov] " + str(msg), RuntimeWarning, stacklevel=2)


def apply_arial():
    """Set matplotlib rcParams to the house Arial standard. Idempotent.

    Optional convenience: ``save_figure`` does NOT force this (so it never
    fights a script's own rcParams), but call it once at the top of a plotting
    script if you want the project default. Symbols must be literal Unicode in
    labels (no mathtext) per the house rule.
    """
    global _ARIAL_APPLIED
    try:
        import matplotlib
        matplotlib.rcParams.update({
            "font.family": "sans-serif",
            "font.sans-serif": _ARIAL_STACK,
            "pdf.fonttype": 42,   # editable TrueType text in PDF (Illustrator)
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": False,
        })
        _ARIAL_APPLIED = True
    except Exception as exc:  # noqa: BLE001
        _warn("apply_arial failed: %s" % exc)


def _sha256_file(path):
    """Streaming sha256 of a file. Returns hex str, or None on any failure."""
    try:
        h = _hashlib.sha256()
        with open(path, "rb") as fh:
            for block in iter(lambda: fh.read(_HASH_CHUNK), b""):
                h.update(block)
        return h.hexdigest()
    except Exception as exc:  # noqa: BLE001 - fail safe
        _warn("sha256 failed for %r: %s" % (path, exc))
        return None


def _abspath(path):
    try:
        return _os.path.abspath(_os.path.expanduser(path))
    except Exception:  # noqa: BLE001
        return path


def _git_commit(start_dir):
    """Return (commit_hash, dirty_bool) or (None, reason_str).

    Relies on ``git -C`` discovery. If the tree is not a git repo (the current
    DHX8-Sakai case) returns (None, 'not-a-git-repo').
    """
    try:
        top = _subprocess.run(
            ["git", "-C", start_dir, "rev-parse", "--show-toplevel"],
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE,
            universal_newlines=True, timeout=10,
        )
        if top.returncode != 0:
            return None, "not-a-git-repo"
        commit = _subprocess.run(
            ["git", "-C", start_dir, "rev-parse", "HEAD"],
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE,
            universal_newlines=True, timeout=10,
        )
        if commit.returncode != 0:
            return None, "no-HEAD"
        sha = commit.stdout.strip()
        status = _subprocess.run(
            ["git", "-C", start_dir, "status", "--porcelain"],
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE,
            universal_newlines=True, timeout=10,
        )
        # figprov 自身の出力（sidecar と _provenance/ のスナップショット）は
        # dirty の判定から外す。1 スクリプトが複数の図を続けて書くと、
        # 2 枚目以降が「1 枚目の sidecar が未コミット」という理由で
        # 常に dirty になっていた（2026-09-03）。
        lines = [ln for ln in status.stdout.splitlines() if ln.strip()]
        own = [ln for ln in lines
               if ln[3:].rstrip("/").endswith(".prov.json") or "/_provenance/" in ln[3:] or ln[3:].startswith("_provenance/")]
        return sha, bool(len(lines) - len(own))
    except FileNotFoundError:
        return None, "git-not-installed"
    except Exception as exc:  # noqa: BLE001
        return None, "git-error:%s" % exc


def _conda_env():
    """Best-effort conda/venv identity. None if plain system python."""
    env = _os.environ.get("CONDA_DEFAULT_ENV")
    prefix = _os.environ.get("CONDA_PREFIX") or _os.environ.get("VIRTUAL_ENV")
    if not env and prefix:
        env = _os.path.basename(prefix)
    return {"name": env, "prefix": prefix} if (env or prefix) else None


def _package_versions():
    """Version strings for packages of interest that are importable."""
    versions = {}
    for name in _PKGS_OF_INTEREST:
        try:
            mod = _sys.modules.get(name)
            if mod is not None:
                versions[name] = getattr(mod, "__version__", "unknown")
                continue
            try:
                from importlib import metadata as _md  # py3.8+
                dist = {
                    "sklearn": "scikit-learn",
                    "matplotlib_venn": "matplotlib-venn",
                }.get(name, name)
                versions[name] = _md.version(dist)
            except Exception:  # noqa: BLE001 - package absent
                continue
        except Exception:  # noqa: BLE001
            continue
    return versions


def _iso_now():
    """Local time with UTC offset (JST on these machines), plus UTC twin."""
    now = _dt.datetime.now().astimezone()
    return {
        "local": now.isoformat(timespec="seconds"),
        "utc": now.astimezone(_dt.timezone.utc).isoformat(timespec="seconds"),
        "tz": str(now.tzinfo),
    }


def _caller_script():
    """Best-effort path to the .py file that called us, if ``script`` omitted."""
    try:
        for frame in _inspect.stack():
            fn = frame.filename
            if (fn and fn != __file__ and fn.endswith(".py")
                    and _os.path.exists(fn)):
                return _abspath(fn)
    except Exception:  # noqa: BLE001
        pass
    if _sys.argv and _sys.argv[0]:
        return _abspath(_sys.argv[0])
    return None


def _snapshot_script(script_path, fig_dir):
    """Copy the generating script into <fig_dir>/_provenance/.

    Returns (snapshot_abspath, snapshot_sha256) or (None, None). The snapshot
    name embeds a 12-char content-hash tag so re-running an edited script keeps
    every distinct version rather than clobbering the previous snapshot.
    """
    try:
        if not script_path or not _os.path.isfile(script_path):
            return None, None
        src_sha = _sha256_file(script_path)
        prov_dir = _os.path.join(fig_dir, "_provenance")
        _os.makedirs(prov_dir, exist_ok=True)
        stem, ext = _os.path.splitext(_os.path.basename(script_path))
        tag = (src_sha or "nohash")[:12]
        snap = _os.path.join(prov_dir, "%s.%s%s" % (stem, tag, ext or ".py"))
        if not _os.path.exists(snap):
            _shutil.copy2(script_path, snap)
        return _abspath(snap), src_sha
    except Exception as exc:  # noqa: BLE001
        _warn("script snapshot failed: %s" % exc)
        return None, None


def provenance_dict(out_pdf, inputs, script, note, extra=None):
    """Build the provenance dictionary (no figure writes). Never raises."""
    try:
        script_abs = _abspath(script) if script else _caller_script()
        git_dir = _os.path.dirname(script_abs) if script_abs else _os.getcwd()
        commit, dirty = _git_commit(git_dir)

        input_records = []
        for p in (inputs or []):
            ap = _abspath(p)
            rec = {"path": ap, "exists": _os.path.isfile(ap)}
            if rec["exists"]:
                rec["sha256"] = _sha256_file(ap)
                try:
                    rec["bytes"] = _os.path.getsize(ap)
                except Exception:  # noqa: BLE001
                    rec["bytes"] = None
            else:
                rec["sha256"] = None
                _warn("input not found (recorded as missing): %r" % ap)
            input_records.append(rec)

        prov = {
            "schema": "figprov/1",
            "generator": "figprov.py",
            "timestamp": _iso_now(),
            "host": _socket.gethostname(),
            "platform": _platform.platform(),
            "sys_platform": _sys.platform,
            "user": _os.environ.get("USER") or _os.environ.get("LOGNAME"),
            "cwd": _os.getcwd(),
            "argv": list(_sys.argv),
            "note": note or "",
            "output": {
                "pdf": _abspath(out_pdf),
                "png": _abspath(_os.path.splitext(_abspath(out_pdf))[0] + ".png"),
                "pdf_sha256": None,   # filled by save_figure after writing
                "png_sha256": None,
            },
            "script": {
                "path": script_abs,
                "sha256": _sha256_file(script_abs) if script_abs else None,
                "snapshot": None,     # filled by save_figure after snapshot
            },
            "inputs": input_records,
            "git": {"commit": commit,
                    "dirty": (dirty if commit else None),
                    "note": (None if commit else dirty)},
            "python": {
                "version": _platform.python_version(),
                "executable": _sys.executable,
                "implementation": _sys.implementation.name,
            },
            "conda_env": _conda_env(),
            "packages": _package_versions(),
        }
        if extra:
            prov["extra"] = extra
        return prov
    except Exception as exc:  # noqa: BLE001 - fully fail-safe
        _warn("provenance_dict failed, returning stub: %s" % exc)
        return {"schema": "figprov/1", "error": str(exc),
                "timestamp": _iso_now(), "note": note or ""}


def _pdf_metadata(prov):
    """Map provenance into matplotlib PDF InfoDict keys (all values str)."""
    inputs = prov.get("inputs", [])
    inp_summary = ";".join(
        "%s@%s" % (_os.path.basename(i.get("path", "?")),
                   (i.get("sha256") or "NA")[:12])
        for i in inputs
    )
    ts = prov.get("timestamp", {}).get("local", "")
    return {
        "Title": _os.path.basename(prov.get("output", {}).get("pdf", "")),
        "Author": prov.get("user") or "",
        "Subject": "figprov provenance; note=%s" % (prov.get("note") or ""),
        "Keywords": "figprov;host=%s;script=%s;inputs=%s" % (
            prov.get("host", ""),
            _os.path.basename((prov.get("script") or {}).get("path") or ""),
            inp_summary,
        ),
        "Creator": "figprov.py (%s)" % ts,
        "Producer": "matplotlib+figprov schema=%s" % prov.get("schema", ""),
    }


def _png_metadata(prov):
    """PNG text chunks (arbitrary keys allowed by the matplotlib PNG backend)."""
    inputs = prov.get("inputs", [])
    return {
        "Software": "figprov.py",
        "Title": _os.path.basename(prov.get("output", {}).get("png", "")),
        "Author": prov.get("user") or "",
        "Description": "figprov provenance; note=%s" % (prov.get("note") or ""),
        "Source-Host": prov.get("host", ""),
        "Source-Script": _os.path.basename(
            (prov.get("script") or {}).get("path") or ""),
        "Source-Inputs": ";".join(_os.path.basename(i.get("path", "?"))
                                  for i in inputs),
        "Creation-Time": prov.get("timestamp", {}).get("local", ""),
    }


def label_points(ax, xs, ys, labels, fontsize=6, color="#555", zorder=4,
                 pad_w=22.0, pad_h=11.0, avoid_points=True,
                 point_pad_w=11.0, point_pad_h=7.0, **kwargs):
    """Annotate points, moving each label until it clears the ones already placed.

    A fixed ``xytext=(3, 3)`` offset silently stacks labels on top of each
    other wherever two cases share a coordinate, which is exactly where a
    reader most wants to tell them apart. This tries a ladder of offsets and
    takes the first that does not collide.

    Call it **after** ``set_xlim``/``set_ylim``: the collision test runs in
    display space, so it needs the final axis ranges to be in place. If a
    later ``tight_layout`` shrinks the axes the labels move closer together
    again, so the default padding carries some slack; widen ``pad_w`` /
    ``pad_h`` if a dense figure still reports overlaps.

    ``avoid_points`` also keeps a label off the plotted markers themselves,
    not only off the other labels. Without it a label clears its neighbours
    and then lands squarely on somebody else's point, which is what a network
    or a dense scatter does most of the time. Widen ``point_pad_w`` /
    ``point_pad_h`` for large markers.
    """
    cand = [(4, 4), (4, -9), (-18, 4), (-18, -9),
            (4, 12), (4, -17), (-18, 12), (-18, -17),
            (4, 20), (4, -25), (-18, 20), (-18, -25),
            (14, 8), (14, -13), (-28, 8), (-28, -13),
            (4, 28), (4, -33), (-18, 28), (-18, -33)]
    try:
        ax.get_figure().canvas.draw()
    except Exception:  # noqa: BLE001
        pass
    pts = []
    if avoid_points:
        for xv, yv in zip(xs, ys):
            try:
                pts.append(ax.transData.transform((xv, yv)))
            except Exception:  # noqa: BLE001
                pass
    placed = []
    for xv, yv, lab in zip(xs, ys, labels):
        try:
            px, py = ax.transData.transform((xv, yv))
        except Exception:  # noqa: BLE001
            px = py = 0.0
        best = cand[-1]
        for dx, dy in cand:
            bx, by = px + dx, py + dy
            if not all(abs(bx - qx) > pad_w or abs(by - qy) > pad_h
                       for qx, qy in placed):
                continue
            if not all(abs(bx - sx) > point_pad_w or abs(by - sy) > point_pad_h
                       for sx, sy in pts):
                continue
            best = (dx, dy)
            break
        placed.append((px + best[0], py + best[1]))
        ax.annotate(str(lab), (xv, yv), fontsize=fontsize, color=color,
                    xytext=best, textcoords="offset points", zorder=zorder,
                    **kwargs)


# --- Layout sanity check (overlaps / tofu glyphs / spill-out). --------------
# Rationale: a figure can pass "the script ran" and still be unreadable.
# In practice the two recurring failures are (a) legends and annotations drawn
# on top of the data, and (b) glyphs the chosen font does not have, which
# render as tofu boxes. Neither shows up unless somebody looks at the picture.
# Since save_figure is the single house entry point for writing figures, the
# check lives here so every project gets it for free. Advisory only: it warns,
# it never blocks the save.

_CHECK_MIN_OVERLAP = 0.35     # fraction of the smaller text box that must be
                              # covered before we call it an overlap
_CHECK_SPILL_PT = 3.0         # points a text may stick out of its axes


def _bbox_overlap_frac(a, b):
    """Intersection area as a fraction of the smaller box (0 if disjoint)."""
    x0, x1 = max(a.x0, b.x0), min(a.x1, b.x1)
    y0, y1 = max(a.y0, b.y0), min(a.y1, b.y1)
    if x1 <= x0 or y1 <= y0:
        return 0.0
    inter = (x1 - x0) * (y1 - y0)
    small = min(max(a.width * a.height, 1e-9), max(b.width * b.height, 1e-9))
    return inter / small


def _visible_texts(fig):
    """(artist, label) for every Text that actually draws something."""
    import matplotlib.text as _mtext
    out = []
    for ax in fig.axes:
        groups = [
            (ax.texts, "text"),
            ([ax.title, ax._left_title, ax._right_title], "title"),
            ([ax.xaxis.label, ax.yaxis.label], "axis-label"),
            (ax.get_xticklabels(), "xtick"),
            (ax.get_yticklabels(), "ytick"),
        ]
        leg = ax.get_legend()
        if leg is not None:
            groups.append((leg.get_texts(), "legend"))
        for artists, kind in groups:
            for t in artists:
                if isinstance(t, _mtext.Text) and t.get_visible() \
                        and (t.get_text() or "").strip():
                    out.append((t, kind, ax))
    for t in fig.texts:
        if t.get_visible() and (t.get_text() or "").strip():
            out.append((t, "figure-text", None))
    for leg in getattr(fig, "legends", []):
        for t in leg.get_texts():
            if t.get_visible() and (t.get_text() or "").strip():
                out.append((t, "figure-legend", None))
    return out


def _text_bbox(t, renderer):
    """Display bbox of the *glyphs* only.

    ``Annotation.get_window_extent`` unions the text box with the arrow, so an
    annotation with a long leader line reports a huge box and collides with
    everything under the arrow. Fall back to the plain Text implementation.
    """
    import matplotlib.text as _mtext
    try:
        if isinstance(t, _mtext.Annotation):
            return _mtext.Text.get_window_extent(t, renderer=renderer)
        return t.get_window_extent(renderer=renderer)
    except Exception:  # noqa: BLE001
        return None


def _short(txt, n=34):
    txt = " ".join((txt or "").split())
    return txt if len(txt) <= n else txt[: n - 1] + "…"


def _check_tofu(fig):
    """Chars the resolved font lacks — these render as tofu boxes."""
    problems = []
    try:
        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            fig.canvas.draw()
        missing = set()
        for w in caught:
            m = str(w.message)
            if "missing from font" in m:
                import re as _re
                hit = _re.search(r"\\N\{([^}]*)\}", m)
                cp = _re.search(r"Glyph (\d+)", m)
                ch = ""
                if cp:
                    try:
                        ch = chr(int(cp.group(1)))
                    except Exception:  # noqa: BLE001
                        ch = ""
                missing.add((ch, hit.group(1) if hit else m[:60]))
        for ch, name in sorted(missing):
            problems.append("tofu: 文字 %r (%s) がフォントに無く豆腐になります"
                            % (ch, name))
    except Exception as exc:  # noqa: BLE001
        _warn("tofu check skipped: %s" % exc)
    return problems


def _check_text_overlap(fig, renderer):
    """Text drawn on top of other text."""
    problems = []
    items = []
    for t, kind, ax in _visible_texts(fig):
        bb = _text_bbox(t, renderer)
        if bb is None or bb.width <= 0 or bb.height <= 0:
            continue
        items.append((t, kind, ax, bb))
    for i in range(len(items)):
        ti, ki, axi, bi = items[i]
        for j in range(i + 1, len(items)):
            tj, kj, axj, bj = items[j]
            # Tick-vs-tick is skipped outright: matplotlib's locator already
            # manages spacing within an axis, and across adjacent panels the
            # last tick of one and the first of the next routinely touch
            # without being a readability problem. Tick vs anything else is
            # still checked — that is how a text block in one panel colliding
            # with the neighbouring panel's labels gets caught.
            if ki in ("xtick", "ytick") and kj in ("xtick", "ytick"):
                continue
            f = _bbox_overlap_frac(bi, bj)
            if f >= _CHECK_MIN_OVERLAP:
                problems.append(
                    "overlap: %s「%s」と %s「%s」が %.0f%% 重なっています"
                    % (ki, _short(ti.get_text()), kj, _short(tj.get_text()),
                       f * 100))
    return problems


def _check_legend_over_data(fig, renderer):
    """Legend boxes sitting on top of plotted points/lines."""
    problems = []
    import numpy as _np
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg is None or not leg.get_visible():
            continue
        try:
            lb = leg.get_window_extent(renderer=renderer)
        except Exception:  # noqa: BLE001
            continue
        covered = 0
        try:
            for ln in ax.lines:
                if not ln.get_visible():
                    continue
                xy = ln.get_xydata()
                if xy is None or len(xy) == 0:
                    continue
                pts = ax.transData.transform(_np.asarray(xy, float))
                pts = pts[_np.isfinite(pts).all(axis=1)]
                covered += int(((pts[:, 0] >= lb.x0) & (pts[:, 0] <= lb.x1) &
                                (pts[:, 1] >= lb.y0) & (pts[:, 1] <= lb.y1)).sum())
            for col in ax.collections:
                if not col.get_visible():
                    continue
                off = col.get_offsets()
                if off is None or len(off) == 0:
                    continue
                pts = col.get_offset_transform().transform(_np.asarray(off, float))
                pts = pts[_np.isfinite(pts).all(axis=1)]
                covered += int(((pts[:, 0] >= lb.x0) & (pts[:, 0] <= lb.x1) &
                                (pts[:, 1] >= lb.y0) & (pts[:, 1] <= lb.y1)).sum())
        except Exception:  # noqa: BLE001
            continue
        if covered:
            problems.append(
                "legend: 凡例がデータ点 %d 個の上に載っています"
                "（軸の外に出すか、空いている象限へ）" % covered)
    return problems


def _check_spill(fig, renderer):
    """Annotations sticking out of their axes — these get clipped or collide."""
    problems = []
    try:
        fb = fig.bbox
    except Exception:  # noqa: BLE001
        return problems
    for t, kind, ax in _visible_texts(fig):
        if ax is None or kind in ("xtick", "ytick", "axis-label", "title"):
            continue
        bb = _text_bbox(t, renderer)
        if bb is None:
            continue
        try:
            ab = ax.get_window_extent(renderer=renderer)
        except Exception:  # noqa: BLE001
            continue
        # Horizontal spill only; putting a legend or caption *below* the axes
        # is the normal fix, not a fault.
        worst = max(bb.x1 - ab.x1, ab.x0 - bb.x0)
        if worst <= _CHECK_SPILL_PT or bb.y0 >= ab.y1 or bb.y1 <= ab.y0:
            continue
        # Spilling into an empty margin is harmless — bbox_inches='tight'
        # keeps it. Only complain when it lands on another panel or leaves
        # the canvas, which is when it actually gets cut or collides.
        onto = None
        for other in fig.axes:
            if other is ax:
                continue
            try:
                ob = other.get_window_extent(renderer=renderer)
            except Exception:  # noqa: BLE001
                continue
            if _bbox_overlap_frac(bb, ob) > 0.05:
                onto = "隣のパネル"
                break
        if onto is None and (bb.x0 < fb.x0 - 1 or bb.x1 > fb.x1 + 1):
            onto = "図の外"
        if onto:
            problems.append(
                "spill: %s「%s」が軸の外に %.0fpt はみ出し %s に掛かっています"
                % (kind, _short(t.get_text()), worst, onto))
    return problems


def check_layout(fig, verbose=True):
    """Look for the failures that only show up when you *look* at the figure.

    Checks, in order of how often they bite:

    1. **tofu** — characters the resolved font does not have. ``≥`` in a
       Japanese-font figure is the classic one; it renders as an empty box.
    2. **legend over data** — a legend placed with ``loc="best"`` or a hard
       coded corner that happens to sit on the points.
    3. **text overlapping text** — annotations colliding with each other,
       with tick labels, or with a title.
    4. **spill** — an annotation sticking out sideways past its axes, which
       gets clipped or lands on the neighbouring panel.

    Returns a list of human-readable strings (empty when the figure is clean).
    Advisory only — ``save_figure`` prints these but still writes the file.
    """
    problems = []
    try:
        problems += _check_tofu(fig)
        renderer = fig.canvas.get_renderer()
        problems += _check_legend_over_data(fig, renderer)
        problems += _check_text_overlap(fig, renderer)
        problems += _check_spill(fig, renderer)
    except Exception as exc:  # noqa: BLE001
        _warn("layout check failed (figure is unaffected): %s" % exc)
        return []
    seen, uniq = set(), []
    for pb in problems:
        if pb not in seen:
            seen.add(pb)
            uniq.append(pb)
    problems = uniq
    if verbose and problems:
        print("[figprov] レイアウト要確認 %d 件:" % len(problems))
        for p in problems:
            print("          - " + p)
    return problems


def save_figure(fig, out_path, inputs=None, script=None, note="",
                dpi=300, png=True, pdf=True, facecolor="white",
                extra=None, savefig_kwargs=None, check=True):
    """Save ``fig`` as PDF+PNG twins with embedded + sidecar provenance.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to write.
    out_path : str
        Output path. Extension is ignored for twin generation: both
        ``<stem>.pdf`` and ``<stem>.png`` are written (project convention).
        May be relative to cwd or absolute.
    inputs : list[str], optional
        Data files the figure was computed from. Each is sha256-hashed and
        recorded. Missing files are recorded ``exists=false`` (warn only) so
        server-vs-local path drift never crashes the plot. Skip the giant count
        matrices you did not directly plot from; pass only true inputs.
    script : str, optional
        Path to the generating script — pass ``__file__``. If omitted, figprov
        infers the caller. Snapshotted into ``<figdir>/_provenance/`` and hashed.
    note : str, optional
        Free-text note (mirrors ``mscript-fig-update --note``).
    dpi : int, optional
        Raster DPI for PNG (default 300, matching manuscript scripts).
    png, pdf : bool, optional
        Toggle each twin. Both default True.
    facecolor : str, optional
        Figure background (default 'white', matching project convention).
    extra : dict, optional
        Extra JSON-serializable metadata for the sidecar (parameters,
        thresholds, random seed, etc.).
    savefig_kwargs : dict, optional
        Extra kwargs forwarded to ``fig.savefig`` (e.g. ``bbox_inches='tight'``).
        ``metadata`` and ``facecolor`` are managed by figprov.
    check : bool, optional
        Run :func:`check_layout` before saving and report tofu glyphs, legends
        covering data, overlapping text and annotations spilling out of their
        axes. Findings are printed and recorded in the sidecar under
        ``layout_check``. Advisory only — the figure is written either way.
        Default True; pass ``check=False`` for a figure you have already
        inspected and know trips a false positive.

    Returns
    -------
    dict
        The provenance dictionary that was written (also useful for logging).

    Notes
    -----
    The figure is *always* written even if provenance capture fails. Provenance
    errors surface only as ``RuntimeWarning``.
    """
    stem, _ext = _os.path.splitext(_os.path.expanduser(out_path))
    pdf_path = stem + ".pdf"
    png_path = stem + ".png"
    fig_dir = _os.path.dirname(_abspath(pdf_path)) or _os.getcwd()
    sk = dict(savefig_kwargs or {})
    sk.pop("metadata", None)   # figprov owns metadata
    sk.pop("facecolor", None)  # managed explicitly below
    sk.setdefault("bbox_inches", "tight")

    # --- Build provenance (best-effort; never fatal). ---
    try:
        prov = provenance_dict(pdf_path, inputs, script, note, extra=extra)
        script_abs = (prov.get("script") or {}).get("path")
        snap_path, _snap_sha = _snapshot_script(script_abs, fig_dir)
        if snap_path:
            prov.setdefault("script", {})["snapshot"] = snap_path
        pdf_md = _pdf_metadata(prov)
        png_md = _png_metadata(prov)
    except Exception as exc:  # noqa: BLE001
        _warn("provenance prep failed; writing figure without metadata: %s"
              % exc)
        prov, pdf_md, png_md = ({"schema": "figprov/1", "error": str(exc)},
                                {}, {})

    # --- Look at the picture before writing it. A figure can be numerically
    #     correct and still unreadable; this is the only automatic pass that
    #     catches tofu glyphs and things drawn on top of the data. ---
    layout_problems = check_layout(fig) if check else []
    if layout_problems:
        try:
            prov["layout_check"] = {"n_problems": len(layout_problems),
                                    "problems": layout_problems}
        except Exception:  # noqa: BLE001
            pass

    # --- Write the figure twins. A failed plot is a real failure the user must
    #     see, so these propagate; but each is attempted independently. ---
    if pdf:
        try:
            _os.makedirs(_os.path.dirname(pdf_path) or ".", exist_ok=True)
            fig.savefig(pdf_path, metadata=(pdf_md or None),
                        facecolor=facecolor, **sk)
        except TypeError:
            _warn("PDF metadata= unsupported; writing PDF without embedded meta")
            fig.savefig(pdf_path, facecolor=facecolor, **sk)
    if png:
        try:
            _os.makedirs(_os.path.dirname(png_path) or ".", exist_ok=True)
            fig.savefig(png_path, dpi=dpi, metadata=(png_md or None),
                        facecolor=facecolor,
                        bbox_inches=sk.get("bbox_inches", "tight"))
        except TypeError:
            _warn("PNG metadata= unsupported; writing PNG without embedded meta")
            fig.savefig(png_path, dpi=dpi, facecolor=facecolor,
                        bbox_inches=sk.get("bbox_inches", "tight"))

    # --- Fill output hashes now that files exist, then write the sidecar. ---
    try:
        out = prov.setdefault("output", {})
        if pdf and _os.path.isfile(pdf_path):
            out["pdf_sha256"] = _sha256_file(pdf_path)
        if png and _os.path.isfile(png_path):
            out["png_sha256"] = _sha256_file(png_path)
        sidecar = pdf_path + ".prov.json"
        with open(sidecar, "w", encoding="utf-8") as fh:
            _json.dump(prov, fh, ensure_ascii=False, indent=2, sort_keys=False)
            fh.write("\n")
        print("[figprov] wrote %s (+ .png) with provenance -> %s"
              % (_os.path.basename(pdf_path), _os.path.basename(sidecar)))
    except Exception as exc:  # noqa: BLE001
        _warn("sidecar write failed (figure is still saved): %s" % exc)

    return prov


# --- Backward-compatible legacy helper (name + outdir signature). ------------
def save_fig(fig, name, outdir="figures", inputs=None, note="", dpi=300,
             **kwargs):
    """Legacy name-based wrapper for older callers.

    ``save_fig(fig, "fig1_overview", outdir="figures", note="...")`` maps onto
    ``save_figure``. New code should call ``save_figure`` with an explicit
    ``out_path`` and ``inputs=``. A trailing ``.pdf``/``.png`` on ``name`` is
    stripped.
    """
    base = name
    for ext in (".pdf", ".png"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    out_path = base if _os.path.isabs(base) else _os.path.join(outdir, base)
    script = kwargs.pop("script", None)
    return save_figure(fig, out_path + ".pdf", inputs=inputs, script=script,
                       note=note, dpi=dpi, **kwargs)


# --- CLI: verify/inspect a sidecar without importing matplotlib. -------------
def _cli(argv):
    import argparse
    ap = argparse.ArgumentParser(
        prog="figprov.py",
        description="Inspect / verify a figprov sidecar (.prov.json), or run a "
                    "self-test.")
    sub = ap.add_subparsers(dest="cmd")

    p_show = sub.add_parser("show", help="pretty-print a sidecar")
    p_show.add_argument("prov_json")

    p_ver = sub.add_parser("verify", help="re-hash inputs/outputs vs sidecar")
    p_ver.add_argument("prov_json")

    sub.add_parser("selftest", help="write a trivial figure with provenance")

    p_bf = sub.add_parser(
        "backfill",
        help="attach provenance to an EXISTING figure without re-rendering "
             "(hashes the current output + script + inputs, snapshots the script)")
    p_bf.add_argument("figure", help="existing figure path (.pdf/.png or stem)")
    p_bf.add_argument("--script", required=True, help="generating script path")
    p_bf.add_argument("--inputs", default="",
                      help="comma-separated data files the figure was computed from")
    p_bf.add_argument("--note", default="")

    args = ap.parse_args(argv)

    if args.cmd == "selftest" or args.cmd is None:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        apply_arial()
        f, a = plt.subplots(figsize=(3, 2))
        a.plot([0, 1, 2], [0, 1, 4])
        a.set_title("figprov self-test (β ± μ)")  # literal Unicode, no mathtext
        rec = save_figure(f, "/tmp/figprov_test/selftest.pdf",
                          inputs=[__file__], script=__file__,
                          note="self-test")
        print(_json.dumps(rec, indent=2, ensure_ascii=False))
        return 0

    if args.cmd == "backfill":
        stem, _ext = _os.path.splitext(_os.path.expanduser(args.figure))
        pdf_path, png_path = stem + ".pdf", stem + ".png"
        inputs = [x for x in args.inputs.split(",") if x]
        prov = provenance_dict(pdf_path, inputs, args.script, args.note)
        fig_dir = _os.path.dirname(_abspath(pdf_path)) or _os.getcwd()
        script_abs = (prov.get("script") or {}).get("path")
        snap, _snap_sha = _snapshot_script(script_abs, fig_dir)
        if snap:
            prov.setdefault("script", {})["snapshot"] = snap
        out = prov.setdefault("output", {})
        if _os.path.isfile(pdf_path):
            out["pdf_sha256"] = _sha256_file(pdf_path)
        if _os.path.isfile(png_path):
            out["png_sha256"] = _sha256_file(png_path)
        prov["backfilled"] = True   # attached post-hoc, not captured at render time
        sidecar = pdf_path + ".prov.json"
        with open(sidecar, "w", encoding="utf-8") as fh:
            _json.dump(prov, fh, ensure_ascii=False, indent=2, sort_keys=False)
            fh.write("\n")
        print("[figprov] backfilled provenance -> %s" % sidecar)
        return 0

    with open(args.prov_json, encoding="utf-8") as fh:
        prov = _json.load(fh)

    if args.cmd == "show":
        print(_json.dumps(prov, ensure_ascii=False, indent=2))
        return 0

    # verify
    problems = 0
    for rec in prov.get("inputs", []):
        p, want = rec.get("path"), rec.get("sha256")
        if not want:
            continue
        got = _sha256_file(p) if p and _os.path.isfile(p) else None
        if got != want:
            problems += 1
        print("%s input  %s" % ("OK " if got == want else "MISS", p))
    out = prov.get("output", {})
    for kind in ("pdf", "png"):
        p, want = out.get(kind), out.get(kind + "_sha256")
        if not (p and want):
            continue
        got = _sha256_file(p) if _os.path.isfile(p) else None
        if got != want:
            problems += 1
        print("%s output %s" % ("OK " if got == want else "MISS", p))
    print("VERIFY: %d mismatch(es)" % problems)
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(_cli(_sys.argv[1:]))
