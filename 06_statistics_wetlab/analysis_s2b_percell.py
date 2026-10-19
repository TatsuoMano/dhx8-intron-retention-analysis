#!/usr/bin/env python3
"""Supplementary Fig. S2B: per-cell statistics of the arsenite time course (GFP-G3BP1 stress granules per cell).

The Prism file stores one row per cell and one column per minute after arsenite (the same cells followed over time).
The P values printed in the original panel (0.0003 / 0.9154) came from an ordinary one-way ANOVA + Dunnett's test of
1 min vs each later time point, which treats the repeated measurements as independent groups (identified 2026-09-28).
This script computes the per-cell analysis adopted on 2026-10-01 (author decision):

  within each group   change from 1 to 20 min per cell; exact two-sided Wilcoxon signed-rank test (paired t for reference)
  between the groups  Welch's t-test on the per-cell changes (Mann-Whitney and a 10,000-fold label permutation for reference)
  for the record      the ANOVA + Dunnett reproduction of the printed values

Outputs (next to this script): s2b_percell_stats.tsv, S2B_STATS.md (legend sentence), s2b_percell_changes.tsv (the per-cell values).
"""
import os, sys, json
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "2026-09-16_figcomp_reproduction"))
from prism_panels import read_pzfx, FFF  # noqa: E402

PZ = os.path.join(FFF, "SupFig02", "SupFig02a. 2SGs定量_230609IC.pzfx")
TABLE = "SGs number"


def fmt_p(p):
    if p < 1e-4:
        m, e = ("%.1e" % p).split("e")
        return "%s × 10^%d^" % (m, int(e))
    return "%.3f" % p if p >= 0.001 else "%.4f" % p


def main():
    t = read_pzfx(PZ)[TABLE]
    xlab = [str(v) for v in t["x"]]
    groups = {}
    for name, subs in t["cols"]:
        a = np.array(subs, dtype=float)                      # cells x time points
        ok = ~np.isnan(a[:, 0]) & ~np.isnan(a[:, -1])
        groups[name.strip()] = a[ok]
    names = list(groups)
    rows, per_cell = [], []
    changes = {}
    for name in names:
        a = groups[name]
        d = a[:, -1] - a[:, 0]
        changes[name] = d
        w = stats.wilcoxon(d, alternative="two-sided", method="exact")
        pt = stats.ttest_rel(a[:, -1], a[:, 0])
        rows.append(dict(group=name, n_cells=len(d), mean_1min=a[:, 0].mean(), mean_20min=a[:, -1].mean(),
                         mean_change=d.mean(), sd_change=d.std(ddof=1), sem_change=d.std(ddof=1) / np.sqrt(len(d)),
                         n_increased=int((d > 0).sum()), n_unchanged=int((d == 0).sum()), n_decreased=int((d < 0).sum()),
                         wilcoxon_exact_p=w.pvalue, paired_t_p=pt.pvalue))
        for i, (v1, v20, dv) in enumerate(zip(a[:, 0], a[:, -1], d)):
            per_cell.append((name, i + 1, v1, v20, dv))
        # ANOVA + Dunnett across time points, each time point as an independent group (the printed values)
        cols = [a[:, j][~np.isnan(a[:, j])] for j in range(a.shape[1])]
        dn = stats.dunnett(*cols[1:], control=cols[0])
        rows[-1]["printed_dunnett_1_vs_20"] = float(dn.pvalue[-1])
        rows[-1]["anova_p_across_timepoints"] = float(stats.f_oneway(*cols).pvalue)
    dA, dB = changes[names[0]], changes[names[1]]
    welch = stats.ttest_ind(dA, dB, equal_var=False)
    mw = stats.mannwhitneyu(dA, dB, alternative="two-sided")
    rng = np.random.default_rng(0)
    pooled = np.concatenate([dA, dB]); obs = dA.mean() - dB.mean(); cnt = 0; N = 10000
    for _ in range(N):
        p = rng.permutation(pooled)
        cnt += abs(p[:len(dA)].mean() - p[len(dA):].mean()) >= abs(obs) - 1e-12
    perm_p = (cnt + 1) / (N + 1)
    # Holm over the two within-group tests
    ps = sorted((r["wilcoxon_exact_p"], r["group"]) for r in rows)
    holm = {ps[0][1]: min(1.0, 2 * ps[0][0]), ps[1][1]: min(1.0, max(ps[1][0], 2 * ps[0][0]))}
    for r in rows:
        r["wilcoxon_holm_p"] = holm[r["group"]]
    between = dict(welch_t_p=welch.pvalue, welch_t=welch.statistic, mann_whitney_p=mw.pvalue, permutation_p=perm_p,
                   n_permutations=N, mean_difference=obs)

    with open(os.path.join(HERE, "s2b_percell_stats.tsv"), "w") as f:
        keys = list(rows[0].keys())
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(("%.6g" % r[k]) if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
        f.write("\n# between groups (per-cell changes 1 -> 20 min, %s vs %s)\n" % (names[0], names[1]))
        for k, v in between.items():
            f.write("%s\t%s\n" % (k, ("%.6g" % v) if isinstance(v, float) else v))
        f.write("\n# time labels in the Prism file: %s\n" % ", ".join(xlab))
    with open(os.path.join(HERE, "s2b_percell_changes.tsv"), "w") as f:
        f.write("group\tcell\tsg_1min\tsg_20min\tchange\n")
        for g, i, v1, v20, dv in per_cell:
            f.write("%s\t%d\t%g\t%g\t%g\n" % (g, i, v1, v20, dv))
    A, B = rows[0], rows[1]
    legend = ("Change from 1 to 20 min: %s +%.1f stress granules per cell (%d/%d cells increased, Wilcoxon signed-rank P = %s), "
              "%s +%.1f (%d/%d cells, P = %s); difference between the groups, Welch's *t*-test on the per-cell changes, P = %s."
              % (A["group"], A["mean_change"], A["n_increased"], A["n_cells"], fmt_p(A["wilcoxon_exact_p"]),
                 B["group"], B["mean_change"], B["n_increased"], B["n_cells"], fmt_p(B["wilcoxon_exact_p"]), fmt_p(between["welch_t_p"])))
    md = ["# S2B per-cell statistics (%s)" % os.path.basename(PZ), "",
          "Time labels in the Prism file: %s" % ", ".join(xlab), "",
          "| group | n | mean 1 min | mean 20 min | mean change (SD) | increased | Wilcoxon exact P | Holm | paired t P | printed Dunnett 1 vs 20 (reproduced) |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append("| %s | %d | %.1f | %.1f | %+.1f (%.1f) | %d/%d | %.2e | %.2e | %.3f | %.4f |" % (
            r["group"], r["n_cells"], r["mean_1min"], r["mean_20min"], r["mean_change"], r["sd_change"], r["n_increased"], r["n_cells"],
            r["wilcoxon_exact_p"], r["wilcoxon_holm_p"], r["paired_t_p"], r["printed_dunnett_1_vs_20"]))
    md += ["", "Between groups (per-cell changes): Welch t = %.2f, P = %.3f; Mann-Whitney P = %.3f; permutation P = %.3f (%d permutations, seed 0); mean difference %.1f."
           % (between["welch_t"], between["welch_t_p"], between["mann_whitney_p"], between["permutation_p"], N, obs), "",
           "Legend sentence:", "", legend, ""]
    open(os.path.join(HERE, "S2B_STATS.md"), "w", encoding="utf-8").write("\n".join(md))
    json.dump({"rows": rows, "between": between, "x_labels": xlab, "legend": legend}, open(os.path.join(HERE, "s2b_percell_stats.json"), "w"), indent=1, default=float)
    print("\n".join(md))


if __name__ == "__main__":
    main()
