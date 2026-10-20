#!/usr/bin/env python3
"""Direction of the significant DHX8 KD SE events and its binomial test (written for the code
release, 2026-10-20).

rMATS was run with --b1 control (siCont) and --b2 knockdown (siDHX8) (_v2_rmats_dhx8.sh), so
IncLevelDifference = PSI(control) - PSI(KD). An SE event is in the skipping direction when
IncLevelDifference > 0 (control PSI > KD PSI). Among SE events with FDR < 0.05 (JCEC, no
|IncLevelDifference| cutoff, as in Fig. 4B), the number in the skipping direction is tested against
p = 0.5 with an exact binomial test. The P value lies far below the smallest double (~1e-308), so
it is computed as a log10 sum of binomial probabilities instead of with scipy.stats.binomtest,
which returns 0.

Run on the server (read-only, streamed):
    ssh ncnp-linux 'python3 -' < 02_splicing/rmats/se_direction_binomial.py
"""
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import binom

D = "/mnt/HDD2/unified_results/rmats/dhx8_v2"
df = pd.read_csv(f"{D}/SE.MATS.JCEC.txt", sep="\t")
sig = df[df["FDR"] < 0.05]
n = len(sig)
k = int((sig["IncLevelDifference"] > 0).sum())
k0 = int((sig["IncLevelDifference"] == 0).sum())
print(f"significant SE events (FDR < 0.05): n = {n}")
print(f"skipping direction (control PSI > KD PSI): k = {k} ({100 * k / n:.1f}%); inclusion direction: {n - k - k0}; zero difference: {k0}")

# exact binomial, p = 0.5, in log space
logpmf = binom.logpmf(np.arange(n + 1), n, 0.5)
log_upper = logsumexp(logpmf[k:])            # P(X >= k)
log_two = np.log(2) + log_upper if k > n / 2 else None   # symmetric null: two-sided = 2 x one tail
print(f"one-sided P(X >= k) = 10^{log_upper / np.log(10):.1f}")
print(f"two-sided P         = 10^{log_two / np.log(10):.1f}")

# cross-check with the rMATS summary (Sample1 = control higher inclusion)
s = pd.read_csv(f"{D}/summary.txt", sep="\t").set_index("EventType")
print("summary.txt SE: SignificantEventsJCEC = {}, Sample1HigherInclusion = {}".format(
    s.loc["SE", "SignificantEventsJCEC"], s.loc["SE", "SigEventsJCECSample1HigherInclusion"]))
