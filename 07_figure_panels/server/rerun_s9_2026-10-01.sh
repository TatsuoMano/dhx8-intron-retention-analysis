#!/bin/bash
# v3 panels, server scripts (ncnp-linux). Working dir: /mnt/HDD2/unified_results/figregen_2026-09-29 (project dir, not /tmp).
# Same rules as regen_local_v3.sh. Outputs: panels/ (copied back to Analysis/2026-09-29_illustrator_v3/panels).
set -u
cd /mnt/HDD2/unified_results/figregen_2026-09-29
export PYTHONPATH=$HOME/bin
UR=/mnt/HDD2/unified_results
R="python3 regen_runner.py"
K="--keep-look --keep-ticks"
run() { $R "$@" 2>&1 | grep -E "regen_runner ->|Error|Traceback|no figure" ; }
run $UR/analysis_sashimi_dhx8_events.py --figsize 3.3 2.2 --cap 7 $K --text-replace "RPKM" "Depth" --text-replace "FDR = 0.0" "FDR < 10⁻¹⁶" --sup-mathtext bold --sup-min-pt 5.5 --sup-raise 0.33 --text-offset-pt "^(42|91|248|245|120|124|8|77|85|4028|2839|23|395)$" 0 2.5 --text-offset-pt "^(93|69)$" -10 1 --text-offset-pt "^573$" 0 6 --text-offset-pt "^CE \(175" 18 0 --text-offset-pt "^Ccdc71 \(chr9" -14 0 --text-offset-pt "^IR \(intron" 0 -2.5 --title-regex "\s+\(" "\n(" --out-dir panels/sashimi   # titles in two lines: at 7 pt the longest is wider than the panel
echo S9_DONE
