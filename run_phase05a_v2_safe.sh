#!/usr/bin/env bash
set +e

echo "========================================"
echo "PHASE 05A-V2 — BASAL TOE NOTCH"
echo "========================================"

python \
  src/build_phase05a_v2_toe_notch.py

RC=$?

echo
echo "Build return code: $RC"

if [ "$RC" -eq 0 ]; then

    echo
    echo "=== OUTPUT CHECK ==="

    wc -l \
      results/phase05a_v2_toe_notch/toe_notch_bands.csv \
      results/phase05a_v2_toe_notch/toe_notch_levels.csv \
      results/phase05a_v2_toe_notch/phase05a_v2_summary.txt

    ls -lh \
      model/phase05a_v2_toe_notch/slope_toe_notch_ready.vtu \
      results/phase05a_v2_toe_notch/figure_01_basal_toe_notch.png \
      results/phase05a_v2_toe_notch/phase05a_v2_summary.txt

    echo
    echo "=== SUMMARY ==="

    cat \
      results/phase05a_v2_toe_notch/phase05a_v2_summary.txt

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
