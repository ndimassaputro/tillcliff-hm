#!/usr/bin/env bash
set +e

echo "========================================"
echo "PHASE 05E-5 — DEFORMATION-MODE AUDIT"
echo "========================================"

echo
echo "=== INPUT CHECK ==="

for CASE in \
  refP_refS \
  dryP_refS \
  wetP_refS
do

    ls -lh \
      results/phase05e4_decomposition/"$CASE"/*.pvd

done

echo
echo "=== RUN EXISTING-OUTPUT ANALYSIS ==="

python \
  src/analyze_phase05e5_mode_audit.py

RC=$?

echo
echo "Analysis return code: $RC"

if [ "$RC" -eq 0 ]; then

    echo
    echo "=== OUTPUT CHECK ==="

    wc -l \
      results/phase05e5_mode_audit/mode_shape_audit.csv \
      results/phase05e5_mode_audit/phase05e5_summary.txt

    ls -lh \
      results/phase05e5_mode_audit/mode_shape_audit.csv \
      results/phase05e5_mode_audit/phase05e5_summary.txt

    echo
    echo "=== FINAL SUMMARY ==="

    cat \
      results/phase05e5_mode_audit/phase05e5_summary.txt

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
