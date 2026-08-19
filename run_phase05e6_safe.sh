#!/usr/bin/env bash
set +e

echo "========================================"
echo "PHASE 05E-6 — SIGNAL-TO-DRIFT AUDIT"
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
  src/analyze_phase05e6_signal_audit.py

RC=$?

echo
echo "Analysis return code: $RC"

if [ "$RC" -eq 0 ]; then

    echo
    echo "=== OUTPUT CHECK ==="

    wc -l \
      results/phase05e6_signal_audit/signal_to_drift_audit.csv \
      results/phase05e6_signal_audit/phase05e6_summary.txt

    ls -lh \
      results/phase05e6_signal_audit/signal_to_drift_audit.csv \
      results/phase05e6_signal_audit/phase05e6_summary.txt

    echo
    echo "=== FINAL SUMMARY ==="

    cat \
      results/phase05e6_signal_audit/phase05e6_summary.txt

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
