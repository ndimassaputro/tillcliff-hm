#!/usr/bin/env bash
set +e

echo "========================================"
echo "PHASE 05D-4 — DISTRIBUTED RESPONSE"
echo "========================================"

echo
echo "=== INPUT CHECK ==="

ls -lh \
  results/phase05d2_mesh_erosion/coarse/*.pvd \
  results/phase05d2_mesh_erosion/medium/*.pvd \
  results/phase05d2_mesh_erosion/fine/*.pvd

echo
echo "=== ANALYZE EXISTING OUTPUT ==="

python \
  src/analyze_phase05d4_distributed_response.py

RC=$?

echo
echo "Analysis return code: $RC"

if [ "$RC" -eq 0 ]; then

    echo
    echo "=== OUTPUT CHECK ==="

    wc -l \
      results/phase05d4_distributed_response/distributed_response_history.csv \
      results/phase05d4_distributed_response/common_window_mesh_comparison.csv \
      results/phase05d4_distributed_response/phase05d4_summary.txt

    ls -lh \
      results/phase05d4_distributed_response/distributed_response_history.csv \
      results/phase05d4_distributed_response/common_window_mesh_comparison.csv \
      results/phase05d4_distributed_response/figure_01_body_rms_vs_E.png \
      results/phase05d4_distributed_response/figure_02_body_rms_vs_removed_area.png \
      results/phase05d4_distributed_response/phase05d4_summary.txt

    echo
    echo "=== FINAL SUMMARY ==="

    cat \
      results/phase05d4_distributed_response/phase05d4_summary.txt

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
