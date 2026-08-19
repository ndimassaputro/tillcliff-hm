#!/usr/bin/env bash
set +e

echo "========================================"
echo "PHASE 05D-3 — LOCALIZATION AUDIT"
echo "========================================"

echo
echo "=== WRITE/INPUT CHECK ==="

ls -lh \
  model/phase05d2_mesh_erosion/coarse/bulk.vtu \
  model/phase05d2_mesh_erosion/medium/bulk.vtu \
  model/phase05d2_mesh_erosion/fine/bulk.vtu

ls -lh \
  results/phase05d2_mesh_erosion/coarse/*.pvd \
  results/phase05d2_mesh_erosion/medium/*.pvd \
  results/phase05d2_mesh_erosion/fine/*.pvd

echo
echo "=== ANALYZE EXISTING OUTPUT ==="

python \
  src/analyze_phase05d3_localization.py

RC=$?

echo
echo "Analysis return code: $RC"

if [ "$RC" -eq 0 ]; then

    echo
    echo "=== OUTPUT CHECK ==="

    wc -l \
      results/phase05d3_localization/localization_history.csv \
      results/phase05d3_localization/mesh_localization_summary.csv \
      results/phase05d3_localization/phase05d3_summary.txt

    ls -lh \
      results/phase05d3_localization/localization_history.csv \
      results/phase05d3_localization/mesh_localization_summary.csv \
      results/phase05d3_localization/figure_01_hotspot_interface_distance.png \
      results/phase05d3_localization/figure_02_integrated_plastic_area.png \
      results/phase05d3_localization/phase05d3_summary.txt

    echo
    echo "=== FINAL SUMMARY ==="

    cat \
      results/phase05d3_localization/phase05d3_summary.txt

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
