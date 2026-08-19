#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

echo "========================================"
echo "PHASE 05A — EROSION-READY TOE MESH"
echo "========================================"

echo
echo "=== BUILD + CHECK ==="

python \
  src/build_phase05a_erosion_mesh.py

RC=$?

echo
echo "Build return code: $RC"

if [ "$RC" -eq 0 ]; then

    echo
    echo "=== OUTPUT CHECK ==="

    wc -l \
      results/phase05a_erosion_mesh/toe_erosion_bands.csv \
      results/phase05a_erosion_mesh/toe_erosion_levels.csv \
      results/phase05a_erosion_mesh/phase05a_summary.txt

    ls -lh \
      model/phase05a_erosion_mesh/slope_erosion_ready.vtu \
      results/phase05a_erosion_mesh/toe_erosion_bands.csv \
      results/phase05a_erosion_mesh/toe_erosion_levels.csv \
      results/phase05a_erosion_mesh/figure_01_toe_erosion_bands.png \
      results/phase05a_erosion_mesh/phase05a_summary.txt

    echo
    echo "=== SUMMARY ==="

    cat \
      results/phase05a_erosion_mesh/phase05a_summary.txt

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
