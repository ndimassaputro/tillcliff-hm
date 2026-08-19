#!/usr/bin/env bash
set +e

echo "========================================"
echo "PHASE 05D-0 — MESH REFINEMENT PREFLIGHT"
echo "========================================"

echo
echo "=== RUN AUDIT ==="

python \
  src/audit_phase05d_mesh_refinement_preflight.py

RC=$?

echo
echo "Audit return code: $RC"

if [ "$RC" -eq 0 ]; then

    echo
    echo "=== OUTPUT CHECK ==="

    wc -l \
      results/phase05d_preflight/erosion_event_audit.csv \
      results/phase05d_preflight/tool_availability.csv \
      results/phase05d_preflight/phase05d_preflight_summary.txt

    ls -lh \
      results/phase05d_preflight/erosion_event_audit.csv \
      results/phase05d_preflight/tool_availability.csv \
      results/phase05d_preflight/phase05d_preflight_summary.txt

    echo
    echo "=== FINAL SUMMARY ==="

    cat \
      results/phase05d_preflight/phase05d_preflight_summary.txt

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
