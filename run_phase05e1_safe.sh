#!/usr/bin/env bash
set +e

echo "========================================"
echo "PHASE 05E-1 — ANTECEDENT STATE AUDIT"
echo "========================================"

echo
echo "=== WRITE CHECK ==="

wc -l \
  src/audit_phase05e1_antecedent_states.py

ls -lh \
  src/audit_phase05e1_antecedent_states.py

echo
echo "=== RUN ==="

python \
  src/audit_phase05e1_antecedent_states.py

RC=$?

echo
echo "Audit return code: $RC"

if [ "$RC" -eq 0 ]; then

    echo
    echo "=== OUTPUT CHECK ==="

    wc -l \
      results/phase05e1_antecedent_audit/antecedent_state_audit.csv \
      results/phase05e1_antecedent_audit/phase05e1_summary.txt

    ls -lh \
      results/phase05e1_antecedent_audit/antecedent_state_audit.csv \
      results/phase05e1_antecedent_audit/phase05e1_summary.txt

    echo
    echo "=== FINAL SUMMARY ==="

    cat \
      results/phase05e1_antecedent_audit/phase05e1_summary.txt

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
