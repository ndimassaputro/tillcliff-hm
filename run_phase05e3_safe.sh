#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

MODEL="$ROOT/model/phase05e3_antecedent_erosion"
RUN="$ROOT/results/phase05e3_antecedent_erosion"

echo "========================================"
echo "PHASE 05E-3 — ANTECEDENT × EROSION"
echo "========================================"

echo
echo "=== BUILD ==="

python \
  src/build_phase05e3_antecedent_erosion.py

BUILD_RC=$?

echo
echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== WRITE CHECK ==="

    for STATE in dry reference wet
    do

        wc -l \
          "$MODEL/$STATE/antecedent_erosion.prj"

        ls -lh \
          "$MODEL/$STATE/bulk.vtu" \
          "$MODEL/$STATE/antecedent_erosion.prj"

    done

    echo
    echo "=== XML CHECK ==="

    for STATE in dry reference wet
    do

        xmllint \
          --noout \
          "$MODEL/$STATE/antecedent_erosion.prj"

        echo \
          "$STATE XML RC: $?"

    done

    echo
    echo "=== CLEAN RUN OUTPUT ==="

    rm -rf "$RUN"
    mkdir -p "$RUN"

    cat > "$RUN/run_status.csv" <<'CSV'
state,return_code
CSV

    echo
    echo "=== RUN DRY / REFERENCE / WET ==="

    for STATE in dry reference wet
    do

        echo
        echo "========================================"
        echo "RUN $STATE"
        echo "========================================"

        STATE_OUT="$RUN/$STATE"

        mkdir -p \
          "$STATE_OUT"

        ogs \
          "$MODEL/$STATE/antecedent_erosion.prj" \
          -m "$MODEL/$STATE" \
          -o "$STATE_OUT" \
          2>&1 \
          | tee \
          "$STATE_OUT/ogs.log"

        OGS_RC=${PIPESTATUS[0]}

        echo
        echo \
          "$STATE OGS RETURN CODE: $OGS_RC"

        echo \
          "$STATE,$OGS_RC" \
          >> "$RUN/run_status.csv"

        if [ "$OGS_RC" -ne 0 ]; then

            echo
            echo \
              "=== $STATE FAILURE TAIL ==="

            tail -n 40 \
              "$STATE_OUT/ogs.log"

        fi

    done

    echo
    echo "=== RUN STATUS ==="

    cat \
      "$RUN/run_status.csv"

    echo
    echo "=== ANALYZE ==="

    python \
      src/analyze_phase05e3_antecedent_erosion.py

    RC=$?

    echo
    echo \
      "Analysis return code: $RC"

    if [ "$RC" -eq 0 ]; then

        echo
        echo "=== FINAL STATE COMPARISON ==="

        cat \
          results/phase05e3_analysis/antecedent_state_comparison.csv

        echo
        echo "=== FINAL SUMMARY ==="

        cat \
          results/phase05e3_analysis/phase05e3_summary.txt

    fi
fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
