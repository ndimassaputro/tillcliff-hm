#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

MODEL="$ROOT/model/phase05e2_mc_antecedent"
RUN="$ROOT/results/phase05e2_mc_antecedent"

echo "========================================"
echo "PHASE 05E-2 — MC ANTECEDENT BRANCHES"
echo "========================================"

echo
echo "=== BUILD ==="

python \
  src/build_phase05e2_mc_antecedent.py

BUILD_RC=$?

echo
echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== FILE CHECK ==="

    for STATE in dry reference wet
    do

        wc -l \
          "$MODEL/$STATE/mc_antecedent_hold.prj"

        ls -lh \
          "$MODEL/$STATE/bulk.vtu" \
          "$MODEL/$STATE/mc_antecedent_hold.prj"

    done

    echo
    echo "=== XML CHECK ==="

    for STATE in dry reference wet
    do

        xmllint \
          --noout \
          "$MODEL/$STATE/mc_antecedent_hold.prj"

        echo \
          "$STATE XML RC: $?"

    done

    echo
    echo "=== CLEAN OUTPUT ==="

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
          "$MODEL/$STATE/mc_antecedent_hold.prj" \
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

            tail -n 45 \
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
      src/analyze_phase05e2_mc_antecedent.py

    ANALYSIS_RC=$?

    echo
    echo \
      "Analysis return code: $ANALYSIS_RC"

    if [ "$ANALYSIS_RC" -eq 0 ]; then

        echo
        echo "=== FINAL SUMMARY ==="

        cat \
          results/phase05e2_analysis/phase05e2_summary.txt

    fi
fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
