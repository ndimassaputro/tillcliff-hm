#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

MODEL="$ROOT/model/phase05d2_mesh_erosion"
RUN="$ROOT/results/phase05d2_mesh_erosion"

echo "========================================"
echo "PHASE 05D-2 — MESH EROSION ROBUSTNESS"
echo "========================================"

echo
echo "=== BUILD ==="

python \
  src/build_phase05d2_mesh_erosion.py

BUILD_RC=$?

echo
echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== WRITE CHECK ==="

    wc -l \
      "$MODEL/coarse/erosion_probe.prj" \
      "$MODEL/medium/erosion_probe.prj" \
      "$MODEL/fine/erosion_probe.prj" \
      "$MODEL/mesh_cases.csv"

    ls -lh \
      "$MODEL/coarse/erosion_probe.prj" \
      "$MODEL/medium/erosion_probe.prj" \
      "$MODEL/fine/erosion_probe.prj" \
      "$MODEL/mesh_cases.csv"

    echo
    echo "=== XML CHECK ==="

    for CASE in coarse medium fine
    do

        xmllint \
          --noout \
          "$MODEL/$CASE/erosion_probe.prj"

        echo \
          "$CASE XML RC: $?"

    done

    echo
    echo "=== CLEAN OUTPUT ==="

    rm -rf "$RUN"
    mkdir -p "$RUN"

    cat > "$RUN/run_status.csv" <<'CSV'
case,return_code
CSV

    echo
    echo "=== RUN COARSE / MEDIUM / FINE ==="

    for CASE in coarse medium fine
    do

        echo
        echo "========================================"
        echo "RUN $CASE"
        echo "========================================"

        CASE_OUT="$RUN/$CASE"

        mkdir -p \
          "$CASE_OUT"

        ogs \
          "$MODEL/$CASE/erosion_probe.prj" \
          -m "$MODEL/$CASE" \
          -o "$CASE_OUT" \
          2>&1 \
          | tee \
          "$CASE_OUT/ogs.log"

        OGS_RC=${PIPESTATUS[0]}

        echo
        echo \
          "$CASE OGS RETURN CODE: $OGS_RC"

        echo \
          "$CASE,$OGS_RC" \
          >> "$RUN/run_status.csv"

        if [ "$OGS_RC" -ne 0 ]; then

            echo
            echo \
              "=== $CASE FAILURE TAIL ==="

            tail -n 35 \
              "$CASE_OUT/ogs.log"

        fi

    done

    echo
    echo "=== RUN STATUS ==="

    cat \
      "$RUN/run_status.csv"

    echo
    echo "=== ANALYZE ==="

    python \
      src/analyze_phase05d2_mesh_erosion.py

    ANALYSIS_RC=$?

    echo
    echo \
      "Analysis return code: $ANALYSIS_RC"

    if [ "$ANALYSIS_RC" -eq 0 ]; then

        echo
        echo "=== FINAL COMPARISON ==="

        cat \
          results/phase05d2_analysis/mesh_erosion_comparison.csv

        echo
        echo "=== FINAL SUMMARY ==="

        cat \
          results/phase05d2_analysis/phase05d2_summary.txt

    fi

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
