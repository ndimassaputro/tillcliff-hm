#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

MODEL="$ROOT/model/phase05c_transition_refinement"
RUN_ROOT="$ROOT/results/phase05c_transition_refinement"

echo "========================================"
echo "PHASE 05C — TRANSITION REFINEMENT"
echo "========================================"

echo
echo "=== BUILD + GEOMETRY AUDIT ==="

python \
  src/build_phase05c_transition_refinement.py

BUILD_RC=$?

echo
echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== WRITE CHECK ==="

    wc -l \
      "$MODEL"/dt0p10.prj \
      "$MODEL"/dt0p05.prj \
      results/phase05c_analysis/deactivation_events_035_060.csv

    ls -lh \
      "$MODEL"/dt0p10.prj \
      "$MODEL"/dt0p05.prj \
      results/phase05c_analysis/deactivation_events_035_060.csv

    echo
    echo "=== XML CHECK ==="

    xmllint \
      --noout \
      "$MODEL/dt0p10.prj"

    echo \
      "dt0p10 XML RC: $?"

    xmllint \
      --noout \
      "$MODEL/dt0p05.prj"

    echo \
      "dt0p05 XML RC: $?"

    echo
    echo "=== CLEAN RUN OUTPUT ==="

    rm -rf \
      "$RUN_ROOT"

    mkdir -p \
      "$RUN_ROOT"

    for CASE in dt0p10 dt0p05
    do

        echo
        echo "========================================"
        echo "RUN $CASE"
        echo "========================================"

        CASE_OUT="$RUN_ROOT/$CASE"

        mkdir -p \
          "$CASE_OUT"

        ogs \
          "$MODEL/$CASE.prj" \
          -m "$MODEL" \
          -o "$CASE_OUT" \
          2>&1 \
          | tee \
          "$CASE_OUT/ogs.log"

        OGS_RC=${PIPESTATUS[0]}

        echo
        echo \
          "$CASE OGS RETURN CODE: $OGS_RC"

        if [ "$OGS_RC" -ne 0 ]; then

            echo
            echo "=== $CASE FAILURE TAIL ==="

            tail -n 35 \
              "$CASE_OUT/ogs.log"

        fi

    done

    echo
    echo "=== ANALYZE ==="

    python \
      src/analyze_phase05c_transition_refinement.py

    ANALYSIS_RC=$?

    echo
    echo \
      "Analysis return code: $ANALYSIS_RC"

    if [ "$ANALYSIS_RC" -eq 0 ]; then

        echo
        echo "=== FINAL SUMMARY ==="

        cat \
          results/phase05c_analysis/phase05c_summary.txt

    fi

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
