#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

MODEL="$ROOT/model/phase05b_v2_moving_front"
RUN="$ROOT/results/phase05b_v2_moving_front"

echo "========================================"
echo "PHASE 05B-V2 — MOVING EROSION FRONT"
echo "========================================"

echo
echo "=== BUILD ==="

python \
  src/build_phase05b_v2_moving_front.py

BUILD_RC=$?

echo
echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== FILE CHECK ==="

    wc -l \
      "$MODEL/reference_moving_front.prj"

    ls -lh \
      "$MODEL/reference_moving_front.prj" \
      "$MODEL/slope_toe_notch_ready.vtu"

    echo
    echo "=== XML CHECK ==="

    xmllint \
      --noout \
      "$MODEL/reference_moving_front.prj"

    XML_RC=$?

    echo \
      "XML return code: $XML_RC"

    if [ "$XML_RC" -eq 0 ]; then

        echo
        echo "=== CLEAN OUTPUT ==="

        rm -rf "$RUN"
        mkdir -p "$RUN"

        echo
        echo "=== RUN OGS ==="

        ogs \
          "$MODEL/reference_moving_front.prj" \
          -m "$MODEL" \
          -o "$RUN" \
          2>&1 \
          | tee \
          "$RUN/ogs.log"

        OGS_RC=${PIPESTATUS[0]}

        echo
        echo "=== OGS RETURN CODE ==="
        echo "$OGS_RC"

        echo
        echo "=== ANALYZE VALID ACCEPTED STATES ==="

        if find "$RUN" \
          -maxdepth 1 \
          -name '*.pvd' \
          | grep -q .
        then

            python \
              src/analyze_phase05b_v2_moving_front.py

            ANALYSIS_RC=$?

            echo
            echo \
              "Analysis return code: $ANALYSIS_RC"

            if [ "$ANALYSIS_RC" -eq 0 ]; then

                echo
                echo "=== FINAL SUMMARY ==="

                cat \
                  results/phase05b_v2_analysis/phase05b_v2_summary.txt

            fi

        else

            echo \
              "NO PVD AVAILABLE FOR ANALYSIS"

        fi

        if [ "$OGS_RC" -ne 0 ]; then

            echo
            echo "=== OGS FAILURE TAIL ==="

            tail -n 60 \
              "$RUN/ogs.log"

        fi
    fi
fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
