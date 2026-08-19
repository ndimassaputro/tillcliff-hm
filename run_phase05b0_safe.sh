#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

MODEL="$ROOT/model/phase05b0_single_erosion"
OUT="$ROOT/results/phase05b0_single_erosion"

echo "========================================"
echo "PHASE 05B-0 — SINGLE EROSION TEST"
echo "========================================"

echo
echo "=== BUILD ==="

python \
  src/build_phase05b0_single_erosion.py

BUILD_RC=$?

echo
echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== FILE CHECK ==="

    wc -l \
      "$MODEL/reference_E0p4.prj"

    ls -lh \
      "$MODEL/reference_E0p4.prj" \
      "$MODEL/slope_toe_notch_ready.vtu"

    echo
    echo "=== XML CHECK ==="

    xmllint \
      --noout \
      "$MODEL/reference_E0p4.prj"

    XML_RC=$?

    echo \
      "XML return code: $XML_RC"

    if [ "$XML_RC" -eq 0 ]; then

        echo
        echo "=== CLEAN OUTPUT ==="

        rm -rf "$OUT"
        mkdir -p "$OUT"

        echo
        echo "=== RUN OGS ==="

        ogs \
          "$MODEL/reference_E0p4.prj" \
          -m "$MODEL" \
          -o "$OUT" \
          2>&1 \
          | tee \
          "$OUT/ogs.log"

        OGS_RC=${PIPESTATUS[0]}

        echo
        echo "=== OGS RETURN CODE ==="
        echo "$OGS_RC"

        if [ "$OGS_RC" -eq 0 ]; then

            echo
            echo "=== ANALYZE ==="

            python \
              src/analyze_phase05b0_single_erosion.py

            ANALYSIS_RC=$?

            echo
            echo \
              "Analysis return code: $ANALYSIS_RC"

            if [ "$ANALYSIS_RC" -eq 0 ]; then

                echo
                echo "=== FINAL SUMMARY ==="

                cat \
                  "$OUT/phase05b0_summary.txt"

            fi

        else

            echo
            echo "=== FAILURE LOG ==="

            tail -n 80 \
              "$OUT/ogs.log"

        fi
    fi
fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
