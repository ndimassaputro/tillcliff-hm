#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

MODEL="$ROOT/model/phase04d_mc_hm"
OUT="$ROOT/results/phase04d_mc_hm"

echo "========================================"
echo "PHASE 04D — COUPLED HM + MOHR-COULOMB"
echo "========================================"

echo
echo "=== BUILD PROJECT ==="

python \
  src/build_phase04d_mc_hm.py

BUILD_RC=$?

echo
echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== COPY VERIFIED SLOPE MESH ==="

    rm -f "$MODEL"/slope*.vtu

    cp \
      model/phase04a_slope/slope*.vtu \
      "$MODEL"/

    echo
    echo "=== FILE CHECK ==="

    ls -lh \
      "$MODEL"/mc_hm_baseline.prj \
      "$MODEL"/slope*.vtu

    echo
    echo "=== XML CHECK ==="

    xmllint \
      --noout \
      "$MODEL/mc_hm_baseline.prj"

    XML_RC=$?

    echo "XML return code: $XML_RC"

    if [ "$XML_RC" -eq 0 ]; then

        echo
        echo "=== CLEAN OUTPUT ==="

        rm -rf "$OUT"
        mkdir -p "$OUT"

        echo
        echo "=== RUN OGS ==="

        ogs \
          "$MODEL/mc_hm_baseline.prj" \
          -m "$MODEL" \
          -o "$OUT" \
          2>&1 \
          | tee "$OUT/ogs.log"

        OGS_RC=${PIPESTATUS[0]}

        echo
        echo "=== OGS RETURN CODE ==="
        echo "$OGS_RC"

        if [ "$OGS_RC" -eq 0 ]; then

            echo
            echo "=== ANALYZE ==="

            python \
              src/analyze_phase04d_mc_hm.py

            ANALYSIS_RC=$?

            echo
            echo "Analysis return code: $ANALYSIS_RC"

            if [ "$ANALYSIS_RC" -eq 0 ]; then

                echo
                echo "=== SUMMARY ==="

                cat \
                  "$OUT/phase04d_summary.txt"

            fi

        else

            echo
            echo "=== OGS FAILURE LOG ==="

            tail -n 80 \
              "$OUT/ogs.log"

        fi
    fi
fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
