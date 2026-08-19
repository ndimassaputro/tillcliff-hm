#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

MODEL="$ROOT/model/phase04a_slope"
OUT="$ROOT/results/phase04a_gravity"

echo "========================================"
echo "PHASE 04A — 2D SLOPE + GRAVITY"
echo "========================================"

echo
echo "=== GENERATE RECTANGULAR PARENT MESH ==="

rm -f "$MODEL"/slope*.vtu

generateStructuredMesh \
  -e quad \
  -o "$MODEL/slope.vtu" \
  --lx 30 \
  --ly 10 \
  --nx 150 \
  --ny 50

MESH_RC=$?

echo "Mesh return code: $MESH_RC"

if [ "$MESH_RC" -ne 0 ]; then
    echo "MESH GENERATION FAILED"
else

    echo
    echo "=== TRANSFORM TO COASTAL SLOPE ==="

    python \
      src/build_phase04a_slope.py

    BUILD_RC=$?

    echo
    echo "Build return code: $BUILD_RC"

    if [ "$BUILD_RC" -eq 0 ]; then

        echo
        echo "=== MESH FILES ==="

        ls -lh \
          "$MODEL"/slope*.vtu

        echo
        echo "=== XML CHECK ==="

        xmllint \
          --noout \
          "$MODEL/gravity_baseline.prj"

        XML_RC=$?

        echo "XML return code: $XML_RC"

        echo
        echo "=== CLEAN OUTPUT ==="

        rm -rf "$OUT"
        mkdir -p "$OUT"

        echo
        echo "=== RUN OGS ==="

        ogs \
          "$MODEL/gravity_baseline.prj" \
          -m "$MODEL" \
          -o "$OUT" \
          2>&1 \
          | tee "$OUT/ogs.log"

        OGS_RC=${PIPESTATUS[0]}

        echo
        echo "OGS return code: $OGS_RC"

        if [ "$OGS_RC" -eq 0 ]; then

            echo
            echo "=== ANALYZE ==="

            python \
              src/analyze_phase04a_gravity.py

            echo
            echo "=== SUMMARY ==="

            cat \
              "$OUT/phase04a_summary.txt"

        else

            echo
            echo "=== OGS FAILURE LOG ==="

            tail -n 60 \
              "$OUT/ogs.log"

        fi
    fi
fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
