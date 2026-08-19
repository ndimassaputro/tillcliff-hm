#!/usr/bin/env bash
set +e

ROOT="$(pwd)"
MODEL="$ROOT/model/phase03d_8m"
OUT="$ROOT/results/phase03d_8m"

echo "========================================"
echo "PHASE 03D — 8 m DOMAIN TEST"
echo "========================================"

echo
echo "=== GENERATE 8 m MESH ==="

rm -f \
  "$MODEL"/column8*.vtu

generateStructuredMesh \
  -e quad \
  -o "$MODEL/column8.vtu" \
  --lx 1 \
  --ly 8 \
  --nx 10 \
  --ny 160

MESH_RC=$?

echo "Mesh return code: $MESH_RC"

echo
echo "=== MESH CHECK ==="

find "$MODEL" \
  -maxdepth 1 \
  -name '*.vtu' \
  -print \
  | sort

if [ "$MESH_RC" -ne 0 ]; then
    echo "MESH GENERATION FAILED"
else

    echo
    echo "=== BUILD PROJECT ==="

    python \
      src/build_phase03d_8m.py

    BUILD_RC=$?

    echo "Build return code: $BUILD_RC"

    if [ "$BUILD_RC" -eq 0 ]; then

        echo
        echo "=== XML CHECK ==="

        xmllint \
          --noout \
          "$MODEL/seasonal_column_8m_3cycle.prj"

        echo
        echo "=== CLEAN OUTPUT ==="

        rm -rf "$OUT"
        mkdir -p "$OUT"

        echo
        echo "=== RUN OGS ==="

        ogs \
          "$MODEL/seasonal_column_8m_3cycle.prj" \
          -m "$MODEL" \
          -o "$OUT" \
          2>&1 \
          | tee "$OUT/ogs.log"

        OGS_RC=${PIPESTATUS[0]}

        echo
        echo "OGS return code: $OGS_RC"

        if [ "$OGS_RC" -eq 0 ]; then

            echo
            echo "=== ANALYZE BOUNDARY SENSITIVITY ==="

            python \
              src/analyze_phase03d_boundary.py

        else

            echo
            echo "=== FAILURE LOG ==="

            tail -n 50 \
              "$OUT/ogs.log"

        fi
    fi
fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
