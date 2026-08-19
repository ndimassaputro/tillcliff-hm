#!/usr/bin/env bash
set +e

ROOT="$(pwd)"
OUT="$ROOT/results/phase03c_three_cycles"

echo "========================================"
echo "PHASE 03C — THREE SEASONAL CYCLES"
echo "========================================"

echo
echo "=== BUILD ==="

python \
  src/build_phase03c_three_cycles.py

BUILD_RC=$?

if [ "$BUILD_RC" -ne 0 ]; then
    echo "BUILD FAILED"
    echo "Terminal remains open."
else

    echo
    echo "=== WRITE CHECK ==="

    wc -l \
      model/seasonal_column_3cycle.prj

    ls -lh \
      model/seasonal_column_3cycle.prj

    echo
    echo "=== CLEAN OUTPUT ==="

    rm -rf "$OUT"
    mkdir -p "$OUT"

    echo
    echo "=== RUN OGS ==="

    ogs \
      model/seasonal_column_3cycle.prj \
      -m model \
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
          src/analyze_phase03c_periodicity.py

    else

        echo
        echo "=== FAILURE LOG ==="

        tail -n 50 \
          "$OUT/ogs.log"

        echo
        echo "Terminal remains open."

    fi
fi
