#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

MODEL="$ROOT/model/phase05e4_decomposition"
RUN="$ROOT/results/phase05e4_decomposition"

echo "========================================"
echo "PHASE 05E-4 — P/SIGMA DECOMPOSITION"
echo "========================================"

echo
echo "=== BUILD ==="

python \
  src/build_phase05e4_decomposition.py

BUILD_RC=$?

echo
echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== WRITE CHECK ==="

    for CASE in \
      refP_refS \
      dryP_refS \
      wetP_refS \
      refP_dryS \
      refP_wetS
    do

        xmllint \
          --noout \
          "$MODEL/$CASE/decomposition.prj"

        echo \
          "$CASE XML RC: $?"

        ls -lh \
          "$MODEL/$CASE/bulk.vtu" \
          "$MODEL/$CASE/decomposition.prj"

    done

    echo
    echo "=== CLEAN OUTPUT ==="

    rm -rf "$RUN"
    mkdir -p "$RUN"

    cat > "$RUN/run_status.csv" <<'CSV'
case,return_code
CSV

    echo
    echo "=== RUN FIVE FACTORIAL CASES ==="

    for CASE in \
      refP_refS \
      dryP_refS \
      wetP_refS \
      refP_dryS \
      refP_wetS
    do

        echo
        echo "========================================"
        echo "RUN $CASE"
        echo "========================================"

        OUT="$RUN/$CASE"

        mkdir -p "$OUT"

        ogs \
          "$MODEL/$CASE/decomposition.prj" \
          -m "$MODEL/$CASE" \
          -o "$OUT" \
          2>&1 \
          | tee "$OUT/ogs.log"

        RC=${PIPESTATUS[0]}

        echo
        echo \
          "$CASE OGS RETURN CODE: $RC"

        echo \
          "$CASE,$RC" \
          >> "$RUN/run_status.csv"

        if [ "$RC" -ne 0 ]; then

            echo
            echo "=== FAILURE TAIL ==="

            tail -n 35 \
              "$OUT/ogs.log"

        fi

    done

    echo
    echo "=== RUN STATUS ==="

    cat \
      "$RUN/run_status.csv"

    echo
    echo "=== ANALYZE ==="

    python \
      src/analyze_phase05e4_decomposition.py

    ANALYSIS_RC=$?

    echo
    echo \
      "Analysis return code: $ANALYSIS_RC"

    if [ "$ANALYSIS_RC" -eq 0 ]; then

        echo
        echo "=== FINAL SUMMARY ==="

        cat \
          results/phase05e4_analysis/phase05e4_summary.txt

    fi
fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
