#!/usr/bin/env bash
set +e

ROOT="$(pwd)"
MODEL="$ROOT/model/phase04b_states"
RESULT="$ROOT/results/phase04b_states"

echo "========================================"
echo "PHASE 04B — ANTECEDENT STATES"
echo "========================================"

echo
echo "=== BUILD PROJECT FILES ==="

python \
  src/build_phase04b_states.py

BUILD_RC=$?

echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== COPY VERIFIED SLOPE MESH ==="

    cp \
      model/phase04a_slope/slope*.vtu \
      "$MODEL"/

    echo
    echo "=== MODEL FILES ==="

    ls -lh \
      "$MODEL"

    echo
    echo "=== RUN STATES ==="

    rm -rf "$RESULT"
    mkdir -p "$RESULT"

    ALL_PASS=1

    for STATE in dry reference wet
    do

        echo
        echo "----------------------------------------"
        echo "RUN STATE: $STATE"
        echo "----------------------------------------"

        mkdir -p \
          "$RESULT/$STATE"

        xmllint \
          --noout \
          "$MODEL/$STATE.prj"

        XML_RC=$?

        echo "XML return code: $XML_RC"

        if [ "$XML_RC" -ne 0 ]; then
            ALL_PASS=0
            continue
        fi

        ogs \
          "$MODEL/$STATE.prj" \
          -m "$MODEL" \
          -o "$RESULT/$STATE" \
          2>&1 \
          | tee \
          "$RESULT/$STATE/ogs.log"

        OGS_RC=${PIPESTATUS[0]}

        echo
        echo "$STATE OGS return code: $OGS_RC"

        if [ "$OGS_RC" -ne 0 ]; then

            ALL_PASS=0

            echo
            echo "=== FAILURE LOG ==="

            tail -n 50 \
              "$RESULT/$STATE/ogs.log"

        fi

    done

    echo
    echo "=== RUN STATUS ==="

    if [ "$ALL_PASS" -eq 1 ]; then

        echo "ALL THREE STATES: PASS"

        echo
        echo "=== ANALYZE CANDIDATE As ==="

        python \
          src/analyze_phase04b_states.py

        echo
        echo "=== FINAL SUMMARY ==="

        cat \
          results/phase04b_analysis/phase04b_summary.txt

    else

        echo "ONE OR MORE STATES FAILED"
        echo "Terminal remains open."

    fi
fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
