#!/usr/bin/env bash
set +e

ROOT="$(pwd)"
MODEL="$ROOT/model/phase04d_v2"
OUT="$ROOT/results/phase04d_v2"

echo "========================================"
echo "PHASE 04D-V2 — INITIALIZED MC-HM"
echo "========================================"

echo
echo "=== BUILD INITIALIZED STATE ==="

python \
  src/build_phase04d_v2_initialized.py

BUILD_RC=$?

echo
echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== XML CHECK ==="

    for CASE in strong screening
    do

        xmllint \
          --noout \
          "$MODEL/$CASE.prj"

        echo \
          "$CASE XML return code: $?"

    done

    echo
    echo "=== CLEAN OUTPUT ==="

    rm -rf "$OUT"
    mkdir -p "$OUT"

    STRONG_PASS=0
    SCREENING_PASS=0

    echo
    echo "========================================"
    echo "RUN STRONG GUARDRAIL"
    echo "========================================"

    mkdir -p \
      "$OUT/strong"

    ogs \
      "$MODEL/strong.prj" \
      -m "$MODEL" \
      -o "$OUT/strong" \
      2>&1 \
      | tee \
      "$OUT/strong/ogs.log"

    STRONG_RC=${PIPESTATUS[0]}

    echo
    echo "STRONG OGS return code: $STRONG_RC"

    if [ "$STRONG_RC" -eq 0 ]; then

        STRONG_PASS=1

        echo \
          "STRONG INITIALIZATION: PASS"

        echo
        echo "========================================"
        echo "RUN SCREENING STRENGTH"
        echo "========================================"

        mkdir -p \
          "$OUT/screening"

        ogs \
          "$MODEL/screening.prj" \
          -m "$MODEL" \
          -o "$OUT/screening" \
          2>&1 \
          | tee \
          "$OUT/screening/ogs.log"

        SCREENING_RC=${PIPESTATUS[0]}

        echo
        echo \
          "SCREENING OGS return code: $SCREENING_RC"

        if [ "$SCREENING_RC" -eq 0 ]; then

            SCREENING_PASS=1

            echo \
              "SCREENING INITIALIZATION: PASS"

        else

            echo \
              "SCREENING INITIALIZATION: REVIEW"

            echo
            echo "=== SCREENING FAILURE TAIL ==="

            tail -n 60 \
              "$OUT/screening/ogs.log"

        fi

    else

        echo \
          "STRONG INITIALIZATION: REVIEW"

        echo
        echo "=== STRONG FAILURE TAIL ==="

        tail -n 60 \
          "$OUT/strong/ogs.log"

    fi

    echo
    echo "========================================"
    echo "DIAGNOSTIC RESULT"
    echo "========================================"

    echo \
      "STRONG_PASS=$STRONG_PASS"

    echo \
      "SCREENING_PASS=$SCREENING_PASS"

    echo
    echo "=== ANALYZE AVAILABLE RESULTS ==="

    python \
      src/analyze_phase04d_v2.py

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
