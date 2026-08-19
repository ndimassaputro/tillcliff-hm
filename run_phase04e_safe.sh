#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

MODEL="$ROOT/model/phase04e_strength_probe"
OUT="$ROOT/results/phase04e_strength_probe"

echo "========================================"
echo "PHASE 04E — STRENGTH REDUCTION PROBE"
echo "========================================"

echo
echo "=== BUILD CASES ==="

python \
  src/build_phase04e_strength_probe.py

BUILD_RC=$?

echo
echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== WRITE / MODEL CHECK ==="

    wc -l \
      "$MODEL"/cases.csv

    ls -lh \
      "$MODEL"/cases.csv \
      "$MODEL"/*.prj

    echo
    echo "=== CLEAN RESULTS ==="

    rm -rf "$OUT"
    mkdir -p "$OUT"

    cat > "$OUT/run_status.csv" <<'CSV'
case,srf,return_code
CSV

    echo
    echo "=== RUN STRENGTH CASES ==="

    while IFS=, read -r \
      CASE \
      SRF \
      COHESION \
      PHI \
      C_RATIO \
      PHI_RATIO
    do

        echo
        echo "----------------------------------------"
        echo "CASE: $CASE"
        echo "SRF : $SRF"
        echo "c   : $COHESION kPa"
        echo "phi : $PHI deg"
        echo "----------------------------------------"

        CASE_OUT="$OUT/$CASE"

        mkdir -p \
          "$CASE_OUT"

        xmllint \
          --noout \
          "$MODEL/$CASE.prj"

        XML_RC=$?

        echo \
          "XML return code: $XML_RC"

        if [ "$XML_RC" -ne 0 ]; then

            echo \
              "$CASE,$SRF,98" \
              >> "$OUT/run_status.csv"

            continue

        fi

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
          "$CASE OGS return code: $OGS_RC"

        echo \
          "$CASE,$SRF,$OGS_RC" \
          >> "$OUT/run_status.csv"

        if [ "$OGS_RC" -ne 0 ]; then

            echo
            echo "=== FAILURE TAIL ==="

            tail -n 40 \
              "$CASE_OUT/ogs.log"

        fi

    done < <(
        tail -n +2 \
          "$MODEL/cases.csv"
    )

    echo
    echo "=== RUN STATUS TABLE ==="

    cat \
      "$OUT/run_status.csv"

    echo
    echo "=== ANALYZE ==="

    python \
      src/analyze_phase04e_strength_probe.py

    ANALYSIS_RC=$?

    echo
    echo \
      "Analysis return code: $ANALYSIS_RC"

    if [ "$ANALYSIS_RC" -eq 0 ]; then

        echo
        echo "=== FINAL SUMMARY ==="

        cat \
          results/phase04e_analysis/phase04e_summary.txt

    fi

fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
