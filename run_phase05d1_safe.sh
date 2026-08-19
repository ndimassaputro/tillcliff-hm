#!/usr/bin/env bash
set +e

ROOT="$(pwd)"

MODEL="$ROOT/model/phase05d1_restart_mesh"
RUN="$ROOT/results/phase05d1_restart_mesh"

echo "========================================"
echo "PHASE 05D-1 — REFINED RESTART CHECK"
echo "========================================"

echo
echo "=== BUILD REFINED BULK MESHES ==="

python \
  src/build_phase05d1_refined_restart_meshes.py

BUILD_RC=$?

echo
echo "Build return code: $BUILD_RC"

if [ "$BUILD_RC" -eq 0 ]; then

    echo
    echo "=== BULK FILE CHECK ==="

    for CASE in coarse medium fine
    do

        ls -lh \
          "$MODEL/$CASE/bulk.vtu" \
          "$MODEL/$CASE/restart_check.prj" \
          "$MODEL/$CASE/erosion_events.csv"

    done

    echo
    echo "=== EXTRACT FULL BOUNDARIES ==="

    for CASE in coarse medium fine
    do

        echo
        echo "--- $CASE ---"

        ExtractBoundary \
          -i "$MODEL/$CASE/bulk.vtu" \
          -o "$MODEL/$CASE/boundary_all.vtu"

        echo \
          "$CASE ExtractBoundary RC: $?"

    done

    echo
    echo "=== SPLIT BOUNDARIES ==="

    python \
      src/split_phase05d1_boundaries.py

    SPLIT_RC=$?

    echo
    echo "Boundary split return code: $SPLIT_RC"

    if [ "$SPLIT_RC" -eq 0 ]; then

        echo
        echo "=== BOUNDARY FILE CHECK ==="

        for CASE in coarse medium fine
        do

            ls -lh \
              "$MODEL/$CASE/slope_left.vtu" \
              "$MODEL/$CASE/slope_right.vtu" \
              "$MODEL/$CASE/slope_top.vtu" \
              "$MODEL/$CASE/slope_bottom.vtu"

        done

        echo
        echo "=== XML CHECK ==="

        for CASE in coarse medium fine
        do

            xmllint \
              --noout \
              "$MODEL/$CASE/restart_check.prj"

            echo \
              "$CASE XML RC: $?"

        done

        echo
        echo "=== CLEAN RUN OUTPUT ==="

        rm -rf "$RUN"
        mkdir -p "$RUN"

        cat > "$RUN/run_status.csv" <<'CSV'
case,return_code
CSV

        echo
        echo "=== RUN INTACT RESTART TESTS ==="

        for CASE in coarse medium fine
        do

            echo
            echo "========================================"
            echo "RUN $CASE"
            echo "========================================"

            CASE_OUT="$RUN/$CASE"

            mkdir -p \
              "$CASE_OUT"

            ogs \
              "$MODEL/$CASE/restart_check.prj" \
              -m "$MODEL/$CASE" \
              -o "$CASE_OUT" \
              2>&1 \
              | tee \
              "$CASE_OUT/ogs.log"

            OGS_RC=${PIPESTATUS[0]}

            echo
            echo \
              "$CASE OGS RETURN CODE: $OGS_RC"

            echo \
              "$CASE,$OGS_RC" \
              >> "$RUN/run_status.csv"

            if [ "$OGS_RC" -ne 0 ]; then

                echo
                echo \
                  "=== $CASE FAILURE TAIL ==="

                tail -n 45 \
                  "$CASE_OUT/ogs.log"

            fi

        done

        echo
        echo "=== RUN STATUS ==="

        cat \
          "$RUN/run_status.csv"

        echo
        echo "=== ANALYZE ==="

        python \
          src/analyze_phase05d1_restart_equilibrium.py

        ANALYSIS_RC=$?

        echo
        echo \
          "Analysis return code: $ANALYSIS_RC"

        if [ "$ANALYSIS_RC" -eq 0 ]; then

            echo
            echo "=== MESH REFINEMENT SUMMARY ==="

            cat \
              results/phase05d1_analysis/mesh_refinement_summary.csv

            echo
            echo "=== FINAL SUMMARY ==="

            cat \
              results/phase05d1_analysis/phase05d1_summary.txt

        fi
    fi
fi

echo
echo "========================================"
echo "SCRIPT FINISHED — TERMINAL STAYS OPEN"
echo "========================================"
