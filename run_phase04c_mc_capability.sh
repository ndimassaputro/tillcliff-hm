#!/usr/bin/env bash
set +e

ROOT="$(pwd)"
OGSROOT="$ROOT/external/ogs"

CASE="$OGSROOT/Tests/Data/Mechanics/MohrCoulombAbboSloan"
PRJ="$CASE/load_test_mc.prj"

OUT="$ROOT/results/phase04c_mc_capability"

echo "========================================"
echo "PHASE 04C — MOHR-COULOMB CAPABILITY"
echo "========================================"

echo
echo "=== ADD OFFICIAL OGS MC BENCHMARK ==="

git -C "$OGSROOT" sparse-checkout add \
    Tests/Data/Mechanics/MohrCoulombAbboSloan

GIT_RC=$?

echo "Sparse-checkout return code: $GIT_RC"

echo
echo "=== INPUT CHECK ==="

if [ -f "$PRJ" ]; then
    echo "FOUND: $PRJ"
else
    echo "MISSING: $PRJ"
fi

grep -n -A20 \
    '<type>MFront</type>' \
    "$PRJ" \
    | sed -n '1,30p'

echo
echo "=== CLEAN OUTPUT ==="

rm -rf "$OUT"
mkdir -p "$OUT"

echo
echo "=== RUN OFFICIAL MOHR-COULOMB TEST ==="

ogs \
    "$PRJ" \
    -m "$CASE" \
    -o "$OUT" \
    2>&1 \
    | tee "$OUT/ogs.log"

OGS_RC=${PIPESTATUS[0]}

echo
echo "=== OGS RETURN CODE ==="
echo "$OGS_RC"

echo
echo "=== OUTPUT COUNTS ==="

PVD_COUNT="$(
    find "$OUT" \
        -maxdepth 1 \
        -name '*.pvd' \
        | wc -l \
        | tr -d ' '
)"

VTU_COUNT="$(
    find "$OUT" \
        -maxdepth 1 \
        -name '*.vtu' \
        | wc -l \
        | tr -d ' '
)"

echo "PVD files: $PVD_COUNT"
echo "VTU files: $VTU_COUNT"

echo
echo "=== LOG TAIL ==="
tail -n 40 "$OUT/ogs.log"

echo
echo "========================================"

if [ "$OGS_RC" -eq 0 ] && \
   [ "$PVD_COUNT" -ge 1 ] && \
   [ "$VTU_COUNT" -ge 1 ]; then

    echo "PHASE 04C MOHR-COULOMB CAPABILITY: PASS"

else

    echo "PHASE 04C MOHR-COULOMB CAPABILITY: REVIEW"
    echo "Terminal remains open."

fi

echo "========================================"
