#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
SRCROOT="$ROOT/external/ogs"
SRC="$SRCROOT/Tests/Data/Parabolic/Richards"
OUT="$ROOT/results/phase01_official_richards"

echo "=== ENVIRONMENT ==="
command -v python
python --version
command -v ogs
ogs --version

echo
echo "=== GET OFFICIAL OGS BENCHMARK ==="

mkdir -p "$ROOT/external" "$ROOT/results"

if [ ! -d "$SRCROOT/.git" ]; then
    git clone \
        --depth 1 \
        --filter=blob:none \
        --sparse \
        https://github.com/ufz/ogs.git \
        "$SRCROOT"
fi

git -C "$SRCROOT" sparse-checkout set \
    Tests/Data/Parabolic/Richards

PRJ="$SRC/RichardsFlow_2d_small.prj"

if [ ! -f "$PRJ" ]; then
    echo "FAIL: benchmark project file not found:"
    echo "$PRJ"
    exit 1
fi

echo "PASS: benchmark input found"
echo "$PRJ"

echo
echo "=== INPUT CHECK ==="
ls -lh "$PRJ"

echo
echo "Relevant Richards files:"
find "$SRC" \
    -maxdepth 1 \
    -type f \
    \( -name '*.prj' -o -name '*.vtu' -o -name '*.gml' \) \
    | sort \
    | sed -n '1,40p'

echo
echo "=== CLEAN OUTPUT ==="
rm -rf "$OUT"
mkdir -p "$OUT"

echo
echo "=== RUN OFFICIAL RICHARDS BENCHMARK ==="

ogs \
    "$PRJ" \
    -m "$SRC" \
    -o "$OUT" \
    2>&1 | tee "$OUT/ogs.log"

echo
echo "=== RESULT CHECK ==="

PVD_COUNT="$(
    find "$OUT" -maxdepth 1 -type f -name '*.pvd' | wc -l | tr -d ' '
)"

VTU_COUNT="$(
    find "$OUT" -maxdepth 1 -type f -name '*.vtu' | wc -l | tr -d ' '
)"

echo "PVD files: $PVD_COUNT"
echo "VTU files: $VTU_COUNT"

if [ "$PVD_COUNT" -lt 1 ]; then
    echo "FAIL: no PVD result file generated"
    exit 1
fi

if [ "$VTU_COUNT" -lt 1 ]; then
    echo "FAIL: no VTU result files generated"
    exit 1
fi

echo
echo "Generated files:"
find "$OUT" \
    -maxdepth 1 \
    -type f \
    | sort \
    | sed -n '1,40p'

echo
echo "=== LOG TAIL ==="
tail -n 30 "$OUT/ogs.log"

echo
echo "========================================"
echo "PHASE 01 RICHARDS BENCHMARK: PASS"
echo "========================================"
