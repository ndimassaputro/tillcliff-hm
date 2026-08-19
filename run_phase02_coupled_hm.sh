#!/usr/bin/env bash
set -euo pipefail

ROOT="$(pwd)"
OGSROOT="$ROOT/external/ogs"
CASE="$OGSROOT/Tests/Data/RichardsMechanics"

OUT_LIA="$ROOT/results/phase02_liakopoulos_hm"
OUT_POR="$ROOT/results/phase02_deformation_porosity"

echo "=== ENVIRONMENT ==="
python --version
ogs --version

echo
echo "=== ADD OFFICIAL RICHARDS-MECHANICS BENCHMARKS ==="

if [ ! -d "$OGSROOT/.git" ]; then
    echo "FAIL: external/ogs repository not found"
    exit 1
fi

git -C "$OGSROOT" sparse-checkout add \
    Tests/Data/RichardsMechanics

LIA="$CASE/LiakopoulosHM/liakopoulos.prj"
POR="$CASE/deformation_dependent_porosity.prj"

echo
echo "=== INPUT CHECK ==="

for f in "$LIA" "$POR"; do
    if [ ! -f "$f" ]; then
        echo "FAIL: missing benchmark:"
        echo "$f"
        exit 1
    fi

    echo "PASS: $f"
done

echo
echo "=== PROCESS-TYPE CHECK ==="

grep -n -A3 -B2 \
    '<type>RICHARDS_MECHANICS</type>' \
    "$LIA" \
    "$POR"

echo
echo "=== DEFORMATION-POROSITY OUTPUT CHECK ==="

grep -n \
    -E 'saturation|porosity|sigma|epsilon|velocity' \
    "$POR" \
    | sed -n '1,40p'

echo
echo "=== CLEAN OUTPUT ==="

rm -rf "$OUT_LIA" "$OUT_POR"
mkdir -p "$OUT_LIA" "$OUT_POR"

echo
echo "========================================"
echo "RUN 1: LIAKOPOULOS COUPLED HM"
echo "========================================"

ogs \
    "$LIA" \
    -m "$(dirname "$LIA")" \
    -o "$OUT_LIA" \
    2>&1 | tee "$OUT_LIA/ogs.log"

echo
echo "=== LIAKOPOULOS OUTPUT CHECK ==="

LIA_VTU="$(
    find "$OUT_LIA" \
        -maxdepth 1 \
        -type f \
        -name '*.vtu' \
        | wc -l \
        | tr -d ' '
)"

LIA_PVD="$(
    find "$OUT_LIA" \
        -maxdepth 1 \
        -type f \
        -name '*.pvd' \
        | wc -l \
        | tr -d ' '
)"

echo "VTU files: $LIA_VTU"
echo "PVD files: $LIA_PVD"

if [ "$LIA_VTU" -lt 2 ]; then
    echo "FAIL: Liakopoulos produced too few VTU files"
    exit 1
fi

if [ "$LIA_PVD" -lt 1 ]; then
    echo "FAIL: Liakopoulos produced no PVD"
    exit 1
fi

echo
echo "========================================"
echo "RUN 2: DEFORMATION-DEPENDENT POROSITY"
echo "========================================"

ogs \
    "$POR" \
    -m "$CASE" \
    -o "$OUT_POR" \
    2>&1 | tee "$OUT_POR/ogs.log"

echo
echo "=== POROSITY BENCHMARK OUTPUT CHECK ==="

POR_VTU="$(
    find "$OUT_POR" \
        -maxdepth 1 \
        -type f \
        -name '*.vtu' \
        | wc -l \
        | tr -d ' '
)"

POR_PVD="$(
    find "$OUT_POR" \
        -maxdepth 1 \
        -type f \
        -name '*.pvd' \
        | wc -l \
        | tr -d ' '
)"

echo "VTU files: $POR_VTU"
echo "PVD files: $POR_PVD"

if [ "$POR_VTU" -lt 2 ]; then
    echo "FAIL: porosity benchmark produced too few VTU files"
    exit 1
fi

if [ "$POR_PVD" -lt 1 ]; then
    echo "FAIL: porosity benchmark produced no PVD"
    exit 1
fi

echo
echo "=== OUTPUT FILE SAMPLE ==="

echo "--- Liakopoulos ---"
find "$OUT_LIA" \
    -maxdepth 1 \
    -type f \
    | sort \
    | sed -n '1,20p'

echo
echo "--- Deformation-dependent porosity ---"
find "$OUT_POR" \
    -maxdepth 1 \
    -type f \
    | sort \
    | sed -n '1,20p'

echo
echo "=== LOG TAILS ==="

echo "--- Liakopoulos ---"
tail -n 15 "$OUT_LIA/ogs.log"

echo
echo "--- Porosity ---"
tail -n 15 "$OUT_POR/ogs.log"

echo
echo "========================================"
echo "PHASE 02 COUPLED HM BENCHMARKS: PASS"
echo "========================================"
