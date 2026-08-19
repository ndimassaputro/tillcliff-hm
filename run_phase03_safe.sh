#!/usr/bin/env bash
set +e

ROOT="$(pwd)"
OUT="$ROOT/results/phase03_seasonal_column"

echo "========================================"
echo "PHASE 03A SAFE RUN"
echo "========================================"

echo
echo "=== APPLY MESH-NAME FIX ==="
python src/fix_phase03_mesh_names.py

echo
echo "=== PROJECT MESH REFERENCES ==="
grep -n '<mesh>' model/seasonal_column.prj

echo
echo "=== REQUIRED FILE CHECK ==="

for f in \
    model/column.vtu \
    model/column_left.vtu \
    model/column_right.vtu \
    model/column_top.vtu \
    model/column_bottom.vtu
do
    if [ -f "$f" ]; then
        echo "FOUND $f"
    else
        echo "MISSING $f"
    fi
done

echo
echo "=== XML CHECK ==="
xmllint --noout model/seasonal_column.prj
XML_RC=$?
echo "XML return code: $XML_RC"

echo
echo "=== CLEAN OLD OUTPUT ==="
rm -rf "$OUT"
mkdir -p "$OUT"

echo
echo "=== RUN OGS ==="

ogs \
    model/seasonal_column.prj \
    -m model \
    -o "$OUT" \
    2>&1 | tee "$OUT/ogs.log"

OGS_RC=${PIPESTATUS[0]}

echo
echo "=== OGS RETURN CODE ==="
echo "$OGS_RC"

echo
echo "=== OUTPUT COUNTS ==="

PVD_COUNT="$(
    find "$OUT" -maxdepth 1 -name '*.pvd' \
    | wc -l | tr -d ' '
)"

VTU_COUNT="$(
    find "$OUT" -maxdepth 1 -name '*.vtu' \
    | wc -l | tr -d ' '
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

    echo "PHASE 03A OWN SEASONAL HM COLUMN: PASS"

else

    echo "PHASE 03A: NEEDS ONE MORE FIX"
    echo "Terminal remains open."

fi

echo "========================================"
