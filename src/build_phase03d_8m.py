from pathlib import Path
import xml.etree.ElementTree as ET

src = Path("model/seasonal_column_3cycle.prj")
dst = Path("model/phase03d_8m/seasonal_column_8m_3cycle.prj")

if not src.exists():
    raise SystemExit(
        "FAIL: model/seasonal_column_3cycle.prj not found"
    )

text = src.read_text(encoding="utf-8")

replacements = {
    "<mesh>column.vtu</mesh>":
        "<mesh>column8.vtu</mesh>",

    "<mesh>column_left.vtu</mesh>":
        "<mesh>column8_left.vtu</mesh>",

    "<mesh>column_right.vtu</mesh>":
        "<mesh>column8_right.vtu</mesh>",

    "<mesh>column_top.vtu</mesh>":
        "<mesh>column8_top.vtu</mesh>",

    "<mesh>column_bottom.vtu</mesh>":
        "<mesh>column8_bottom.vtu</mesh>",

    "<mesh>column_left</mesh>":
        "<mesh>column8_left</mesh>",

    "<mesh>column_right</mesh>":
        "<mesh>column8_right</mesh>",

    "<mesh>column_top</mesh>":
        "<mesh>column8_top</mesh>",

    "<mesh>column_bottom</mesh>":
        "<mesh>column8_bottom</mesh>",

    "<prefix>seasonal_column_3cycle</prefix>":
        "<prefix>seasonal_column_8m_3cycle</prefix>",
}

for old, new in replacements.items():
    if old not in text:
        raise SystemExit(
            f"FAIL: expected text not found:\n{old}"
        )

    text = text.replace(
        old,
        new,
        1,
    )

dst.write_text(
    text,
    encoding="utf-8",
)

ET.parse(dst)

print("PASS:", dst)
print("Domain height: 8 m")
print("Analysis zone will remain: 0-4 m depth")
