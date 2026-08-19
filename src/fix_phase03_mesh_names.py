from pathlib import Path

path = Path("model/seasonal_column.prj")
text = path.read_text(encoding="utf-8")

replacements = {
    "<mesh>left.vtu</mesh>":
        "<mesh>column_left.vtu</mesh>",
    "<mesh>right.vtu</mesh>":
        "<mesh>column_right.vtu</mesh>",
    "<mesh>top.vtu</mesh>":
        "<mesh>column_top.vtu</mesh>",
    "<mesh>bottom.vtu</mesh>":
        "<mesh>column_bottom.vtu</mesh>",

    "<mesh>left</mesh>":
        "<mesh>column_left</mesh>",
    "<mesh>right</mesh>":
        "<mesh>column_right</mesh>",
    "<mesh>top</mesh>":
        "<mesh>column_top</mesh>",
    "<mesh>bottom</mesh>":
        "<mesh>column_bottom</mesh>",
}

for old, new in replacements.items():
    text = text.replace(old, new)

path.write_text(text, encoding="utf-8")

print("PASS: Phase 03 mesh names updated")
