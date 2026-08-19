from pathlib import Path
import re
import xml.etree.ElementTree as ET

src = Path("model/seasonal_column.prj")
dst = Path("model/seasonal_column_3cycle.prj")

text = src.read_text(encoding="utf-8")

# ------------------------------------------------------------
# 30 d initialization + 3 x 365 d seasonal cycles
# total = 1125 d
# ------------------------------------------------------------

text, n1 = re.subn(
    r"<t_end>\s*34128000\s*</t_end>",
    "<t_end>97200000</t_end>",
    text,
    count=1,
)

text, n2 = re.subn(
    r"<repeat>\s*395\s*</repeat>\s*"
    r"<delta_t>\s*86400\s*</delta_t>",
    "<repeat>1125</repeat>\n"
    "                            <delta_t>86400</delta_t>",
    text,
    count=1,
)

text = text.replace(
    "<prefix>seasonal_column</prefix>",
    "<prefix>seasonal_column_3cycle</prefix>",
    1,
)

curve_pattern = re.compile(
    r"(<curve>\s*"
    r"<name>seasonal_suction</name>.*?"
    r"<coords>)(.*?)(</coords>.*?"
    r"<values>)(.*?)(</values>)",
    re.DOTALL,
)

coords_days = [
    0,
    30,
    90,
    180,
    270,
    360,
    395,
    455,
    545,
    635,
    725,
    760,
    820,
    910,
    1000,
    1090,
    1125,
]

values = [
    1.00,
    1.00,
    0.30,
    0.05,
    0.50,
    1.20,
    1.00,
    0.30,
    0.05,
    0.50,
    1.20,
    1.00,
    0.30,
    0.05,
    0.50,
    1.20,
    1.00,
]

coords_seconds = [
    int(day * 86400)
    for day in coords_days
]

coords_text = "\n                " + "\n                ".join(
    str(v) for v in coords_seconds
) + "\n            "

values_text = "\n                " + "\n                ".join(
    f"{v:.2f}" for v in values
) + "\n            "


def replace_curve(match):
    return (
        match.group(1)
        + coords_text
        + match.group(3)
        + values_text
        + match.group(5)
    )


text, n3 = curve_pattern.subn(
    replace_curve,
    text,
    count=1,
)

if n1 != 1:
    raise SystemExit(
        f"FAIL: t_end replacement count = {n1}"
    )

if n2 != 1:
    raise SystemExit(
        f"FAIL: timestep replacement count = {n2}"
    )

if n3 != 1:
    raise SystemExit(
        f"FAIL: seasonal curve replacement count = {n3}"
    )

dst.write_text(
    text,
    encoding="utf-8",
)

# Hard XML validation.
ET.parse(dst)

print("PASS: wrote", dst)
print("Simulation duration: 1125 days")
print("Seasonal cycles: 3")
