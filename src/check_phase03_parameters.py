import numpy as np

SR = 0.05
M = 0.235
PB = 123000.0

def saturation(pc):
    n = 1.0 / (1.0 - M)
    se = (1.0 + (pc / PB) ** n) ** (-M)
    return SR + (1.0 - SR) * se

print("=== VAN GENUCHTEN SCREENING CHECK ===")

for pc_kpa in [5, 27.5, 30, 50, 100, 108, 120]:
    sr = saturation(pc_kpa * 1000.0)
    print(
        f"pc = {pc_kpa:6.1f} kPa"
        f" -> Sr = {sr:.5f}"
    )

print()
print("Reference Avedore 1 m observations:")
print("~27.5 kPa suction -> Sr ~0.971")
print("~108  kPa suction -> Sr ~0.873")
