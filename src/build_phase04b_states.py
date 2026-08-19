from pathlib import Path
import re
import xml.etree.ElementTree as ET

SRC = Path(
    "model/phase04a_slope/gravity_baseline.prj"
)

OUT = Path(
    "model/phase04b_states"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

if not SRC.exists():
    raise SystemExit(
        f"FAIL: missing {SRC}"
    )

base = SRC.read_text(
    encoding="utf-8"
)

states = {
    "dry": 1.20,        # -120 kPa
    "reference": 1.00,  # -100 kPa
    "wet": 0.30,        #  -30 kPa
}

for state, multiplier in states.items():

    text = base

    # --------------------------------------------------------
    # 120 d total:
    # 0-30 d gravity/reference equilibration
    # 30-60 d hydraulic transition
    # 60-120 d hold target state
    # --------------------------------------------------------

    text, n = re.subn(
        r"<t_end>\s*2592000\s*</t_end>",
        "<t_end>10368000</t_end>",
        text,
        count=1,
    )

    if n != 1:
        raise SystemExit(
            f"FAIL: t_end replacement for {state}"
        )

    old_steps = r"""
                    <timesteps>

                        <pair>
                            <repeat>1</repeat>
                            <delta_t>3600</delta_t>
                        </pair>

                        <pair>
                            <repeat>5</repeat>
                            <delta_t>18000</delta_t>
                        </pair>

                        <pair>
                            <repeat>29</repeat>
                            <delta_t>86400</delta_t>
                        </pair>

                    </timesteps>
"""

    new_steps = r"""
                    <timesteps>

                        <pair>
                            <repeat>1</repeat>
                            <delta_t>3600</delta_t>
                        </pair>

                        <pair>
                            <repeat>5</repeat>
                            <delta_t>18000</delta_t>
                        </pair>

                        <pair>
                            <repeat>119</repeat>
                            <delta_t>86400</delta_t>
                        </pair>

                    </timesteps>
"""

    if old_steps not in text:
        raise SystemExit(
            f"FAIL: timestep block not found for {state}"
        )

    text = text.replace(
        old_steps,
        new_steps,
        1,
    )

    text = text.replace(
        "<prefix>\n                phase04a_gravity\n            </prefix>",
        (
            "<prefix>\n"
            f"                phase04b_{state}\n"
            "            </prefix>"
        ),
        1,
    )

    # --------------------------------------------------------
    # Replace constant surface pressure with curve-scaled BC.
    # Reference magnitude remains -100 kPa.
    # --------------------------------------------------------

    old_parameter = r"""
        <parameter>
            <name>surface_pressure</name>
            <type>Constant</type>
            <value>-100000</value>
        </parameter>
"""

    new_parameter = r"""
        <parameter>
            <name>surface_pressure_scale</name>
            <type>Constant</type>
            <value>-100000</value>
        </parameter>

        <parameter>
            <name>surface_pressure</name>
            <type>CurveScaled</type>
            <curve>antecedent_state_curve</curve>
            <parameter>surface_pressure_scale</parameter>
        </parameter>
"""

    if old_parameter not in text:
        raise SystemExit(
            f"FAIL: surface pressure parameter not found for {state}"
        )

    text = text.replace(
        old_parameter,
        new_parameter,
        1,
    )

    curves = f"""
    <curves>

        <curve>
            <name>antecedent_state_curve</name>

            <!--
            0-30 d:
                common -100 kPa reference state

            30-60 d:
                smooth linear transition

            60-120 d:
                hold final antecedent state
            -->

            <coords>
                0
                2592000
                5184000
                10368000
            </coords>

            <values>
                1.0
                1.0
                {multiplier:.8f}
                {multiplier:.8f}
            </values>

        </curve>

    </curves>

"""

    marker = "    <process_variables>\n"

    if marker not in text:
        raise SystemExit(
            f"FAIL: process_variables marker absent for {state}"
        )

    text = text.replace(
        marker,
        curves + marker,
        1,
    )

    dst = OUT / f"{state}.prj"

    dst.write_text(
        text,
        encoding="utf-8",
    )

    ET.parse(dst)

    print(
        f"PASS {state:9s} -> "
        f"target surface suction "
        f"{100*multiplier:.1f} kPa"
    )

print()
print("PHASE 04B PROJECT BUILD: PASS")
