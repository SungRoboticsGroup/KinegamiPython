#!/usr/bin/env python3
"""
Generate an Arduino .ino file from a Kinegami .session file.

Reads the saved configuration sequence and durations, identifies the real
(non-Waypoint, non-Tip) joints, and emits Arduino code that interpolates
between configurations and drives 270° servos via the Servo library.

Usage:
    python generateArduino.py save/dance_party_study/beach_pepe.session
    python generateArduino.py save/dance_party_study/beach_pepe.session -o pepe_dance.ino
"""

import sys, os, argparse, textwrap
import dill

# Ensure project root is importable
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Joint import Waypoint, Tip, Revolute, Prismatic


def real_joint_info(tree):
    """Return list of (tree_index, joint_type_name) for actuated joints."""
    info = []
    for i, j in enumerate(tree.Joints):
        if isinstance(j, (Waypoint, Tip)):
            continue
        info.append((i, type(j).__name__))
    return info


def generate_arduino(session_path: str, output_path: str | None = None,
                     default_pin_start: int = 3):
    """Load a .session file and write an Arduino .ino controlling the dance."""

    # ── Load session ─────────────────────────────────────────────────────
    with open(session_path, "rb") as f:
        state = dill.load(f)

    tree = state.get("tree")
    if tree is None:
        raise ValueError(f"No tree found in {session_path}")

    configs = state.get("saved_configurations", [])
    durations = state.get("config_durations", [])

    if len(configs) < 2:
        raise ValueError("Need at least 2 saved configurations for a sequence")

    joints = real_joint_info(tree)
    num_servos = len(joints)
    num_configs = len(configs)
    num_segments = num_configs - 1

    # Pad durations if shorter than segments
    while len(durations) < num_segments:
        durations.append(1.0)

    # ── Derive output filename ───────────────────────────────────────────
    if output_path is None:
        base = os.path.splitext(os.path.basename(session_path))[0]
        output_path = base + ".ino"

    # ── Format helpers ───────────────────────────────────────────────────
    def fmt_row(values):
        """Format a config row as C array literal contents."""
        return ", ".join(f"{v:8.2f}" for v in values)

    # ── Build Arduino source ─────────────────────────────────────────────
    lines = []
    w = lines.append  # shorthand

    w("// ===================================================================")
    w(f"// Auto-generated from: {os.path.basename(session_path)}")
    w(f"// Servos: {num_servos}   Configurations: {num_configs}")
    w("//")
    w("// 270-degree servos: physical range [0, 270]")
    w("//   Servo library signal range [0, 180] maps linearly to [0, 270].")
    w("//   Joint states are expressed in [-135, 135] (0 = center).")
    w("//   Conversion:  signal = (SIGN * state + 135) * 180.0 / 270.0")
    w("// ===================================================================")
    w("")
    w('#include <Servo.h>')
    w("")

    # ── Constants ────────────────────────────────────────────────────────
    w(f"#define NUM_SERVOS  {num_servos}")
    w(f"#define NUM_CONFIGS {num_configs}")
    w(f"#define NUM_SEGMENTS (NUM_CONFIGS - 1)")
    w("")

    # Pin assignments
    w("// -- Pin assignments (edit these to match your wiring) ----------------")
    w("const int SERVO_PINS[NUM_SERVOS] = {")
    pin_entries = []
    for idx, (tree_idx, jtype) in enumerate(joints):
        pin = default_pin_start + idx
        comma = "," if idx < num_servos - 1 else ""
        pin_entries.append(f"  {pin:3d}{comma}  // servo {idx}: joint {tree_idx} ({jtype})")
    w("\n".join(pin_entries))
    w("};")
    w("")

    # Sign per servo (+1 or -1 to reverse direction)
    w("// -- Direction signs (+1 normal, -1 reversed) -----------------------")
    w("const int SERVO_SIGN[NUM_SERVOS] = {")
    sign_entries = []
    for idx in range(num_servos):
        comma = "," if idx < num_servos - 1 else ""
        sign_entries.append(f"   1{comma}  // servo {idx}")
    w("\n".join(sign_entries))
    w("};")
    w("")

    # Configuration table (degrees, range [-135, 135])
    w("// -- Configuration sequence (degrees, range [-135, 135]) -------------")
    w("// Each row is one configuration; columns are servos 0 .. N-1.")
    w("// Edit these values to adjust the dance.")
    w("const float CONFIGS[NUM_CONFIGS][NUM_SERVOS] = {")
    for ci, cfg in enumerate(configs):
        # Pad or truncate to num_servos
        row = list(cfg[:num_servos])
        while len(row) < num_servos:
            row.append(0.0)
        comma = "," if ci < num_configs - 1 else ""
        w(f"  {{ {fmt_row(row)} }}{comma}  // config {ci}")
    w("};")
    w("")

    # Durations (seconds per segment, stored as milliseconds in Arduino)
    w("// -- Segment durations (milliseconds) --------------------------------")
    w("// Duration between config[i] and config[i+1].")
    w("const unsigned long SEGMENT_MS[NUM_SEGMENTS] = {")
    dur_entries = []
    for si in range(num_segments):
        ms = int(round(durations[si] * 1000))
        comma = "," if si < num_segments - 1 else ""
        dur_entries.append(f"  {ms:6d}{comma} // segment {si} -> {si+1}")
    w("\n".join(dur_entries))
    w("};")
    w("")

    # ── Servo objects ────────────────────────────────────────────────────
    w("// -- Servo objects ----------------------------------------------------")
    w("Servo servos[NUM_SERVOS];")
    w("")

    # ── Conversion function ──────────────────────────────────────────────
    w("// -- Convert joint state [-135, 135] -> servo signal [0, 180] --------")
    w("int stateToSignal(int servoIdx, float stateDeg) {")
    w("  // Apply direction sign")
    w("  float effective = SERVO_SIGN[servoIdx] * stateDeg;")
    w("  // Shift from [-135,135] to [0,270], then scale to [0,180]")
    w("  float signal = (effective + 135.0) * 180.0 / 270.0;")
    w("  // Clamp to valid servo range")
    w("  if (signal < 0.0)   signal = 0.0;")
    w("  if (signal > 180.0) signal = 180.0;")
    w("  return (int)(signal + 0.5);  // round to nearest int")
    w("}")
    w("")

    # ── setup() ──────────────────────────────────────────────────────────
    w("void setup() {")
    w("  Serial.begin(9600);")
    w("  for (int i = 0; i < NUM_SERVOS; i++) {")
    w("    servos[i].attach(SERVO_PINS[i]);")
    w("  }")
    w("  // Move to initial configuration")
    w("  for (int i = 0; i < NUM_SERVOS; i++) {")
    w("    servos[i].write(stateToSignal(i, CONFIGS[0][i]));")
    w("  }")
    w("  delay(500);  // let servos reach start position")
    w("}")
    w("")

    # ── loop() ───────────────────────────────────────────────────────────
    w("void loop() {")
    w("  for (int seg = 0; seg < NUM_SEGMENTS; seg++) {")
    w("    unsigned long segStart = millis();")
    w("    unsigned long duration = SEGMENT_MS[seg];")
    w("")
    w("    while (true) {")
    w("      unsigned long elapsed = millis() - segStart;")
    w("      if (elapsed >= duration) break;")
    w("")
    w("      float t = (float)elapsed / (float)duration;  // 0.0 -> 1.0")
    w("")
    w("      for (int i = 0; i < NUM_SERVOS; i++) {")
    w("        float val = CONFIGS[seg][i] * (1.0 - t) + CONFIGS[seg + 1][i] * t;")
    w("        servos[i].write(stateToSignal(i, val));")
    w("      }")
    w("")
    w("      delay(20);  // ~50 Hz update rate")
    w("    }")
    w("")
    w("    // Snap to exact end-of-segment position")
    w("    for (int i = 0; i < NUM_SERVOS; i++) {")
    w("      servos[i].write(stateToSignal(i, CONFIGS[seg + 1][i]));")
    w("    }")
    w("  }")
    w("  // Loop back to start (sequence repeats)")
    w("}")
    w("")

    # ── Write file ───────────────────────────────────────────────────────
    source = "\n".join(lines)
    with open(output_path, "w", encoding="ascii") as f:
        f.write(source)

    print(f"Generated {output_path}")
    print(f"  {num_servos} servos, {num_configs} configs, {num_segments} segments")
    print(f"  Total duration: {sum(durations[:num_segments]):.1f}s")
    for idx, (tree_idx, jtype) in enumerate(joints):
        print(f"  servo {idx}: pin {default_pin_start + idx}, "
              f"joint[{tree_idx}] ({jtype}), sign +1")


# ── CLI ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Arduino .ino from a Kinegami .session file")
    parser.add_argument("session", help="Path to .session file")
    parser.add_argument("-o", "--output", default=None,
                        help="Output .ino file path (default: <session_name>.ino)")
    parser.add_argument("--pin-start", type=int, default=3,
                        help="First servo pin number (default: 3)")
    args = parser.parse_args()

    generate_arduino(args.session, args.output, args.pin_start)
