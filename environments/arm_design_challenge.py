from environment_shapes import *


def get_items():
    INCH = 25.4  # mm
    # Listed first so the translucent base is drawn over it.
    # Lowered slightly below z=0 to avoid flickering against the base's
    # bottom face and outline, which also sit at z=0.
    board = flat_rectangle(0, 0, sx=16 * INCH, sy=24 * INCH, z=-0.5,
                           color=(0.93, 0.85, 0.72, 1.0))
    # 8 cubes in a ring inside the start circle, resting on z=0
    START_X, START_Y, START_R = 320, 200, 100
    CUBE = 25
    RING_R = 60  # cube corners reach 60 + 25/sqrt(2) ~ 78 < START_R
    cubes = [solid_box(START_X + RING_R * np.cos(a),
                       START_Y + RING_R * np.sin(a),
                       CUBE / 2, size=CUBE)
             for a in np.linspace(0, 2 * np.pi, 8, endpoint=False)]
    # 2 more cubes side by side near the center, 10mm apart
    cubes += [solid_box(START_X + dx, START_Y, CUBE / 2, size=CUBE)
              for dx in (-17.5, 17.5)]
    return [board,
            circle_outline(cx=START_X, cy=START_Y, radius=START_R, color='blue'), #start
            circle_outline(cx=320, cy=-200, radius=100, color='green'), #goal
            *cubes,
            #base
            solid_cylinder(0, 0, 0, radius=100, length=33,
                           color=(0.5, 0.5, 0.5, 0.5), cols=64,
                           glOptions='translucent'),
            circle_outline(0, 0, radius=100, z=0, color=(0, 0, 0, 1)),
            circle_outline(0, 0, radius=100, z=33, color=(0, 0, 0, 1))]
