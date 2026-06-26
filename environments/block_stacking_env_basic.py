from environment_shapes import *


def get_items():
    return [circle_outline(cx=300, cy=0, radius=100, color='blue'), #start
            circle_outline(cx=300, cy=200, radius=20, color='green'), #goal
            solid_box(300, 50, 10, size=20)]
