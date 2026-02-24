# -*- coding: utf-8 -*-
"""
Dance Party Viewer
==================
Loads all .session files for a chosen environment (ballroom, beach, or club)
from save/dance_party_study/, spaces the robots out so they don't overlap,
and displays them together in a single pyqtgraph OpenGL window.

Usage:
    python dancePartyViewer.py                   # prompts for environment
    python dancePartyViewer.py --env ballroom    # loads ballroom directly
    python dancePartyViewer.py --env beach
    python dancePartyViewer.py --env club
"""

import sys, os, glob, copy, argparse, math
import numpy as np
import dill

# Ensure the project root is on the path so local imports work
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QWidget,
                              QComboBox, QLabel, QHBoxLayout, QPushButton)
from PyQt5.QtCore import Qt, QTimer, QTime
from PyQt5.QtGui import QSurfaceFormat
import pyqtgraph as pg
import pyqtgraph.opengl as gl

from spatialmath import SE3
from KinematicTree import KinematicTree
from Joint import Prismatic, Revolute, Waypoint, Tip
from style import *

# ── Helpers ──────────────────────────────────────────────────────────────────

ENVIRONMENTS = ["ballroom", "beach", "club"]

def session_dir():
    """Return the absolute path to the dance_party_study folder."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "save", "dance_party_study")

def list_sessions(env: str) -> list[str]:
    """Return sorted list of .session file paths matching the environment."""
    pattern = os.path.join(session_dir(), f"{env}_*.session")
    return sorted(glob.glob(pattern))

def load_session(filepath: str) -> dict:
    """Load a .session file (dill-serialized dict) and return the state dict."""
    with open(filepath, "rb") as f:
        return dill.load(f)

def robot_name_from_path(filepath: str) -> str:
    """Extract a human-readable name from a session filename."""
    base = os.path.splitext(os.path.basename(filepath))[0]
    # Remove the env prefix (e.g. "club_alvin" -> "alvin")
    parts = base.split("_", 1)
    return parts[1] if len(parts) > 1 else base


# ── Viewer widget (wraps a GLViewWidget so tree.addToWidget works) ───────────

class SimpleTreeViewer:
    """
    Minimal wrapper that exposes a `plot_widget` attribute (a GLViewWidget)
    so that KinematicTree.addToWidget(self, ...) can call
    self.plot_widget.addItem(...) as expected.
    """
    def __init__(self, gl_widget: gl.GLViewWidget):
        self.plot_widget = gl_widget


# ── Per-robot animation state ────────────────────────────────────────────────

def real_joint_indices(tree: KinematicTree) -> list[int]:
    """Return indices of joints that are neither Waypoint nor Tip."""
    return [i for i, j in enumerate(tree.Joints)
            if not isinstance(j, (Waypoint, Tip))]


class RobotAnimator:
    """Tracks animation state for a single robot."""
    def __init__(self, tree: KinematicTree,
                 saved_configs: list, config_durations: list):
        self.tree = tree
        self.saved_configs = saved_configs          # display-unit values
        self.config_durations = config_durations    # seconds per segment
        self.joint_indices = real_joint_indices(tree)
        self.current_segment = 0
        self.segment_start_time = 0.0

    @property
    def num_segments(self):
        return max(0, len(self.saved_configs) - 1)

    def can_animate(self):
        return self.num_segments >= 1

    def reset(self, now: float):
        self.current_segment = 0
        self.segment_start_time = now

    def step(self, now: float):
        """Advance animation, looping. Returns True if tree was updated."""
        if not self.can_animate():
            return False

        elapsed = now - self.segment_start_time
        dur = self._segment_duration(self.current_segment)
        t = elapsed / dur if dur > 0 else 1.0

        # Advance segments, carrying over leftover time
        while t >= 1.0:
            leftover = (t - 1.0) * dur  # excess time in seconds
            self.current_segment += 1
            if self.current_segment >= self.num_segments:
                self.current_segment = 0  # loop
            dur = self._segment_duration(self.current_segment)
            self.segment_start_time = now - leftover
            t = leftover / dur if dur > 0 else 0.0

        self._interpolate(self.current_segment, t)
        return True

    # ── internals ──

    def _segment_duration(self, seg_idx: int) -> float:
        if seg_idx < len(self.config_durations):
            return max(self.config_durations[seg_idx], 0.01)
        return 1.0  # fallback

    def _interpolate(self, seg_idx: int, t: float):
        """Set tree joint states to interpolation between config seg_idx and seg_idx+1."""
        config_a = self.saved_configs[seg_idx]
        config_b = self.saved_configs[seg_idx + 1]
        for i, joint_idx in enumerate(self.joint_indices):
            if i >= len(config_a) or i >= len(config_b):
                break
            if joint_idx >= len(self.tree.Joints):
                continue
            val = config_a[i] * (1.0 - t) + config_b[i] * t
            joint = self.tree.Joints[joint_idx]
            if isinstance(joint, Revolute):
                val = math.radians(val)
            self.tree.setJointState(joint_idx, val)


# ── Main window ──────────────────────────────────────────────────────────────

class DancePartyWindow(QMainWindow):
    def __init__(self, env: str | None = None):
        super().__init__()
        self.setWindowTitle("Dance Party Viewer")
        self.resize(1200, 800)

        # ── Central widget / layout ──
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # ── Environment selector row ──
        selector_row = QHBoxLayout()
        selector_row.setContentsMargins(6, 2, 6, 2)
        selector_row.addWidget(QLabel("Environment:"))
        self.env_combo = QComboBox()
        self.env_combo.addItems(ENVIRONMENTS)
        self.env_combo.currentTextChanged.connect(self._on_env_changed)
        selector_row.addWidget(self.env_combo)
        self.info_label = QLabel("")
        selector_row.addWidget(self.info_label)
        selector_row.addStretch()

        # Play / Pause button
        self.play_button = QPushButton("▶  Play")
        self.play_button.setFixedWidth(80)
        self.play_button.clicked.connect(self._toggle_animation)
        selector_row.addWidget(self.play_button)

        layout.addLayout(selector_row, 0)  # stretch=0: selector takes minimum height

        # ── GL view ──
        fmt = QSurfaceFormat()
        fmt.setDepthBufferSize(24)
        fmt.setProfile(QSurfaceFormat.CompatibilityProfile)
        QSurfaceFormat.setDefaultFormat(fmt)

        self.gl_widget = gl.GLViewWidget()
        self.gl_widget.setBackgroundColor(backgroundColorDefault)
        layout.addWidget(self.gl_widget, 1)  # stretch=1: GL view fills remaining space

        # Ground-plane grid
        self.grid = gl.GLGridItem()
        self.grid.setSize(600, 600, 600)
        self.grid.setSpacing(10, 10, 10)
        self.grid.setColor(gridColorDefault)
        self.gl_widget.addItem(self.grid)

        # The wrapper object that tree.addToWidget expects
        self.viewer = SimpleTreeViewer(self.gl_widget)

        # Storage
        self.trees: list[KinematicTree] = []
        self.tree_names: list[str] = []
        self.animators: list[RobotAnimator] = []
        self.is_animating = False

        # Animation timer (~60 fps)
        self.animation_timer = QTimer()
        self.animation_timer.timeout.connect(self._animation_step)

        # Load initial environment
        if env and env in ENVIRONMENTS:
            self.env_combo.setCurrentText(env)
        else:
            self._on_env_changed(self.env_combo.currentText())

    # ── Environment switching ────────────────────────────────────────────

    def _on_env_changed(self, env: str):
        """Load and display all robots for the selected environment."""
        # Stop any running animation
        self._stop_animation()
        self._clear_scene()

        session_files = list_sessions(env)
        if not session_files:
            self.info_label.setText(f"No sessions found for '{env}'")
            return

        # Load trees + animation data from session files
        trees = []
        names = []
        session_configs = []   # parallel list of saved_configurations
        session_durations = [] # parallel list of config_durations
        for fpath in session_files:
            try:
                state = load_session(fpath)
                tree = state.get("tree", None)
                if tree is not None:
                    trees.append(copy.deepcopy(tree))
                    names.append(robot_name_from_path(fpath))
                    session_configs.append(state.get("saved_configurations", []))
                    session_durations.append(state.get("config_durations", []))
            except Exception as e:
                print(f"Warning: could not load {fpath}: {e}")

        if not trees:
            self.info_label.setText(f"No valid trees in '{env}' sessions")
            return

        # Recompute bounding balls and figure out spacing
        for t in trees:
            t.recomputeBoundingBall()

        # Arrange trees in a line along the X axis, spaced so bounding
        # spheres don't overlap (with a small gap of 2× radius between them).
        # First, shift each tree so its bounding-ball center is at the origin.
        for t in trees:
            center = t.boundingBall.c
            t.transformAll(SE3.Rt(np.eye(3), -center))
            t.recomputeBoundingBall()

        # Compute cumulative X positions
        GAP_FACTOR = 1.2  # extra spacing multiplier
        x_positions = []
        current_x = 0.0
        for i, t in enumerate(trees):
            r_i = t.boundingBall.r
            if i == 0:
                current_x = 0.0
            else:
                r_prev = trees[i - 1].boundingBall.r
                current_x += (r_prev + r_i) * GAP_FACTOR
            x_positions.append(current_x)

        # Center the whole group around x=0
        total_span = x_positions[-1] if len(x_positions) > 1 else 0.0
        x_offset = total_span / 2.0

        # Translate each tree to its position
        for t, xpos in zip(trees, x_positions):
            t.transformAll(SE3.Rt(np.eye(3), np.array([xpos - x_offset, 0.0, 0.0])))
            t.recomputeBoundingBall()

        # Add name labels (rendered as tiny text items near each robot)
        for t, name, xpos in zip(trees, names, x_positions):
            label_pos = t.boundingBall.c.copy()
            label_pos[2] += t.boundingBall.r + 2.0  # above the robot
            text_item = gl.GLTextItem(pos=label_pos, text=name, color=(0, 0, 0, 255))
            self.gl_widget.addItem(text_item)

        # Render trees into the GL widget
        for t in trees:
            t.addToWidget(self.viewer,
                          showJointSurface=True,
                          showLinkSurface=True,
                          showLinkPath=False,
                          showJointPoses=False,
                          showSpheres=False)

        self.trees = trees
        self.tree_names = names

        # Create per-robot animators
        self.animators = []
        anim_count = 0
        for t, cfgs, durs in zip(trees, session_configs, session_durations):
            anim = RobotAnimator(t, cfgs, durs)
            self.animators.append(anim)
            if anim.can_animate():
                anim_count += 1

        # Set camera to see everything
        self._fit_camera()
        self.info_label.setText(
            f"Loaded {len(trees)} robots for '{env}' "
            f"({anim_count} with dance sequences)"
        )

    def _clear_scene(self):
        """Remove all items from the GL widget except the grid."""
        items = list(self.gl_widget.items)
        for item in items:
            if item is not self.grid:
                self.gl_widget.removeItem(item)
        self.trees = []
        self.tree_names = []
        self.animators = []

    # ── Animation ────────────────────────────────────────────────────────

    def _toggle_animation(self):
        if self.is_animating:
            self._stop_animation()
        else:
            self._start_animation()

    def _start_animation(self):
        if not self.animators:
            return
        self.is_animating = True
        now = QTime.currentTime().msecsSinceStartOfDay() / 1000.0
        for anim in self.animators:
            anim.reset(now)
        self.play_button.setText("⏸ Pause")
        self.animation_timer.start(200)  # ~5 fps (full scene rebuild each frame)

    def _stop_animation(self):
        self.is_animating = False
        self.animation_timer.stop()
        self.play_button.setText("▶  Play")

    def _animation_step(self):
        """Called by timer – advance every robot and re-render."""
        if not self.is_animating:
            return
        now = QTime.currentTime().msecsSinceStartOfDay() / 1000.0
        any_updated = False
        for anim in self.animators:
            if anim.step(now):
                any_updated = True
        if any_updated:
            self._redraw_trees()

    def _redraw_trees(self):
        """Remove old tree visuals and re-add them (keeps grid and labels)."""
        # Remove everything except grid and text labels
        items_to_remove = [item for item in self.gl_widget.items
                           if item is not self.grid
                           and not isinstance(item, gl.GLTextItem)]
        for item in items_to_remove:
            self.gl_widget.removeItem(item)
        # Re-add trees
        for t in self.trees:
            t.addToWidget(self.viewer,
                          showJointSurface=True,
                          showLinkSurface=True,
                          showLinkPath=False,
                          showJointPoses=False,
                          showSpheres=False)

    # ── Camera ───────────────────────────────────────────────────────────

    def _fit_camera(self):
        """Set camera distance and center so all robots are visible."""
        if not self.trees:
            return
        # Compute a bounding ball of all bounding balls
        all_centers = np.array([t.boundingBall.c for t in self.trees])
        all_radii = np.array([t.boundingBall.r for t in self.trees])

        # Simple encompassing sphere: center at mean, radius to farthest edge
        center = np.mean(all_centers, axis=0)
        max_dist = max(np.linalg.norm(c - center) + r
                       for c, r in zip(all_centers, all_radii))

        self.gl_widget.opts["center"] = pg.Vector(center[0], center[1], center[2])
        self.gl_widget.opts["distance"] = max_dist * 3.0
        self.gl_widget.opts["elevation"] = 25
        self.gl_widget.opts["azimuth"] = 45


# ── Entry point ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Dance Party Viewer")
    parser.add_argument("--env", choices=ENVIRONMENTS, default=None,
                        help="Environment to load (ballroom, beach, or club)")
    args = parser.parse_args()

    app = QApplication(sys.argv)
    window = DancePartyWindow(env=args.env)
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
