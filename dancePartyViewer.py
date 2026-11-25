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

import sys, os, glob, copy, argparse, math, time
import numpy as np
import dill

# Ensure the project root is on the path so local imports work
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QWidget,
                              QComboBox, QLabel, QHBoxLayout, QPushButton, QSpinBox,
                              QFileDialog)
from PyQt5.QtCore import Qt, QTimer, QTime
from PyQt5.QtGui import QSurfaceFormat, QVector3D, QImage
import json
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from OpenGL.GL import (glDisable, glEnable, GL_DEPTH_TEST,
                        glClearColor, glClear, GL_COLOR_BUFFER_BIT, GL_DEPTH_BUFFER_BIT,
                        glGenTextures, glDeleteTextures, glBindTexture,
                        glTexImage2D, glTexParameteri,
                        GL_TEXTURE_2D, GL_RGBA, GL_UNSIGNED_BYTE, GL_LINEAR,
                        GL_TEXTURE_MIN_FILTER, GL_TEXTURE_MAG_FILTER,
                        glMatrixMode, GL_PROJECTION, GL_MODELVIEW,
                        glPushMatrix, glPopMatrix, glLoadIdentity,
                        glOrtho, glBegin, glEnd, glVertex2f, glTexCoord2f,
                        GL_QUADS, glColor4f, glViewport)

from spatialmath import SE3
from KinematicTree import KinematicTree
from Joint import Prismatic, Revolute, Waypoint, Tip
from IntersectionHelper import (compute_cylinder_intersection,
                                 compute_torus_intersection,
                                 compute_closest_point_on_axis,
                                 compute_plane_intersection,
                                 compute_sphere_intersection)
from style import *

# ── Helpers ──────────────────────────────────────────────────────────────────

ENVIRONMENTS = ["ballroom", "beach", "club"]

# Colour palette for per-tree joint / link colouring
COLOR_PALETTE = {
    "Blue":        (0.0, 0.0, 1.0, 1.0),
    "Red":         (1.0, 0.0, 0.0, 1.0),
    "Green":       (0.0, 0.7, 0.0, 1.0),
    "Orange":      (1.0, 0.5, 0.0, 1.0),
    "Purple":      (0.6, 0.0, 0.8, 1.0),
    "Cyan":        (0.0, 0.8, 0.8, 1.0),
    "Yellow":      (0.9, 0.9, 0.0, 1.0),
    "Pink":        (1.0, 0.4, 0.7, 1.0),
    "Dark Gray":   (0.3, 0.3, 0.3, 1.0),
    "Light Gray":  (0.7, 0.7, 0.7, 1.0),
    "White":       (1.0, 1.0, 1.0, 1.0),
    "Black":       (0.0, 0.0, 0.0, 1.0),
    "Sky Blue":    (0.4, 0.7, 1.0, 1.0),
    "Lime":        (0.5, 1.0, 0.2, 1.0),
    "Coral":       (1.0, 0.5, 0.5, 1.0),
    "Lavender":    (0.7, 0.5, 1.0, 1.0),
    "Mint":        (0.4, 1.0, 0.7, 1.0),
    "Gold":        (1.0, 0.84, 0.0, 1.0),
    "Hot Pink":    (1.0, 0.2, 0.6, 1.0),
    "Peach":       (1.0, 0.8, 0.6, 1.0),
    "Aqua":        (0.0, 1.0, 1.0, 1.0),
    "Salmon":      (1.0, 0.6, 0.4, 1.0),
}
COLOR_NAMES = list(COLOR_PALETTE.keys())
DEFAULT_JOINT_COLOR = "Blue"
DEFAULT_LINK_COLOR = "Dark Gray"

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


# ── Overlay line (drawn on top of depth buffer) ─────────────────────────────

class OverlayLine(gl.GLLinePlotItem):
    """A GLLinePlotItem that renders on top of everything (no depth test)."""
    def paint(self):
        glDisable(GL_DEPTH_TEST)
        super().paint()
        glEnable(GL_DEPTH_TEST)


# ── Interactive GL widget with per-tree translate / rotate gizmos ────────────

class DancePartyGLWidget(gl.GLViewWidget):
    """GLViewWidget subclass that lets the user click to select a robot,
    then drag translation arrows or rotation rings to reposition it."""

    GIZMO_ARROW_PX = 80     # arrow length in screen pixels
    GIZMO_THICK_PX = 10     # hit-test thickness in screen pixels
    AXIS_COLORS = [(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)]
    SELECTED_COLOR = (1, 1, 0, 1)

    def __init__(self, parent_window, **kwargs):
        super().__init__(**kwargs)
        self.parent_window = parent_window
        self.orbit_speed = 0.3

        # Interaction state
        self.selected_tree_index = -1   # index into parent_window.trees
        self.control_mode = "Translate"  # "Translate" or "Rotate"
        self._is_dragging = False
        self._drag_start_pos = None
        self._selected_axis = None       # QVector3D – translation axis
        self._selected_torus = None      # QVector3D – rotation axis normal
        self._drag_prev_vector = None
        self._facing_same_dir = False
        self._gizmo_items: list = []     # overlay GL items for the gizmo

        # Background image state
        self._bg_texture_id = None
        self._pending_bg_path = ...  # sentinel: no pending change

    # ── Background image ─────────────────────────────────────────────────

    def set_background_image(self, image_path):
        """Queue an image to be uploaded as a GL texture on the next paint.
        Pass None to remove the background image."""
        self._pending_bg_path = image_path
        self.update()

    def _upload_bg_texture(self, image_path):
        """Upload *image_path* as an OpenGL texture (called inside paintGL)."""
        # Delete old texture
        if self._bg_texture_id is not None:
            glDeleteTextures([self._bg_texture_id])
            self._bg_texture_id = None
        if image_path is None:
            return
        img = QImage(image_path)
        if img.isNull():
            return
        img = img.convertToFormat(QImage.Format_RGBA8888).mirrored()
        w, h = img.width(), img.height()
        ptr = img.bits()
        ptr.setsize(w * h * 4)
        data = bytes(ptr)
        tex_id = glGenTextures(1)
        glBindTexture(GL_TEXTURE_2D, tex_id)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0,
                     GL_RGBA, GL_UNSIGNED_BYTE, data)
        glBindTexture(GL_TEXTURE_2D, 0)
        self._bg_texture_id = tex_id

    def _render_bg_quad(self):
        """Draw a fullscreen textured quad with the background image."""
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        glOrtho(0, 1, 0, 1, -1, 1)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        glEnable(GL_TEXTURE_2D)
        glBindTexture(GL_TEXTURE_2D, self._bg_texture_id)
        glColor4f(1, 1, 1, 1)
        glBegin(GL_QUADS)
        glTexCoord2f(0, 0); glVertex2f(0, 0)
        glTexCoord2f(1, 0); glVertex2f(1, 0)
        glTexCoord2f(1, 1); glVertex2f(1, 1)
        glTexCoord2f(0, 1); glVertex2f(0, 1)
        glEnd()
        glDisable(GL_TEXTURE_2D)
        glBindTexture(GL_TEXTURE_2D, 0)
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()

    def paintGL(self, region=None, viewport=None, useItemNames=False):
        # Process any pending background-image change
        if self._pending_bg_path is not ...:
            self._upload_bg_texture(self._pending_bg_path)
            self._pending_bg_path = ...

        # If no background image, use default pipeline
        if self._bg_texture_id is None:
            super().paintGL(region, viewport, useItemNames)
            return

        # Replicate pyqtgraph's paintGL with background image inserted
        if hasattr(self, 'prepareForPaint'):
            self.prepareForPaint()
        if viewport is not None:
            glViewport(*viewport)
        else:
            self.setProjection(region=region)
        self.setModelview()
        bg = self.opts.get('backgroundColor', self.opts.get('bgcolor', (0, 0, 0, 1)))
        if hasattr(bg, 'getRgbF'):
            bgcolor = bg.getRgbF()
        else:
            bgcolor = bg
        glClearColor(*bgcolor)
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT)
        self._render_bg_quad()
        glClear(GL_DEPTH_BUFFER_BIT)  # clear depth so 3D items render on top
        self.drawItemTree(useItemNames=useItemNames)

    # ── Gizmo drawing ────────────────────────────────────────────────────

    def _gizmo_center(self) -> np.ndarray | None:
        """Return the world-space center of the selected tree's bounding ball."""
        if self.selected_tree_index < 0:
            return None
        trees = self.parent_window.trees
        if self.selected_tree_index >= len(trees):
            return None
        return trees[self.selected_tree_index].boundingBall.c.copy()

    def draw_gizmo(self):
        """Draw translation arrows or rotation rings at the selected tree's center."""
        self.clear_gizmo()
        center = self._gizmo_center()
        if center is None:
            return
        rad = self._world_len(self.GIZMO_ARROW_PX)
        axes = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])]

        if self.control_mode == "Translate":
            for i, a in enumerate(axes):
                col = self.AXIS_COLORS[i]
                pts = np.array([center, center + rad * a])
                item = OverlayLine(pos=pts, color=col, width=8, antialias=True)
                self.addItem(item)
                self._gizmo_items.append(item)
        else:  # Rotate
            thick = self._world_len(self.GIZMO_THICK_PX)
            for i, axis in enumerate(axes):
                helper = np.array([1.0, 0.0, 0.0])
                if abs(np.dot(axis, helper)) > 0.9:
                    helper = np.array([0.0, 1.0, 0.0])
                u = np.cross(axis, helper)
                u /= np.linalg.norm(u)
                v = np.cross(axis, u)
                pts = np.array([
                    center + (rad + thick) * (u * np.cos(t) + v * np.sin(t))
                    for t in np.linspace(0, 2 * math.pi, 64)
                ])
                col = self.AXIS_COLORS[i]
                item = OverlayLine(pos=pts, color=col, width=8, antialias=True)
                self.addItem(item)
                self._gizmo_items.append(item)

    def clear_gizmo(self):
        for item in self._gizmo_items:
            try:
                self.removeItem(item)
            except ValueError:
                pass
        self._gizmo_items.clear()

    # ── Raycasting helpers ───────────────────────────────────────────────

    def _world_len(self, px: float) -> float:
        dist = self.opts['distance']
        fov = math.radians(self.opts.get('fov', 60))
        h = self.height()
        if h == 0:
            return 1.0
        return dist * math.tan(fov * (px / h))

    def _get_ray(self, event):
        """Return (origin: QVector3D, direction: QVector3D) for a mouse event."""
        pos = event.localPos() if hasattr(event, 'localPos') else event.position()
        ndc_x = (2.0 * pos.x()) / self.width() - 1.0
        ndc_y = 1.0 - (2.0 * pos.y()) / self.height()
        view = self.viewMatrix()
        proj = self.projectionMatrix()
        inv = (proj * view).inverted()[0]
        camera_pos = self.cameraPosition()
        world_pt = inv.map(QVector3D(ndc_x, ndc_y, 0.0))
        direction = world_pt - camera_pos
        direction.normalize()
        return camera_pos, direction

    def _hit_test_trees(self, origin, direction) -> int:
        """Return index of the closest tree whose root joint sphere (5× enlarged) is hit, or -1."""
        best_dist = float('inf')
        best_idx = -1
        for i, tree in enumerate(self.parent_window.trees):
            if not tree.Joints:
                continue
            root_ball = tree.Joints[0].boundingBall()
            c = root_ball.c
            r = root_ball.r * 5.0
            d = compute_sphere_intersection(
                [origin.x(), origin.y(), origin.z()],
                [direction.x(), direction.y(), direction.z()],
                c, r)
            if d < best_dist:
                best_dist = d
                best_idx = i
        return best_idx

    def _hit_test_translate(self, origin, direction, center_q):
        """Test ray against the three translation arrows. Returns axis index or -1."""
        rad = self._world_len(self.GIZMO_ARROW_PX)
        thick = self._world_len(self.GIZMO_THICK_PX)
        axes = [QVector3D(1, 0, 0), QVector3D(0, 1, 0), QVector3D(0, 0, 1)]
        best = float('inf')
        best_idx = -1
        for i, ax in enumerate(axes):
            d = compute_cylinder_intersection(origin, direction, center_q, ax, thick, rad)
            if d < best:
                best = d
                best_idx = i
        if best_idx >= 0:
            return best_idx, axes[best_idx]
        return -1, None

    def _hit_test_rotate(self, origin, direction, center_q):
        """Test ray against the three rotation tori. Returns axis index or -1."""
        rad = self._world_len(self.GIZMO_ARROW_PX)
        thick = self._world_len(self.GIZMO_THICK_PX)
        axes = [QVector3D(1, 0, 0), QVector3D(0, 1, 0), QVector3D(0, 0, 1)]
        best = float('inf')
        best_idx = -1
        for i, ax in enumerate(axes):
            d = compute_torus_intersection(origin, direction, center_q, ax,
                                           major_radius=rad, minor_radius=thick)
            if d < best:
                best = d
                best_idx = i
        if best_idx >= 0:
            return best_idx, axes[best_idx]
        return -1, None

    # ── Mouse events ─────────────────────────────────────────────────────

    def mousePressEvent(self, event):
        self._drag_start_pos = event.pos()

        if event.button() == Qt.LeftButton:
            origin, direction = self._get_ray(event)
            self._is_dragging = False
            self._selected_axis = None
            self._selected_torus = None

            # If a tree is selected, test gizmo hit first
            center = self._gizmo_center()
            if center is not None:
                center_q = QVector3D(center[0], center[1], center[2])
                if self.control_mode == "Translate":
                    idx, axis = self._hit_test_translate(origin, direction, center_q)
                    if idx >= 0:
                        self._selected_axis = axis
                        self._is_dragging = True
                        return
                else:
                    idx, axis = self._hit_test_rotate(origin, direction, center_q)
                    if idx >= 0:
                        self._selected_torus = axis
                        self._is_dragging = True
                        # Initialise rotation tracking
                        npos = compute_plane_intersection(origin, direction, axis, center_q)
                        if npos is not None:
                            pv = npos - center_q
                            pv.normalize()
                            self._drag_prev_vector = pv
                        dot = QVector3D.dotProduct(direction, axis)
                        self._facing_same_dir = (dot > 0)
                        return

            # No gizmo hit — test for tree selection
            hit = self._hit_test_trees(origin, direction)
            if hit != self.selected_tree_index:
                self.selected_tree_index = hit
                self.draw_gizmo()
                self.parent_window._sync_color_combos()
                self.update()

        # Fall through to default camera controls for middle button / unhandled
        if event.button() != Qt.LeftButton or not self._is_dragging:
            super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._is_dragging and self._selected_axis:
            # Translation drag
            origin, direction = self._get_ray(event)
            center = self._gizmo_center()
            if center is None:
                return
            center_q = QVector3D(center[0], center[1], center[2])
            new_pt = compute_closest_point_on_axis(origin, direction, center_q, self._selected_axis)
            delta = new_pt - center_q
            T = SE3.Trans(delta.x(), delta.y(), delta.z())
            tree = self.parent_window.trees[self.selected_tree_index]
            tree.transformAll(T)
            self.parent_window._redraw_trees()
            return

        if self._is_dragging and self._selected_torus:
            # Rotation drag
            origin, direction = self._get_ray(event)
            center = self._gizmo_center()
            if center is None:
                return
            center_q = QVector3D(center[0], center[1], center[2])
            npos = compute_plane_intersection(origin, direction, self._selected_torus, center_q)
            if npos is None:
                return
            plane_vector = npos - center_q
            plane_vector.normalize()
            prev = self._drag_prev_vector
            if prev is None:
                self._drag_prev_vector = plane_vector
                return
            d = QVector3D.dotProduct(plane_vector, prev)
            prod = plane_vector.length() * prev.length()
            if prod == 0:
                return
            d = max(-1.0, min(1.0, d / prod))
            angle = math.acos(d)
            cross = QVector3D.crossProduct(prev, plane_vector)
            normal = QVector3D.dotProduct(direction, self._selected_torus) * self._selected_torus
            normal.normalize()
            if QVector3D.dotProduct(cross, normal) < 0:
                angle = -angle
            if not self._facing_same_dir:
                angle = -angle
            self._drag_prev_vector = plane_vector

            # Build rotation about the tree center
            axis_np = np.array([self._selected_torus.x(),
                                self._selected_torus.y(),
                                self._selected_torus.z()])
            center_np = self._gizmo_center()
            T_to_origin = SE3.Trans(*(-center_np))
            T_rotate = SE3.AngleAxis(math.degrees(angle), axis_np, unit='deg')
            T_back = SE3.Trans(*center_np)
            T = T_back @ T_rotate @ T_to_origin

            tree = self.parent_window.trees[self.selected_tree_index]
            tree.transformAll(T)
            self.parent_window._redraw_trees()
            return

        # Default camera orbit / pan
        if event.buttons() == Qt.LeftButton and self._drag_start_pos is not None:
            curr = event.position() if hasattr(event, 'position') else event.localPos()
            prev = self._last_drag_pos if self._last_drag_pos is not None else self._drag_start_pos
            diff = curr - prev
            self._last_drag_pos = curr
            if event.modifiers() & Qt.ShiftModifier:
                self.pan(diff.x(), diff.y(), 0, relative='view')
            else:
                self.orbit(-diff.x() * self.orbit_speed, diff.y() * self.orbit_speed)
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        was_dragging_gizmo = self._is_dragging and (self._selected_axis or self._selected_torus)
        self._is_dragging = False
        self._selected_axis = None
        self._selected_torus = None
        self._drag_prev_vector = None
        self._last_drag_pos = None
        if was_dragging_gizmo:
            # Redraw gizmo at new position
            self.draw_gizmo()
            self.update()
        else:
            super().mouseReleaseEvent(event)

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_T:
            self.parent_window._set_translate_mode()
        elif event.key() == Qt.Key_R:
            self.parent_window._set_rotate_mode()
        elif event.key() == Qt.Key_Escape:
            self.selected_tree_index = -1
            self.clear_gizmo()
            self.parent_window._sync_color_combos()
            self.update()
        else:
            super().keyPressEvent(event)


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
        # Forced-loop state
        self._original_configs = list(saved_configs)
        self._original_durations = list(config_durations)
        self._forced_loop = False

    def naturally_loops(self, tol: float = 1e-6) -> bool:
        """Return True if the first and last configs are (nearly) identical."""
        if len(self.saved_configs) < 2:
            return True
        first, last = self.saved_configs[0], self.saved_configs[-1]
        return all(abs(a - b) < tol for a, b in zip(first, last))

    def set_forced_loop(self, enabled: bool):
        """When enabled, append a reversed copy of the dance so non-looping
        dances return smoothly to the start configuration."""
        if enabled == self._forced_loop:
            return
        self._forced_loop = enabled
        if enabled and not self.naturally_loops():
            # Append reversed configs (skip the duplicate endpoints)
            reversed_configs = list(reversed(self._original_configs[:-1]))
            reversed_durations = list(reversed(self._original_durations))
            self.saved_configs = self._original_configs + reversed_configs
            self.config_durations = self._original_durations + reversed_durations
        else:
            self.saved_configs = list(self._original_configs)
            self.config_durations = list(self._original_durations)

    @property
    def num_segments(self):
        return max(0, len(self.saved_configs) - 1)

    def can_animate(self):
        return self.num_segments >= 1

    def total_duration(self) -> float:
        """Total time (seconds) for one complete loop."""
        return sum(self._segment_duration(i) for i in range(self.num_segments))

    def set_time(self, elapsed: float, lightweight: bool = False):
        """Set the robot to the pose at *elapsed* seconds from the start, looping."""
        if not self.can_animate():
            return
        total = self.total_duration()
        if total <= 0:
            return
        elapsed = elapsed % total  # loop
        t_acc = 0.0
        for seg_idx in range(self.num_segments):
            dur = self._segment_duration(seg_idx)
            if t_acc + dur > elapsed:
                t = (elapsed - t_acc) / dur
                self._interpolate(seg_idx, t, lightweight=lightweight)
                return
            t_acc += dur
        self._interpolate(self.num_segments - 1, 1.0, lightweight=lightweight)  # very end

    def reset(self, now: float):
        self.current_segment = 0
        self.segment_start_time = now

    def step(self, now: float, lightweight: bool = False):
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

        self._interpolate(self.current_segment, t, lightweight=lightweight)
        return True

    # ── internals ──

    def _segment_duration(self, seg_idx: int) -> float:
        if seg_idx < len(self.config_durations):
            return max(self.config_durations[seg_idx], 0.01)
        return 1.0  # fallback

    def _interpolate(self, seg_idx: int, t: float, lightweight: bool = False):
        """Set tree joint states to interpolation between config seg_idx and seg_idx+1.
        When lightweight=True, skip link rebuilds and collision (for animation fast-path)."""
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
            self.tree.setJointState(joint_idx, val, lightweight=lightweight)


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

        # Export Video button
        self.export_button = QPushButton("Export Video")
        self.export_button.setFixedWidth(100)
        self.export_button.clicked.connect(self._export_video)
        selector_row.addWidget(self.export_button)

        # Debug Export button (low-res, 5 frames)
        self.debug_export_button = QPushButton("Debug Export")
        self.debug_export_button.setFixedWidth(100)
        self.debug_export_button.clicked.connect(lambda: self._export_video(debug=True))
        selector_row.addWidget(self.debug_export_button)

        # Show Names checkbox
        from PyQt5.QtWidgets import QCheckBox
        self.axes_checkbox = QCheckBox("Axes")
        self.axes_checkbox.setChecked(False)
        self.axes_checkbox.toggled.connect(lambda: self._redraw_trees())
        selector_row.addWidget(self.axes_checkbox)

        self.names_checkbox = QCheckBox("Names")
        self.names_checkbox.setChecked(False)
        self.names_checkbox.toggled.connect(self._toggle_names)
        selector_row.addWidget(self.names_checkbox)

        self.loop_checkbox = QCheckBox("Loop All")
        self.loop_checkbox.setChecked(True)
        self.loop_checkbox.toggled.connect(self._toggle_loop)
        selector_row.addWidget(self.loop_checkbox)

        selector_row.addWidget(QLabel("Export Iterations:"))
        self.export_iterations_spin = QSpinBox()
        self.export_iterations_spin.setMinimum(1)
        self.export_iterations_spin.setMaximum(999)
        self.export_iterations_spin.setValue(1)
        self.export_iterations_spin.setFixedWidth(60)
        selector_row.addWidget(self.export_iterations_spin)

        # Background selector
        selector_row.addWidget(QLabel("Background:"))
        self.bg_combo = QComboBox()
        self.bg_combo.addItems(["White", "Green Screen", "Environment Image"])
        self.bg_combo.currentTextChanged.connect(self._on_bg_changed)
        selector_row.addWidget(self.bg_combo)

        # Translate / Rotate mode buttons
        self.translate_button = QPushButton("Translate (T)")
        self.translate_button.setFixedWidth(100)
        self.translate_button.setCheckable(True)
        self.translate_button.setChecked(True)
        self.translate_button.clicked.connect(self._set_translate_mode)
        selector_row.addWidget(self.translate_button)

        self.rotate_button = QPushButton("Rotate (R)")
        self.rotate_button.setFixedWidth(100)
        self.rotate_button.setCheckable(True)
        self.rotate_button.clicked.connect(self._set_rotate_mode)
        selector_row.addWidget(self.rotate_button)

        # Joint / Link colour selectors (apply to the selected tree)
        selector_row.addWidget(QLabel("Joint:"))
        self.joint_color_combo = QComboBox()
        self.joint_color_combo.addItems(COLOR_NAMES)
        self.joint_color_combo.setCurrentText(DEFAULT_JOINT_COLOR)
        self.joint_color_combo.currentTextChanged.connect(self._on_joint_color_changed)
        selector_row.addWidget(self.joint_color_combo)

        selector_row.addWidget(QLabel("Link:"))
        self.link_color_combo = QComboBox()
        self.link_color_combo.addItems(COLOR_NAMES)
        self.link_color_combo.setCurrentText(DEFAULT_LINK_COLOR)
        self.link_color_combo.currentTextChanged.connect(self._on_link_color_changed)
        selector_row.addWidget(self.link_color_combo)

        # Save / Load scene buttons
        self.save_scene_button = QPushButton("Save Scene")
        self.save_scene_button.setFixedWidth(90)
        self.save_scene_button.clicked.connect(self._save_scene)
        selector_row.addWidget(self.save_scene_button)

        self.load_scene_button = QPushButton("Load Scene")
        self.load_scene_button.setFixedWidth(90)
        self.load_scene_button.clicked.connect(self._load_scene)
        selector_row.addWidget(self.load_scene_button)

        # Profile button
        self.profile_button = QPushButton("Profile")
        self.profile_button.setFixedWidth(70)
        self.profile_button.clicked.connect(self._profile_redraw)
        selector_row.addWidget(self.profile_button)

        layout.addLayout(selector_row, 0)  # stretch=0: selector takes minimum height

        # ── GL view ──
        fmt = QSurfaceFormat()
        fmt.setDepthBufferSize(24)
        fmt.setProfile(QSurfaceFormat.CompatibilityProfile)
        QSurfaceFormat.setDefaultFormat(fmt)

        self.gl_widget = DancePartyGLWidget(parent_window=self)
        self.gl_widget.setBackgroundColor(backgroundColorDefault)
        layout.addWidget(self.gl_widget, 1)  # stretch=1: GL view fills remaining space

        # The wrapper object that tree.addToWidget expects
        self.viewer = SimpleTreeViewer(self.gl_widget)

        # Storage
        self.trees: list[KinematicTree] = []
        self.tree_names: list[str] = []
        self.tree_joint_colors: list[str] = []   # per-tree joint colour name
        self.tree_link_colors: list[str] = []    # per-tree link colour name
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

        # Render trees into the GL widget
        show_axes = self.axes_checkbox.isChecked()
        for i, t in enumerate(trees):
            jc = COLOR_PALETTE[DEFAULT_JOINT_COLOR]
            lc = COLOR_PALETTE[DEFAULT_LINK_COLOR]
            t.addToWidget(self.viewer,
                          showJointSurface=True,
                          showLinkSurface=True,
                          showLinkPath=False,
                          showJointPoses=False,
                          showJointAxis=show_axes,
                          showSpheres=False,
                          jointColor=jc,
                          linkColor=lc,
                          linkOpacity=lc[3])

        self.trees = trees
        self.tree_names = names
        self.tree_joint_colors = [DEFAULT_JOINT_COLOR] * len(trees)
        self.tree_link_colors = [DEFAULT_LINK_COLOR] * len(trees)
        self.name_labels: list[gl.GLTextItem] = []

        # Create per-robot animators
        self.animators = []
        anim_count = 0
        for t, cfgs, durs in zip(trees, session_configs, session_durations):
            anim = RobotAnimator(t, cfgs, durs)
            anim.set_forced_loop(self.loop_checkbox.isChecked())
            self.animators.append(anim)
            if anim.can_animate():
                anim_count += 1

        # Set camera to see everything
        self._fit_camera()
        # Refresh background in case "Environment Image" is selected
        self._on_bg_changed()
        self.info_label.setText(
            f"Loaded {len(trees)} robots for '{env}' "
            f"({anim_count} with dance sequences)"
        )

    def _clear_scene(self):
        """Remove all items from the GL widget."""
        items = list(self.gl_widget.items)
        for item in items:
            self.gl_widget.removeItem(item)
        self.trees = []
        self.tree_names = []
        self.tree_joint_colors = []
        self.tree_link_colors = []
        self.name_labels = []
        self.animators = []

    # ── Background ────────────────────────────────────────────────────

    def _on_bg_changed(self, text=None):
        if text is None:
            text = self.bg_combo.currentText()
        if text == "White":
            self.gl_widget.set_background_image(None)
            self.gl_widget.setBackgroundColor((255, 255, 255, 255))
        elif text == "Green Screen":
            self.gl_widget.set_background_image(None)
            self.gl_widget.setBackgroundColor((0, 177, 64, 255))
        elif text == "Environment Image":
            env = self.env_combo.currentText()
            img_path = os.path.join(session_dir(), f"{env}.jpg")
            if os.path.isfile(img_path):
                self.gl_widget.set_background_image(img_path)
            else:
                self.gl_widget.set_background_image(None)
                self.gl_widget.setBackgroundColor((255, 255, 255, 255))
        self.gl_widget.update()

    # ── Per-tree colours ──────────────────────────────────────────────

    def _sync_color_combos(self):
        """Update the colour combo boxes to reflect the selected tree."""
        idx = self.gl_widget.selected_tree_index
        if 0 <= idx < len(self.trees):
            self.joint_color_combo.blockSignals(True)
            self.joint_color_combo.setCurrentText(self.tree_joint_colors[idx])
            self.joint_color_combo.blockSignals(False)
            self.link_color_combo.blockSignals(True)
            self.link_color_combo.setCurrentText(self.tree_link_colors[idx])
            self.link_color_combo.blockSignals(False)

    def _on_joint_color_changed(self, name: str):
        idx = self.gl_widget.selected_tree_index
        if idx < 0 or idx >= len(self.trees):
            return
        self.tree_joint_colors[idx] = name
        self._redraw_trees()

    def _on_link_color_changed(self, name: str):
        idx = self.gl_widget.selected_tree_index
        if idx < 0 or idx >= len(self.trees):
            return
        self.tree_link_colors[idx] = name
        self._redraw_trees()

    # ── Transform modes ────────────────────────────────────────────────

    def _set_translate_mode(self):
        self.gl_widget.control_mode = "Translate"
        self.translate_button.setChecked(True)
        self.rotate_button.setChecked(False)
        self.gl_widget.draw_gizmo()
        self.info_label.setText("Mode: Translate")

    def _set_rotate_mode(self):
        self.gl_widget.control_mode = "Rotate"
        self.translate_button.setChecked(False)
        self.rotate_button.setChecked(True)
        self.gl_widget.draw_gizmo()
        self.info_label.setText("Mode: Rotate")

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
        self.animation_timer.start(80)  # ~12 fps

    def _stop_animation(self):
        self.is_animating = False
        self.animation_timer.stop()
        self.play_button.setText("▶  Play")
        # Resync link geometry after lightweight animation
        for t in self.trees:
            t.resyncFromLightweight()

    def _animation_step(self):
        """Called by timer – advance every robot and re-render."""
        if not self.is_animating:
            return
        now = QTime.currentTime().msecsSinceStartOfDay() / 1000.0
        any_updated = False
        for anim in self.animators:
            if anim.step(now, lightweight=True):
                any_updated = True
        if any_updated:
            self._fast_update_trees()

    def _fast_update_trees(self):
        """Try to update GL items via model-matrix transforms (fast path).
        Falls back to full redraw if any tree lacks a valid GL cache."""
        for t in self.trees:
            if not t.updateGLTransforms():
                # Cache miss — do a full redraw (which re-populates caches)
                self._redraw_trees()
                return
        self.gl_widget.draw_gizmo()
        self.gl_widget.update()

    def _redraw_trees(self):
        """Remove old tree visuals and re-add them (full rebuild)."""
        items_to_remove = list(self.gl_widget.items)
        for item in items_to_remove:
            self.gl_widget.removeItem(item)
        # Clear caches so they get rebuilt
        for t in self.trees:
            t.clearAllGLCaches()
        # Re-add trees
        show_axes = self.axes_checkbox.isChecked()
        for i, t in enumerate(self.trees):
            jc = COLOR_PALETTE.get(
                self.tree_joint_colors[i] if i < len(self.tree_joint_colors) else DEFAULT_JOINT_COLOR,
                COLOR_PALETTE[DEFAULT_JOINT_COLOR])
            lc = COLOR_PALETTE.get(
                self.tree_link_colors[i] if i < len(self.tree_link_colors) else DEFAULT_LINK_COLOR,
                COLOR_PALETTE[DEFAULT_LINK_COLOR])
            t.addToWidget(self.viewer,
                          showJointSurface=True,
                          showLinkSurface=True,
                          showLinkPath=False,
                          showJointPoses=False,
                          showJointAxis=show_axes,
                          showSpheres=False,
                          jointColor=jc,
                          linkColor=lc,
                          linkOpacity=lc[3])
        # Re-add name labels if enabled
        if self.names_checkbox.isChecked():
            self._add_name_labels()
        # Re-add gizmo overlay if a tree is selected
        self.gl_widget.draw_gizmo()

    def _add_name_labels(self):
        """Create and add GLTextItem labels above each robot."""
        self._remove_name_labels()
        for t, name in zip(self.trees, self.tree_names):
            label_pos = t.boundingBall.c.copy()
            label_pos[2] += t.boundingBall.r + 2.0
            text_item = gl.GLTextItem(pos=label_pos, text=name, color=(0, 0, 0, 255))
            self.gl_widget.addItem(text_item)
            self.name_labels.append(text_item)

    def _remove_name_labels(self):
        """Remove any existing name label items from the scene."""
        for item in self.name_labels:
            try:
                self.gl_widget.removeItem(item)
            except ValueError:
                pass
        self.name_labels.clear()

    def _toggle_names(self, checked: bool):
        """Show or hide robot name labels."""
        if checked:
            self._add_name_labels()
        else:
            self._remove_name_labels()

    def _toggle_loop(self, checked: bool):
        """Enable or disable forced looping for non-looping dances."""
        for anim in self.animators:
            anim.set_forced_loop(checked)

    # ── Video export ─────────────────────────────────────────────────────

    def _grab_frame(self, target_width: int = 0, target_height: int = 0) -> np.ndarray:
        """Capture the GL widget's current framebuffer and return as an RGB
        numpy array.

        If target_width/target_height are non-zero the image is scaled;
        otherwise the native framebuffer resolution is used (pixel-perfect).
        """
        qimg = self.gl_widget.grabFramebuffer()          # exact on-screen content
        if target_width > 0 and target_height > 0:
            qimg = qimg.scaled(target_width, target_height,
                               Qt.IgnoreAspectRatio, Qt.SmoothTransformation)
        qimg = qimg.convertToFormat(qimg.Format_RGB888)
        w, h = qimg.width(), qimg.height()
        ptr = qimg.constBits()
        ptr.setsize(h * w * 3)
        arr = np.array(ptr, dtype=np.uint8).reshape(h, w, 3).copy()
        return arr

    def _export_video(self, *, debug: bool = False):
        """Render the dance animation offline and save as an MP4 file.

        Parameters
        ----------
        debug : bool
            If True, render only 5 frames at 320×240 for fast iteration.
        """
        try:
            import imageio
        except ImportError:
            QtWidgets.QMessageBox.warning(
                self, "Missing dependency",
                "Please install imageio and imageio-ffmpeg:\n"
                "  pip install imageio imageio-ffmpeg")
            return

        animatable = [a for a in self.animators if a.can_animate()]
        if not self.trees or not animatable:
            QtWidgets.QMessageBox.information(self, "Nothing to export",
                                              "Load an environment with dance sequences first.")
            return

        # Stop live playback
        was_animating = self.is_animating
        if was_animating:
            self._stop_animation()

        # Choose save path
        default_name = f"dance_{self.env_combo.currentText()}.mp4"
        default_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     default_name)
        filepath, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export Video", default_path, "MP4 Video (*.mp4)")
        if not filepath:
            return

        # Parameters
        fps = 30

        if debug:
            width, height = 320, 240
            num_frames = 5
            export_iterations = 1
            longest_loop_duration = 0.0
        else:
            # Use native framebuffer resolution (pixel-perfect, no scaling).
            # Grab one test frame to get the actual pixel dimensions.
            test_img = self._grab_frame()
            height, width = test_img.shape[:2]
            # H.264 prefers even dimensions; crop by 1px if needed.
            if width % 2 != 0:
                width -= 1
            if height % 2 != 0:
                height -= 1
            export_iterations = self.export_iterations_spin.value()
            longest_loop_duration = max(a.total_duration() for a in animatable)
            num_frames = int(math.ceil(longest_loop_duration * export_iterations * fps))

        # Progress dialog
        if debug:
            progress_text = f"Rendering {num_frames} frames at {width}×{height} …"
        else:
            progress_text = (
                f"Rendering {num_frames} frames at {width}×{height} "
                f"({export_iterations}× longest loop) …"
            )
        progress = QtWidgets.QProgressDialog(
            progress_text,
            "Cancel", 0, num_frames, self)
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)
        progress.show()
        QApplication.processEvents()

        writer = imageio.get_writer(filepath, fps=fps,
                                     codec='libx264',
                                     quality=9,
                                     output_params=['-pix_fmt', 'yuv420p'])
        t_start = time.perf_counter()
        try:
            for frame_idx in range(num_frames):
                if progress.wasCanceled():
                    break

                t = frame_idx / fps
                for anim in self.animators:
                    anim.set_time(t, lightweight=True)

                self._fast_update_trees()
                QApplication.processEvents()

                frame = self._grab_frame(width, height) if debug else self._grab_frame()
                # Crop to even dimensions if native size is odd
                frame = frame[:height, :width, :]
                writer.append_data(frame)

                progress.setValue(frame_idx + 1)
        finally:
            writer.close()
            progress.close()

        elapsed = time.perf_counter() - t_start
        self.info_label.setText(
            f"Exported {num_frames} frames in {elapsed:.1f}s → {os.path.basename(filepath)}")

        # Restore the scene to current time
        self._redraw_trees()

    def _profile_redraw(self):
        """Profile a single redraw cycle: measure time for each phase and
        report vertex/face counts per GL item type."""
        if not self.trees:
            print("No trees loaded.")
            return

        was_animating = self.is_animating
        if was_animating:
            self._stop_animation()

        print("\n" + "=" * 70)
        print("REDRAW PROFILE")
        print("=" * 70)

        # ── Phase 1: remove old items ──
        t0 = time.perf_counter()
        items_to_remove = list(self.gl_widget.items)
        num_removed = len(items_to_remove)
        for item in items_to_remove:
            self.gl_widget.removeItem(item)
        t_remove = time.perf_counter() - t0
        print(f"\nPhase 1 – Remove {num_removed} items: {t_remove*1000:.1f} ms")

        # ── Phase 2: addToWidget per tree ──
        tree_times = []
        for i, t in enumerate(self.trees):
            t1 = time.perf_counter()
            before_count = len(self.gl_widget.items)
            jc = COLOR_PALETTE.get(
                self.tree_joint_colors[i] if i < len(self.tree_joint_colors) else DEFAULT_JOINT_COLOR,
                COLOR_PALETTE[DEFAULT_JOINT_COLOR])
            lc = COLOR_PALETTE.get(
                self.tree_link_colors[i] if i < len(self.tree_link_colors) else DEFAULT_LINK_COLOR,
                COLOR_PALETTE[DEFAULT_LINK_COLOR])
            t.addToWidget(self.viewer,
                          showJointSurface=True,
                          showLinkSurface=True,
                          showLinkPath=False,
                          showJointPoses=False,
                          showJointAxis=False,
                          showSpheres=False,
                          jointColor=jc,
                          linkColor=lc,
                          linkOpacity=lc[3])
            after_count = len(self.gl_widget.items)
            dt = time.perf_counter() - t1
            tree_times.append(dt)
            name = self.tree_names[i] if i < len(self.tree_names) else f"tree_{i}"
            print(f"\nPhase 2 – addToWidget '{name}': {dt*1000:.1f} ms  "
                  f"({after_count - before_count} new items)")

        t_add_total = sum(tree_times)
        print(f"\nPhase 2 total: {t_add_total*1000:.1f} ms")

        # ── Phase 3: inventory all GL items ──
        print(f"\n{'Type':<30} {'Count':>6} {'Verts':>10} {'Faces':>10}")
        print("-" * 60)
        type_stats = {}  # type_name -> (count, total_verts, total_faces)
        for item in self.gl_widget.items:
            type_name = type(item).__name__
            obj_name = item.objectName() if hasattr(item, 'objectName') else ""
            key = f"{type_name}[{obj_name}]" if obj_name else type_name

            verts = 0
            faces = 0
            # Try to extract mesh data
            md = None
            if hasattr(item, 'opts'):
                md = item.opts.get('meshdata', None)
            if md is None and hasattr(item, 'meshDataChanged'):
                try:
                    item.meshDataChanged()
                    if hasattr(item, 'opts'):
                        md = item.opts.get('meshdata', None)
                except Exception:
                    pass
            if md is not None:
                try:
                    v = md.vertexes()
                    if v is not None:
                        verts = len(v)
                except Exception:
                    pass
                try:
                    f = md.faces()
                    if f is not None:
                        faces = len(f)
                except Exception:
                    pass
            # For line items, count points
            if hasattr(item, 'pos') and verts == 0:
                try:
                    p = item.pos
                    if isinstance(p, np.ndarray):
                        verts = len(p)
                except Exception:
                    pass

            if key not in type_stats:
                type_stats[key] = [0, 0, 0]
            type_stats[key][0] += 1
            type_stats[key][1] += verts
            type_stats[key][2] += faces

        # Sort by total faces descending (most geometry first)
        for key, (count, verts, faces) in sorted(type_stats.items(),
                                                  key=lambda x: -x[1][2]):
            print(f"{key:<30} {count:>6} {verts:>10} {faces:>10}")

        total_items = sum(s[0] for s in type_stats.values())
        total_verts = sum(s[1] for s in type_stats.values())
        total_faces = sum(s[2] for s in type_stats.values())
        print("-" * 60)
        print(f"{'TOTAL':<30} {total_items:>6} {total_verts:>10} {total_faces:>10}")
        print(f"\nFull redraw time: {(t_remove + t_add_total)*1000:.1f} ms")

        # ── Phase 4: detailed breakdown of Joint mesh items ──
        print(f"\n{'─' * 70}")
        print("JOINT MESH DETAIL  (individual GLMeshItem[Joint] items)")
        print(f"{'─' * 70}")
        print(f"{'Verts':>8} {'Faces':>8}  Description")
        print(f"{'─'*8} {'─'*8}  {'─'*40}")

        # Histogram: verts -> (count, faces)
        joint_histogram = {}
        joint_items_detail = []
        for item in self.gl_widget.items:
            obj_name = item.objectName() if hasattr(item, 'objectName') else ""
            if obj_name != "Joint":
                continue
            md = None
            if hasattr(item, 'opts'):
                md = item.opts.get('meshdata', None)
            verts = 0
            faces = 0
            if md is not None:
                try:
                    v = md.vertexes()
                    if v is not None:
                        verts = len(v)
                except Exception:
                    pass
                try:
                    f = md.faces()
                    if f is not None:
                        faces = len(f)
                except Exception:
                    pass
            # Guess what it is based on geometry
            if verts == 8 and faces == 12:
                desc = "Box (servo body)"
            elif faces == 2:
                desc = "Bracket quad"
            elif verts > 0 and faces > 0 and faces == verts - 2:
                desc = f"End cap (fan, {verts} pts)"
            elif verts > 0 and faces > 0:
                # Likely a cylinder
                # numCircles * numPointsPerCircle = verts
                # (numCircles-1) * numPointsPerCircle * 2 = faces
                # Solve: if faces = 2*(numCircles-1)*numPts and verts = numCircles*numPts
                # Try numPts = 32 (default)
                for npts in [32, 16, 10, 8]:
                    if verts % npts == 0:
                        nc = verts // npts
                        expected_faces = 2 * (nc - 1) * npts
                        if expected_faces == faces:
                            desc = f"Cylinder ({npts} pts/circle × {nc} circles)"
                            break
                else:
                    desc = f"Mesh ({verts}v, {faces}f)"
            else:
                desc = f"Unknown ({verts}v, {faces}f)"

            joint_items_detail.append((verts, faces, desc))
            key = (verts, faces, desc)
            joint_histogram[key] = joint_histogram.get(key, 0) + 1

        # Print histogram sorted by total faces contribution descending
        print(f"\n{'Count':>6} {'Verts/ea':>9} {'Faces/ea':>9} {'Tot Faces':>10}  Description")
        print(f"{'─'*6} {'─'*9} {'─'*9} {'─'*10}  {'─'*40}")
        for (v, f, desc), count in sorted(joint_histogram.items(),
                                           key=lambda x: -(x[0][1] * x[1])):
            print(f"{count:>6} {v:>9} {f:>9} {f*count:>10}  {desc}")

        tot_joint_items = len(joint_items_detail)
        tot_joint_verts = sum(v for v, _, _ in joint_items_detail)
        tot_joint_faces = sum(f for _, f, _ in joint_items_detail)
        print(f"{'─'*6} {'─'*9} {'─'*9} {'─'*10}  {'─'*40}")
        print(f"{tot_joint_items:>6} {'':>9} {'':>9} {tot_joint_faces:>10}  TOTAL Joint meshes ({tot_joint_verts} verts)")

        # ── Phase 5: profile fast path (updateGLTransforms) ──
        print(f"\n{'─' * 70}")
        print("FAST PATH PROFILE  (setTransform model-matrix update)")
        print(f"{'─' * 70}")
        # Advance animation by one step to change poses, then time the fast update
        if self.animators:
            # Time the animation interpolation (setJointState math)
            # Time lightweight interpolation
            t_interp_start = time.perf_counter()
            for anim in self.animators:
                anim.set_time(0.5, lightweight=True)  # arbitrary time to change poses
            t_interp = time.perf_counter() - t_interp_start
            print(f"Animation interpolation (lightweight setJointState): {t_interp*1000:.3f} ms")
            
            # Time just the GL transform update
            t_fast_start = time.perf_counter()
            fast_ok = True
            for ti, t_tree in enumerate(self.trees):
                if not t_tree.updateGLTransforms():
                    fast_ok = False
                    break
            t_gl = time.perf_counter() - t_fast_start
            if fast_ok:
                self.gl_widget.update()
            t_fast = t_interp + t_gl
            if fast_ok:
                print(f"GL updateGLTransforms: {t_gl*1000:.3f} ms")
                print(f"Total fast path: {t_fast*1000:.3f} ms")
                speedup = (t_remove + t_add_total) / t_fast if t_fast > 0 else float('inf')
                print(f"Speedup vs full redraw: {speedup:.0f}x")
            else:
                print("Fast path FAILED (cache miss) — would fall back to full redraw")
            # Reset to original time (lightweight to stay consistent)
            for anim in self.animators:
                anim.set_time(0.0, lightweight=True)
            self._fast_update_trees()

        print("=" * 70)

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

    # ── Save / Load scene ────────────────────────────────────────────────

    def _gather_scene_state(self) -> dict:
        """Collect all scene state into a serializable dict."""
        # Per-tree: capture the 4x4 pose of the root joint as the world transform
        tree_states = []
        for i, t in enumerate(self.trees):
            root_pose = t.Joints[0].Pose.A.tolist() if t.Joints else np.eye(4).tolist()
            tree_states.append({
                "root_pose": root_pose,
                "joint_color": self.tree_joint_colors[i] if i < len(self.tree_joint_colors) else DEFAULT_JOINT_COLOR,
                "link_color": self.tree_link_colors[i] if i < len(self.tree_link_colors) else DEFAULT_LINK_COLOR,
            })
        # Camera
        center = self.gl_widget.opts["center"]
        camera = {
            "center": [center.x(), center.y(), center.z()],
            "distance": float(self.gl_widget.opts["distance"]),
            "elevation": float(self.gl_widget.opts["elevation"]),
            "azimuth": float(self.gl_widget.opts["azimuth"]),
        }
        return {
            "environment": self.env_combo.currentText(),
            "background": self.bg_combo.currentText(),
            "trees": tree_states,
            "camera": camera,
        }

    def _save_scene(self):
        """Save the current scene arrangement to a JSON file."""
        if not self.trees:
            QtWidgets.QMessageBox.information(self, "Nothing to save",
                                              "Load an environment first.")
            return
        default_name = f"scene_{self.env_combo.currentText()}.dpscene"
        default_path = os.path.join(session_dir(), default_name)
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Save Scene", default_path,
            "Dance Party Scene (*.dpscene)")
        if not filepath:
            return
        state = self._gather_scene_state()
        with open(filepath, "w") as f:
            json.dump(state, f, indent=2)
        self.info_label.setText(f"Scene saved to {os.path.basename(filepath)}")

    def _load_scene(self):
        """Load a scene arrangement from a JSON file."""
        default_path = session_dir()
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Load Scene", default_path,
            "Dance Party Scene (*.dpscene)")
        if not filepath:
            return
        with open(filepath, "r") as f:
            state = json.load(f)

        env = state.get("environment", "")
        if env not in ENVIRONMENTS:
            QtWidgets.QMessageBox.warning(self, "Invalid scene",
                                          f"Unknown environment '{env}'.")
            return

        # Load the environment (resets trees to default layout)
        self.env_combo.blockSignals(True)
        self.env_combo.setCurrentText(env)
        self.env_combo.blockSignals(False)
        self._on_env_changed(env)

        saved_trees = state.get("trees", [])
        if len(saved_trees) != len(self.trees):
            QtWidgets.QMessageBox.warning(
                self, "Mismatch",
                f"Scene has {len(saved_trees)} trees but environment "
                f"loaded {len(self.trees)}. Applying what matches.")

        # Apply saved transforms and colors
        count = min(len(saved_trees), len(self.trees))
        for i in range(count):
            ts = saved_trees[i]
            # Compute delta from current root pose to saved root pose
            saved_pose = SE3(np.array(ts["root_pose"]))
            current_pose = self.trees[i].Joints[0].Pose if self.trees[i].Joints else SE3()
            delta = saved_pose * current_pose.inv()
            self.trees[i].transformAll(delta)
            self.trees[i].recomputeBoundingBall()
            # Colors
            if i < len(self.tree_joint_colors):
                self.tree_joint_colors[i] = ts.get("joint_color", DEFAULT_JOINT_COLOR)
            if i < len(self.tree_link_colors):
                self.tree_link_colors[i] = ts.get("link_color", DEFAULT_LINK_COLOR)

        # Restore camera
        cam = state.get("camera", {})
        if cam:
            c = cam.get("center", [0, 0, 0])
            self.gl_widget.opts["center"] = pg.Vector(c[0], c[1], c[2])
            self.gl_widget.opts["distance"] = cam.get("distance", 500)
            self.gl_widget.opts["elevation"] = cam.get("elevation", 25)
            self.gl_widget.opts["azimuth"] = cam.get("azimuth", 45)

        # Restore background
        bg = state.get("background", "White")
        self.bg_combo.setCurrentText(bg)

        # Sync color combos and redraw
        self._sync_color_combos()
        self._redraw_trees()
        self.info_label.setText(f"Scene loaded from {os.path.basename(filepath)}")


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
