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
                              QComboBox, QLabel, QHBoxLayout, QPushButton, QSpinBox)
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

        self.gl_widget = gl.GLViewWidget()
        self.gl_widget.setBackgroundColor(backgroundColorDefault)
        layout.addWidget(self.gl_widget, 1)  # stretch=1: GL view fills remaining space

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

        # Render trees into the GL widget
        show_axes = self.axes_checkbox.isChecked()
        for t in trees:
            t.addToWidget(self.viewer,
                          showJointSurface=True,
                          showLinkSurface=True,
                          showLinkPath=False,
                          showJointPoses=False,
                          showJointAxis=show_axes,
                          showSpheres=False)

        self.trees = trees
        self.tree_names = names
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
        self.name_labels = []
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
        for t in self.trees:
            t.addToWidget(self.viewer,
                          showJointSurface=True,
                          showLinkSurface=True,
                          showLinkPath=False,
                          showJointPoses=False,
                          showJointAxis=show_axes,
                          showSpheres=False)
        # Re-add name labels if enabled
        if self.names_checkbox.isChecked():
            self._add_name_labels()

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
            t.addToWidget(self.viewer,
                          showJointSurface=True,
                          showLinkSurface=True,
                          showLinkPath=False,
                          showJointPoses=False,
                          showJointAxis=False,
                          showSpheres=False)
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
