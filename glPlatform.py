"""
Chooses PyOpenGL's platform on Linux so it matches the OpenGL context Qt creates. Import this before anything
imports OpenGL (Joint.py imports pyqtgraph.opengl, for example), as the GUI scripts do on their first line.

From version 3.1.6, PyOpenGL on Linux picks EGL whenever WAYLAND_DISPLAY is set, as it is under WSLg and many
desktop sessions. But PyQt5 opens its windows through X11 (xcb) by default, which creates a GLX context, so
PyOpenGL then finds no current EGL context and every pyqtgraph item fails to draw ("Attempt to retrieve context
when no valid context"). This sets PYOPENGL_PLATFORM=glx unless it's already set or Qt is told to use Wayland
(QT_QPA_PLATFORM=wayland, in which case Qt creates an EGL context and PyOpenGL's guess is right).
Windows and macOS are unaffected.
"""
import os
import sys

if sys.platform.startswith('linux') and 'PYOPENGL_PLATFORM' not in os.environ:
    qtPlatform = os.environ.get('QT_QPA_PLATFORM', '').split(';')[0].split(':')[0].strip().lower()
    if not qtPlatform.startswith('wayland'):
        os.environ['PYOPENGL_PLATFORM'] = 'glx'
