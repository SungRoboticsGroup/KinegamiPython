# KinegamiPython
[![Powered by the Spatial Math Toolbox](https://github.com/bdaiinstitute/spatialmath-python/raw/master/.github/svg/sm_powered.min.svg)](https://github.com/bdaiinstitute/spatialmath-python)

This branch builds on the code from:

Daniel A. Feshbach, Wei-Hsi Chen, Daniel E. Koditschek, Cynthia R. Sung. “Kinegami: Open-source Software for Creating Kinematic Chains from Tubular Origami.” In Origami 8: Eighth International Meeting on Origami in Science, Mathematics and Education, 2024.

This is a Python+OpenSCAD repository for creating and modifying kinematic chains and trees made of tubular origami and/or 3D printing. Examples of its usage can be found in the examples folder.

Python requirements (which you can install all at once with `pip install -r requirements.txt`):
- numpy
- scipy
- matplotlib
- spatialmath-python
- ezdxf
- pyswarms

We developed and tested this code in python 3.12.

OpenSCAD installers: https://openscad.org/downloads.html

Once you install OpenSCAD, edit the file openscad.bat to set it to the path to your openscad executable.

Note (current as of 11/26/2024): The development release of OpenSCAD has a new option to use the manifold library as its backend, which speeds up 3D printing file generation by 2 orders of magnitude. To use this option with our code, make sure you have the Nightly build for OpenSCAD (under Development Snapshots in https://openscad.org/downloads.html). Then set the optional paramter `manifold=True` when calling the method `export3DKinematicTree` on a `KinematicTree[PrintedJoint]` object.

If you're on WSL, you may need to do use apt to install qtbase5-dev (if apt says it cannot locate that package, try running sudo apt update and then retrying).