# KinegamiPython
[![Powered by the Spatial Math Toolbox](https://github.com/bdaiinstitute/spatialmath-python/raw/master/.github/svg/sm_powered.min.svg)](https://github.com/bdaiinstitute/spatialmath-python)

Daniel A. Feshbach, Wei-Hsi Chen, Daniel E. Koditschek, Cynthia R. Sung. “Kinegami: Open-source Software for Creating Kinematic Chains from Tubular Origami.” Submitted to Origami 8: Eigth International Meeting on Origami in Science, Mathematics and Education, 2024 (under review).

This is a python repository for creating and modifying kinematic chains made of tubular origami. Examples of its usage can be found in the examples folder. 
It requires the following libraries:
- numpy
- scipy
- matplotlib
- roboticstoolbox-python
- ezdxf (installed as pip3 install ezdxf, see https://ezdxf.mozman.at/docs/setup.html)
We developed and tested this code in python 3.12.

Abstract: Arms, legs, and fingers of animals and robots are all examples of “kinematic chains” - mechanisms with sequences of joints connected by rigid links. Lightweight kinematic chains can be manufactured quickly and cheaply using origami patterns rolled up into tubes. In recent work [[Chen et al. 22](https://repository.upenn.edu/entities/publication/66231fb3-9120-40e0-87c8-0d4f5c058983)], we demonstrated that origami patterns for kinematic chains with arbitrary numbers of degrees of freedom can be constructed algorithmically from a minimal kinematic specification (the axes that joints rotate about or translate along). The work was founded on a catalog of tubular crease patterns for revolute joints (rotation about an axis), prismatic joints (translation along an axis), and links, which compose together to form kinematic chains. In this paper, we implement these patterns and algorithms in python as a user-friendly, open-source tool for creating kinematic chains out of origami with little origami design expertise. Users can specify each joint by its axis only or by its full location and orientation. Our software uses this information to construct a single crease pattern for the corresponding chain. The software also includes a visualization tool so users can check that the chain can achieve their desired configurations, and methods to try modifications of an existing chain. This paper provides a detailed explanation of the code and its usage, including an explanation of our proposed representation for tubular crease patterns. We include a number of examples to illustrate the software’s capabilities and potential for robot and mechanism design.

## Code Structure Overview
The central class for users to interact with the code is `KinematicChain`, which inherits from `KinematicTree` (which is where most of the methods are defined). A chain stores a list of `Joint` objects, which each store their pose. When adding new joints, `KinematicChain` also computes links in the shape of [CSC Dubins Paths](https://en.wikipedia.org/wiki/Dubins_path), and stores them as a list of `LinkCSC` objects. `KinematicChain` has a `show` method for visualization. To make a chain, initialize it with a starting joint and then append further joints using the `append` method on `KinematicChain`. This method has options to place the new joint at its stored pose or to use a joint placement algorithm to find a pose on its axis of motion. You can also edit an existing `KinematicChain` by moving or deleting joints using the methods `translateJointAlongAxisOfMotion`, `rotateJointAboutAxisOfMotion`, `transformJoint`, and `delete`, and can change joint states using the `setJointState` method to visualize different configurations of the same chain. Once you are satisfied with a chain design, you can use the `creasePattern` method to compute its crease pattern as a `TubularPattern` object, which has its own methods to `show` a pattern or `save` it as a DXF.

## System Overview and Definitions
This code enables users to create kinematic chains made of tubular origami. 
The chains have revolute and prismatic joints connected by links following CSC Dubins paths. 
To iteratively construct a chain, users append each joint either by specifying its exact pose or by specifying only its axis of motion and letting an algorithm place the joint along the axis. Users can input poses and axes either in global coordinates or relative to the previous joint's pose.
Additionally, users can modify existing chains by moving or deleting joints. Joint movements can be given as translations along the axis of motion, rotations about the axis of motion, or arbitrary rigid transformations. Joint movements can be set either to propagate to all subsequent joints or to apply only to the given joint. 

If a proposed modification cannot generate feasible link paths, the system will reject the change and warn the user. 
When a user has a candidate design, they can specify the state of each joint to visualize different configurations of the chain, letting them check whether it can do what they want and modify the design accordingly. Here are some example chains visualized in our system (the code for these is in the `examples/` folder): ![Gallery of example kinematic chains](figures/gallery.png)
