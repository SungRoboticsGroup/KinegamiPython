# Kinegami8OSME
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
