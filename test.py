import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)

# Example 1A
from KinematicTree import *
from testqtgraph import *
from makeKinematicTree import *

tree = loadKinematicTree("unchangedHexapod")

tree.show()

tree = loadKinematicTree("optimizedHexapod9")

tree.show()
# #tree.show()
# #tree.show(showSpecificCapsules=([(23, 0)],[(27, 5)]))

# for i in range(0, len(tree.Joints)):
#     specificJointIndices = [i]
#     capsules = tree.selectCollisionCapsules(specificJointIndices=specificJointIndices, ignoreLater=False)
#     #print(f"checking {i}")
#     if tree.detectCollisionsWithCapsules(specificJointIndices, capsules, show=True, debug=True) > 0:
#         print(i)

# # print(separatingAxisTheorem(tree.Joints[23].collisionCapsules[0].box, tree.Links[27].collisionCapsules[5].box))
# # print(tree.Joints[23].collisionCapsules[0].box.points, tree.Links[27].collisionCapsules[5].box.points)
# # if tree.Joints[23].collisionCapsules[0].collidesWith(tree.Links[27].collisionCapsules[5])[0]:
# #     print("SUCCESS")


# # specificJointIndices = [7]
# # capsules = tree.selectCollisionCapsules(specificJointIndices=specificJointIndices, ignoreLater=True)
# # print(tree.detectCollisionsWithCapsules(specificJointIndices, capsules, show=True))
# #tree.show()
# # print(tree.Joints[7].collisionCapsules[0].collidesWith(tree.Joints[2].collisionCapsules[0]))
# # print(tree.Joints[2].collisionCapsules[0].collidesWith(tree.Joints[7].collisionCapsules[0]))


# # for i in range(0, len(tree.Links[7].collisionCapsules)):
# #     for j in range(0, len(tree.Joints[2].collisionCapsules)):
# #         if tree.Links[7].collisionCapsules[i].collidesWith(tree.Joints[2].collisionCapsules[j])[0]:
# #             print(i, j)
# #         if tree.Joints[2].collisionCapsules[j].collidesWith(tree.Links[7].collisionCapsules[i])[0]:
# #             print(i, j)

# #tree.show(addCapsules=[tree.Joints[4].collisionCapsules[0], tree.Joints[3].collisionCapsules[0]])