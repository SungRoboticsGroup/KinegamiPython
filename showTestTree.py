from KinematicTree import *
from optimizationFunctions import *
import json

title = "5 Joint Generalized Gimbal Chains Cube Size 10"
path = "sim_results/" + title + "/0.tree"
construct = loadKinematicTree(path)
construct.show(block=False, showScaleBar=False, jointAxisScale=40)

path = "sim_results/" + title + "/0/Quadratic/final.tree"
sg = loadKinematicTree(path)
sg.show(block=True, showScaleBar=False)


"""
optimizations = [partial(squaredOptimize, childParentRatio=0,streamline=True,guarantee=True),
                partial(squaredOptimize, childParentRatio=0,streamline=True,guarantee=False),
                partial(squaredOptimize, childParentRatio=0,streamline=False,guarantee=False),
                partial(squaredOptimize, childParentRatio=0,streamline=False,guarantee=False,resetOnFail=False),
                partial(squaredOptimize, childParentRatio=0,streamline=False,guarantee=True),
                partial(linearOptimize, childParentRatio=0, streamline=False,guarantee=False)]
labels = ["Streamline + Guarantee (SG)", 
          "Streamline No Guarantee (SNG)",
          "No Streamline No Guarantee (NSNG)",
          "NSNG, No Reset on Fail",
          "No Streamline Guarantee (NSG)",
          "Linear (L)"]


results = [[],[],[],[],[],[],[]]
i=6
for no, f in enumerate(optimizations):
    print(f"\nTrying loss function {no}")
    direc = "sim_results/" + title + "/" + str(i) + "/" + labels[no] + "/"
    os.makedirs(direc, exist_ok=True)
    optimized, times, losses = f(construct, showSteps=False, parallelize=True, evaluate=True, verbose=False, directory=direc)
    results[i].append((times, losses))

with open("sim_results/" + title + "/random_results_chkpt" + str(i) + ".json", "w") as file:
    json.dump(results, file)

optimized, times, losses = f(construct, showSteps=False, parallelize=True, evaluate=True, verbose=False)
print(optimized)
print(times)
print(losses)"
"""