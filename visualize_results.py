import json
import os
from matplotlib import pyplot as plt

title = "10 Joint Trees Dense"

folder = "sim_results/" + title + "/images/"

# Ensure the folder exists
if not os.path.exists(folder):
    os.makedirs(folder)

with open("sim_results/" + title + "/random_results.json", 'r') as file:
    results = json.load(file)

labels = ["Streamline + Guarantee (SG)", 
          "Streamline No Guarantee (SNG)",
          "No Streamline No Guarantee (NSNG)",
          "NSNG, No Reset on Fail",
          "No Streamline Guarantee (NSG)",
          "Linear (L)"]
for i, result in enumerate(results):
    plt.figure(figsize=(8, 5)) 
    
    idx = 0
    for x, y in result:
        plt.plot(x, y, marker='o', linestyle='-', label=labels[idx])
        idx += 1

    plt.xlabel('Time')
    plt.ylabel('Loss')
    plt.title('Loss vs Time: ' + str(title))
    plt.legend()

    plt.grid(True)
    plt.savefig(folder + str(i) + ".png")