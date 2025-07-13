import matplotlib.pyplot as plt
import numpy as np

stats = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/analysis/beta_finger_HMM_refinement/beta_finger_refinement.stats"

with open(stats) as input:
    y = list()
    x = list()

    for line in input:
        data = line.split()
        x.append(data[0])
        y.append(int(data[1]))

plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False


plt.scatter(x,y, marker="x", linewidths=1, color="black")
plt.plot(x,y, linewidth=0.5, color='black')
plt.xlabel("Iteration")
plt.ylabel("# of Sequences")
plt.savefig("/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/analysis/beta_finger_HMM_refinement/beta_finger_iterative_refinement.png")
plt.show()