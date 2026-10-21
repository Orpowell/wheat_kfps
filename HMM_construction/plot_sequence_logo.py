import logomaker
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({'font.size': 14})
data = "/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/analysis/beta_finger_HMM_refinement/skylign_matrix.txt"

df = pd.read_csv(data, skiprows=4, sep="\t", header=None)

df.columns = ["pos", 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'Expected Insert Length', 'Insert Probability', 'Delete Probability', 'Model Mask']

info: pd.Series = df[['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']]

info.index = info.index + 1

logo = logomaker.Logo(info, font_name='Helevetica', color_scheme='skylign_protein')

logo.ax.set_ylabel('Information (bits)')
logo.ax.set_xlabel('Sequence Position')
logo.ax.set_yticks([1,2,3,4,5])
logo.ax.set_ylim(0, 5.01)
logo.ax.set_xlim(0, 31.5)
logo.style_spines(visible=False)
logo.style_spines(spines=['left', 'bottom'], visible=True)
logo.ax.set_xticks(np.arange(1, 32, 1))

plt.tight_layout()
plt.savefig("/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/analysis/beta_finger_HMM_refinement/refined_beta_finger.png", dpi=300)
plt.savefig("/Users/powellor/Desktop/kfp_annotation/KFP_annotation_v2/beta_finger_HMM/analysis/beta_finger_HMM_refinement/refined_beta_finger.svg")
plt.show()
