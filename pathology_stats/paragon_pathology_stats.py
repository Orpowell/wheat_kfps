import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from statsmodels.stats.multicomp import pairwise_tukeyhsd

#DATA

Paragon	= [2,1,1,2,1,1]
Kronos= [5, 5, 6, 6, 6]
J3_48 = [6, 5, 6, 6, 6]
J5_48 = [4, 3, 3, 5, 6, 4]


data = [Paragon, Kronos, J3_48, J5_48]
labels = ['Paragon', 'Kronos', 'J3-48', 'J5-48']

groups = dict(zip(labels, data))

groups = {k : np.asarray(v) for k,v in groups.items()}

names = list(groups.keys())
means = [g.mean() for g in groups.values()]
sems  = [stats.sem(g) for g in groups.values()]
labels = sum([[k] * len(v) for k,v in groups.items()], [])

# Tukey HSD
all_data = np.concatenate(list(groups.values()))

tukey = pairwise_tukeyhsd(all_data, labels, alpha=0.05)

fig, ax = plt.subplots(figsize=(7, 6))
x_pos = np.arange(len(names))


for i, (name, data) in enumerate(groups.items()):
    #jitter = np.random.uniform(-0.12, 0.12, size=len(data))
    #ax.scatter(np.full_like(data, i) + jitter, data, s=18,
    #           color='black', alpha=0.45, zorder=5)
    sns.stripplot(x=i, y=data, color='#000000', ax=ax, s=3)

ax.bar(x_pos, means, yerr=sems, width=0.55, capsize=5, edgecolor='black', lw=0.8, color=['#bbbbbb' for i in range(0, len(groups))])

print(tukey)

# Significance brackets for significant pairs
y_max = max(m + s for m, s in zip(means, sems)) + 0.5
step = 0.45
sig_pairs = [(r[0], r[1], r[3], r[6]) for r in tukey._results_table.data[1:] if r[1] == 'Paragon']  # reject = True

for idx, (g1, g2, pval, bool_diff) in enumerate(sig_pairs):
    x = names.index(g1) - 0.005
    y = max(groups[g1]) + 0.25
    #y = y_max

    print(idx,(g1, g2, pval, bool_diff) )

    if bool_diff:
        print('Ahoy!')
        plt.text(x, y, 'b')

    else:
        plt.text(x, y, 'a')

plt.text(0, max(groups['Paragon']) + 0.25, 'a')

ax.set_xticks(x_pos)
ax.set_xticklabels(names)
ax.set_ylim(0,6.3)
plt.xticks(rotation=45, ha='right')
ax.set_ylabel('Disease Scores')
ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout()
plt.savefig('Paragon_deletion_mutant_leaf_scores.svg', dpi=300)