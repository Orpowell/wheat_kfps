import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from statsmodels.stats.multicomp import pairwise_tukeyhsd

def p_to_stars(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

def add_bracket(ax, x1, x2, y, h, text):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.2, color='black')
    ax.text((x1 + x2) / 2, y + h, text, ha='center', va='bottom', fontsize=10)

#DATA
Cadenza = [1,3,3,3,2,3,1,1,1,1]
Kronos = [5,5,5,5,6,6,6,5,6,6]
Cadenza1637 = [6,5,4,5,5,4,2]
Cadenza0723 = [4,3,4,6,6,6,5,5,4,5]
Cadenza1617 = [5,3,2,1,1,4,2,4,1,1]
Cadenza0911 = [1,1,4,2,5,6]
Cadenza2030 = [5,1,6,5,6,4,4,4,2,3]
Cadenza0346 = [6,5,5,6,6,6,6,5,6,5]
Cadenza1269 = [5,2,2,1,5,4,5,5,5,4]
Cadenza0896 = [5,5,5,4,3,6,4,6,6,4]


data = [Cadenza,Kronos,Cadenza1637,Cadenza0723,Cadenza1617,Cadenza0911,Cadenza2030,Cadenza0346,Cadenza1269,Cadenza0896]
labels = ['Cadenza','Kronos','Cadenza1637','Cadenza0723','Cadenza1617','Cadenza0911','Cadenza2030','Cadenza0346','Cadenza1269','Cadenza0896']

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
sig_pairs = [(r[0], r[1], r[3], r[6]) for r in tukey._results_table.data[1:] if r[0] == 'Cadenza']  # reject = True

for idx, (g1, g2, pval, bool_diff) in enumerate(sig_pairs):
    x = names.index(g2) - 0.1
    y = max(groups[g2]) + 0.25
    #y = y_max

    if bool_diff:
        plt.text(x, y, 'b')

    else:
        plt.text(x, y, 'a')

plt.text(0-0.1, max(groups['Cadenza']) + 0.25, 'a')

   

ax.set_xticks(x_pos)
ax.set_xticklabels(names)
ax.set_ylim(0,6.05)
plt.xticks(rotation=45, ha='right')
ax.set_ylabel('Disease Scores')
ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout()
plt.savefig('Cadenza_deletion_mutant_leaf_scores.svg', dpi=300)