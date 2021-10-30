'''
Probability density distribution for MgSiO3 or max of any molecule.
Probability to get certain minimum mixture ratio for one molecule or any molecule.
Used to analyze how realistic the mixture ratios being generated are.
Ideally get lots of data where at least one molecule has a relatively large mixture ratio.
'''

import h5py
import matplotlib.pyplot as plt
import numpy as np

#   Data file.
filename        = '../Data/data.h5'
main_group_name = 'main'
#   Dataset to use.
data_name       = 2
#   Number of bins for histogram.
n_bins          = 200
#   Show figures? True. Save figures? False.
show            = False

fp                 = h5py.File(filename, 'r')
main_group         = fp[main_group_name]
group              = main_group[str(data_name)]
mixture_ratios     = group['in'][:, 3]
mixture_ratios_max = np.max(group['in'][:, 3 :], axis=1)
rng_param_a        = group.attrs['RNG param a']
rng_param_b        = group.attrs['RNG param b']
fp.close()

n = len(mixture_ratios)

print('% with at least 20% mixture ratio for MgSiO3:',       np.round(100 * mixture_ratios[mixture_ratios > 0.2].shape[0] / n, 1))
print('% with at least 20% mixture ratio for one molecule:', np.round(100 * mixture_ratios_max[mixture_ratios_max > 0.2].shape[0] / n, 1))

ax = (fig1 := plt.figure(figsize=(4, 3))).add_subplot(111)
ax.hist(mixture_ratios, bins=n_bins, range=(0, 1), density=True)
ax.set_xlim([0, 1])
ax.set_xlabel('MgSiO3 mixture ratio')
ax.set_ylabel('occurence (%)')
plt.title(f'MgSiO3 mixture ratio (tan $\pi/2$ {rng_param_a} to {rng_param_b} RNG)')
if show:
    plt.show()

ax = (fig2 := plt.figure(figsize=(4, 3))).add_subplot(111)
ax.hist(mixture_ratios_max, bins=n_bins, range=(0, 1), density=True)
ax.set_xlim([0, 1])
ax.set_xlabel('max mixture ratio of any molecule')
ax.set_ylabel('occurence (%)')
plt.title(f'max mixture ratio of any molecule (tan $\pi/2$ {rng_param_a} to {rng_param_b} RNG)')
plt.show() if show else fig1.savefig('fig1.pdf', bbox_inches='tight'), fig2.savefig('fig2.pdf', bbox_inches='tight')
