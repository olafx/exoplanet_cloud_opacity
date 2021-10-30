'''
Unstructured surface plot of random data.
Intended for datasets of a single molecule.
'''

#   SETTINGS

filename         = '../Data/main.h5'
main_group_name  = 'main'
data_name        = 2 + 1
#   Plot scattering? True. Plot absorption? False.
plot_scattering  = False
#   Show figures? True. Save figures? False.
show             = False
#   Azimuthal and polar angle of view.
azimuthal, polar = -140, 20

import h5py
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.tri as tri

#   Reading data.
fp         = h5py.File(filename, 'r')
main_group = fp[main_group_name]
group      = main_group[str(data_name)]
data_in    = group['in']
data_out   = group['out']

#   Filter out data with sufficiently small dust size standard deviation (3% of max).
index_filtered    = np.arange(data_in.shape[0])[data_in[:, 2] < 0.03 * group.attrs['dust size standard deviation maximum']]
data_in_filtered  = data_in[index_filtered]
data_out_filtered = data_out[index_filtered]

print('filtered', data_in_filtered.shape[0])

fp.close()

wavelengths = data_in_filtered[:, 0]
sizes       = data_in_filtered[:, 1]
scattering  = data_out_filtered[:, 0]
absorption  = data_out_filtered[:, 1]

ax = (fig := plt.figure(figsize=(6, 5))).add_subplot(111, projection='3d')
grid = tri.Triangulation(wavelengths, sizes)
ax.plot_trisurf(grid, absorption, cmap='jet', alpha=0.9, linewidth=0, antialiased=True)
ax.set_xlabel('wavelength (μm)')
ax.set_ylabel('size (μm)')
ax.set_zlabel('relative scattering coefficient' if plot_scattering else 'relative absorption coefficient')
ax.azim, ax.elev = azimuthal, polar
plt.show() if show else fig.savefig('fig.pdf', bbox_inches='tight')
