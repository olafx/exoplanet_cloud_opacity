'''
Grid surface plot of grid molecule pair datasets.
Not so important, since grid datasets are not so important.
'''

#   SETTINGS

filename         = '../Data/grid.h5'
main_group_name  = 'main'
data_name        = 2
#   Number of entries in the grid. (Resolution.)
entries          = 40
#   Number of wavelengths in grid.
n_wavelengths    = 20
#   Controls which molecule to look at (ratio). 0 to only look at first molecule, entries - 1 to only
#   look at second molecule.
offset           = 0
#   Plot scattering? True. Plot absorption? False.
plot_scattering  = False
#   Plot gradient instead? True. No? False.
plot_gradient    = False
#   Show figures? True. Save figures? False.
show             = False
#   Azimuthal and polar angle of view.
azimuthal, polar = -140, 20

import h5py
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fp         = h5py.File(filename, 'r')
main_group = fp[main_group_name]
group      = main_group[str(data_name)]
data_in    = group['in']
data_out   = group['out']
n          = data_in.shape[0]

#   every n_wavelengths * entries,   return to the same mixture with same size but different size sd
#   every n_wavelengths * entries^2, return to the same mixture with different size but same size sd

sizes        = np.empty([entries])
wavelengths  = np.empty([n_wavelengths])
compositions = np.empty([entries])
values       = np.empty([entries, n_wavelengths])
fraction     = n_wavelengths * entries

for i in range(0, n_wavelengths):
    wavelengths[i] = data_in[i, 0]

j = 0
for i in range(0, n, fraction * entries):
    i += n_wavelengths * offset
    sizes[j]         = data_in[i, 1]
    compositions[j]  = data_in[i, 3]
    values[j, :]     = data_out[i : i + n_wavelengths, 0 if plot_scattering else 1]
    j += 1

fp.close()

if plot_gradient:
    values = np.gradient(values)[0]

ax = (fig := plt.figure(figsize=(6, 5))).add_subplot(111, projection='3d')
mesh_wavelengths, mesh_sizes = np.meshgrid(wavelengths, sizes)
surf = ax.plot_surface(mesh_wavelengths, mesh_sizes, values, cmap='jet', alpha=0.9, linewidth=0, antialiased=True)
ax.set_xlabel('wavelength (μm)')
ax.set_ylabel('dust size (μm)')
ax.set_zlabel('relative scattering coefficient' if plot_scattering else 'relative absorption coefficient')
ax.azim, ax.elev = azimuthal, polar
plt.show() if show else fig.savefig('fig.pdf', bbox_inches='tight')
