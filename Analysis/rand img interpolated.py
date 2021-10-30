'''
Grid image plot of random data through interpolation.
Intended for datasets of a single molecule.
'''

#   SETTINGS

filename        = '../Data/main.h5'
main_group_name = 'main'
data_name       = 2 + 1
#   Plot scattering? True. Plot absorption? False.
plot_scattering = False
#   Grid resolution for plotting.
grid_n          = 40
#   Show figures? True. Save figures? False.
show            = True

import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

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
values      = data_out_filtered[:, 0 if plot_scattering else 1]

#   Interpolation.
grid_wavelengths, grid_sizes = np.meshgrid(np.linspace(np.min(wavelengths), np.max(wavelengths), grid_n), np.linspace(np.min(sizes), np.max(sizes), grid_n))
grid_interpolated = griddata(points=(wavelengths, sizes), values=values, xi=(grid_wavelengths, grid_sizes), method='linear', fill_value=0)

ax = (fig := plt.figure(figsize=(6, 5))).add_subplot(111)
img = plt.pcolormesh(grid_wavelengths, grid_sizes, grid_interpolated, shading='auto')
cb = fig.colorbar(img)
cb.set_label('relative scattering coefficient' if plot_scattering else 'relative absorption coefficient')
ax.set_xlabel('wavelength (micron)'), ax.set_ylabel('size (micron)')
plt.show() if show else fig.savefig('fig.pdf', bbox_inches='tight')
