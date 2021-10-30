'''
Residual grid image plot of residual between model and random dataset.
Intended for datasets of a single molecule.
'''

#   SETTINGS

#   Data file.
filename               = '../Data/main.h5'
main_group_name        = 'main'
#   Dataset to use.
data_name              = 2 + 16 + 2
#   File for models (min max properties are stored here).
filename_models        = '../Data/models.h5'
models_main_group_name = 'models'
#   Folder for models.
foldername_models      = '../Data/Models'
#   Model file.
filename_model         = '2.h5'
#   Plot scattering? True. Plot absorption? False.
plot_scattering        = False
#   Grid resolution for plotting.
grid_n                 = 100
#   Show figures? True. Save figures? False.
show                   = True

import h5py
import tensorflow as tf
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.tri as tri

#   Getting min max properties.
fp         = h5py.File(filename_models, 'r')
main_group = fp[models_main_group_name]
group      = main_group[filename_model[:-3]]
model_range_wavelength = group.attrs['range_wavelength']
model_range_size       = group.attrs['range_size']
model_max_size_sd      = group.attrs['max_size_sd']
fp.close()

#   Loading model.
model = tf.keras.models.load_model(foldername_models + '/' + filename_model)

#   Setting up to get data.
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

#   Getting data.
wavelengths = data_in_filtered[:, 0]
sizes       = data_in_filtered[:, 1]
sizes_sd    = data_in_filtered[:, 2]
ratios      = data_in_filtered[:, 3 :]
scattering  = data_out_filtered[:, 0]
absorption  = data_out_filtered[:, 1]

#   Preprocessing input data.
wavelengths_model = (np.log(wavelengths)               - np.log(model_range_wavelength[0])) /\
                    (np.log(model_range_wavelength[1]) - np.log(model_range_wavelength[0]))
sizes_model       = (np.log(sizes)                     - np.log(model_range_size[0])) /\
                    (np.log(model_range_size[1])       - np.log(model_range_size[0]))
sizes_sd_model    = sizes_sd / model_max_size_sd
ratios_model      = ratios[:,:-1]

#   Getting output from model.
input_physical       = np.empty([wavelengths.shape[0], 3])
input_physical[:, 0] = wavelengths_model
input_physical[:, 1] = sizes_model
input_physical[:, 2] = sizes_sd_model
model_output         = model([input_physical, ratios_model])
values               = model_output[:, 0 if plot_scattering else 1]

#   Residual.
values -= scattering if plot_scattering else absorption

#   Make linear grid and interpolate.
grid_wavelengths, grid_sizes = np.meshgrid(np.linspace(np.min(wavelengths), np.max(wavelengths), grid_n), np.linspace(np.min(sizes), np.max(sizes), grid_n))
grid_interpolated = griddata(points=(wavelengths, sizes), values=values, xi=(grid_wavelengths, grid_sizes), method='linear', fill_value=0)

ax = (fig := plt.figure(figsize=(6, 5))).add_subplot(111)
img = plt.pcolormesh(grid_wavelengths, grid_sizes, grid_interpolated)
cb = fig.colorbar(img)
cb.set_label('relative scattering coefficient' if plot_scattering else 'relative absorption coefficient')
ax.set_xlabel('wavelength (μm)')
ax.set_ylabel('size (μm)')
plt.show() if show else fig.savefig('fig.pdf', bbox_inches='tigh')
