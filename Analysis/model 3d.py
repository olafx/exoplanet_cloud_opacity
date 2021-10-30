'''
Grid surface plot of model with controllable range.
'''

#   SETTINGS

#   Input range to plot.
range_wavelength       = [0.01, 1.0]
range_size             = [0.001, 1.0]
#   Pick one dust size standard deviation.
dust_size_sd           = 0.01
mixture_ratios         = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
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
grid_n                 = 40
#   Show figures? True. Save figures? False.
show                   = True
#   Azimuthal and polar angle of view.
azimuthal, polar       = -140, 20

import h5py
import tensorflow as tf
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

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

#   Wavelengths and sizes are exponentially distributed, so linearly for the neural network.
wavelengths = np.exp(np.linspace(np.log(range_wavelength[0]), np.log(range_wavelength[1]), grid_n))
sizes       = np.exp(np.linspace(np.log(range_size[0]),       np.log(range_size[1]),       grid_n))

#   Minmax and logarithmic scaling.
wavelengths_model  = (np.log(wavelengths)               - np.log(model_range_wavelength[0])) /\
                     (np.log(model_range_wavelength[1]) - np.log(model_range_wavelength[0]))
sizes_model        = (np.log(sizes)                     - np.log(model_range_size[0])) /\
                     (np.log(model_range_size[1])       - np.log(model_range_size[0]))
size_sd_model      = dust_size_sd / model_max_size_sd

#   Creating data for the grid. Will be resized later; needs to go through model first.
in_physical, in_chemical = np.empty([grid_n**2, 3]), np.empty([grid_n**2, 15])
for i in range(len(in_physical)):
    in_chemical[i, :] = mixture_ratios
for i in range(0, len(in_physical), grid_n):
    in_physical[i : i + grid_n, 0] = wavelengths_model
for i in range(grid_n):
    in_physical[i * grid_n : (i + 1) * grid_n, 1] = sizes_model[i]
in_physical[:, 2] = size_sd_model

#   Compute values using model and resize for the grid.
values = model([in_physical, in_chemical])
values = np.resize(values[:, 0 if plot_scattering else 1], (grid_n, grid_n))

ax = (fig := plt.figure(figsize=(6, 5))).add_subplot(111, projection='3d')
mesh_wavelengths, mesh_sizes = np.meshgrid(wavelengths, sizes)
surf = ax.plot_surface(mesh_wavelengths, mesh_sizes, values, cmap='jet', alpha=0.9, rcount=grid_n, ccount=grid_n, linewidth=0, antialiased=True)
ax.set_xlabel('wavelength (μm)')
ax.set_ylabel('size (μm)')
ax.set_zlabel('relative scattering coefficient' if plot_scattering else 'relative absorption coefficient')
ax.azim, ax.elev = azimuthal, polar
plt.show() if show else fig.savefig('fig.pdf', bbox_inches='tight')
