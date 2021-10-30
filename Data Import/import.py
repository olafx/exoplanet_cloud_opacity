'''
Add numerical data to existing HDF5 data file.
'''

import pytz

#   Export data file.
filename          = '../Data/data.h5'
main_group_name   = 'main'
#   Import data settings file.
settings_filename = '../Data Generation 2/data.txt'
#   Import data input and output.
input_filename    = '../Data Generation 2/input_set.txt'
output_filename   = '../Data Generation 2/training_set.txt'
#   Computer timezone. (Needed to keep consistent timestamps.)
source_zone       = pytz.timezone('Europe/Brussels')

import h5py
from datetime import datetime
import numpy as np

target_zone = pytz.timezone('UTC')

fp = h5py.File(filename, 'a')

#   Verifying export data file.
if main_group_name not in fp.keys():
    raise Exception(f"group '{main_group_name}' does not exist in file, or file does not exist")
main_group = fp[main_group_name]
n_groups = len(main_group)

#   Verify consistency of export data file; previous dataset should exist and the new one shouldn't
#   override anything.
if 0 < n_groups:
    if str(n_groups) not in main_group:
        raise Exception(f"'{n_groups}' not in '{main_group_name}'")
    if str(n_groups + 1) in main_group:
        raise Exception(f"'{n_groups + 1}' already in '{main_group_name}'")

#   Read.
#   The Fortran floats are read as strings, then modified so that numpy can read them, then read as
#   actual floats.
data_in           = np.genfromtxt(input_filename,                                dtype=np.float64)
data_out          = np.genfromtxt(output_filename,                               dtype=np.float64)
rng_param_a       = np.genfromtxt(settings_filename, skip_header=1,  max_rows=1, dtype=np.str).item()
rng_param_a       = np.asfarray(  rng_param_a[:-2],                              dtype=np.float64)
rng_param_b       = np.genfromtxt(settings_filename, skip_header=3,  max_rows=1, dtype=np.str).item()
rng_param_b       = np.asfarray(  rng_param_b[:-2],                              dtype=np.float64)
data_gen_func     = np.genfromtxt(settings_filename, skip_header=5,  max_rows=1, dtype=np.str).item()
n_training        = np.genfromtxt(settings_filename, skip_header=7,  max_rows=1, dtype=np.int64).item()
dust_size_min     = np.genfromtxt(settings_filename, skip_header=9,  max_rows=1, dtype=np.str).item()
dust_size_min     = np.asfarray(  dust_size_min[:-2],                            dtype=np.float64).item()
dust_size_max     = np.genfromtxt(settings_filename, skip_header=11, max_rows=1, dtype=np.str).item()
dust_size_max     = np.asfarray(  dust_size_max[:-2],                            dtype=np.float64).item()
dust_size_sd_max  = np.genfromtxt(settings_filename, skip_header=13, max_rows=1, dtype=np.str).item()
dust_size_sd_max  = np.asfarray(  dust_size_sd_max[:-2],                         dtype=np.float64).item()
wavelength_min    = np.genfromtxt(settings_filename, skip_header=15, max_rows=1, dtype=np.str).item()
wavelength_min    = np.asfarray(  wavelength_min[:-2],                           dtype=np.float64).item()
wavelength_max    = np.genfromtxt(settings_filename, skip_header=17, max_rows=1, dtype=np.str).item()
wavelength_max    = np.asfarray(  wavelength_max[:-2],                           dtype=np.float64).item()
mixing_ratios_max = np.genfromtxt(settings_filename, skip_header=20, max_rows=1, dtype=np.str)
mixing_ratios_max = np.asfarray(  mixing_ratios_max,                             dtype=np.float64)

if data_in.shape[0] != data_out.shape[0]:
    raise Exception(f'input and output data have different dimensions ({data_in.shape[0]} vs {data_out.shape[0]})')
n_data = data_in.shape[0]

#   Create new group.
group = main_group.create_group(str(n_groups + 1))

#   Write group attributes.
group.attrs       ['creation date'] =                      str(source_zone.localize(datetime.now()).astimezone(target_zone))
group.attrs       ['data generation function'] =           data_gen_func
group.attrs.create('training set size',                    n_training,        dtype=np.uint32)
group.attrs.create('dust size minimum',                    dust_size_min,     dtype=np.float64)
group.attrs.create('dust size maximum',                    dust_size_max,     dtype=np.float64)
group.attrs.create('dust size standard deviation maximum', dust_size_sd_max,  dtype=np.float64)
group.attrs.create('wavelength minimum',                   wavelength_min,    dtype=np.float64)
group.attrs.create('wavelength maximum',                   wavelength_max,    dtype=np.float64)
group.attrs.create('mixing ratios maximum',                mixing_ratios_max, dtype=np.float64)
group.attrs.create('RNG param a',                          rng_param_a,       dtype=np.float64)
group.attrs.create('RNG param b',                          rng_param_b,       dtype=np.float64)

#   Write group data.
#   Sometimes this doesn't work because szip is proprietary. Make sure HDF5 was built with NASA szip.
group.create_dataset('in',  data=data_in,  dtype=np.float64, compression='szip');
group.create_dataset('out', data=data_out, dtype=np.float64, compression='szip');

fp.close()
