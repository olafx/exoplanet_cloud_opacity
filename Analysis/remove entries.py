'''
Remove entries in a dataset by replacing them with a copy of the previous.
Useful to fix results where calculation failed, similar to fix NaN or inf.
'''

#   SETTINGS

filename        = '../Data/data.h5'
main_group_name = 'main'
group           = 9999
#   Indexes of points to remove.
bad_points = []

import h5py
import numpy as np
import math

fp         = h5py.File(filename, 'r+')
main_group = fp[main_group_name]
group      = main_group[str(group)]
data_in    = group['in']
data_out   = group['out']

print(bad_points)

for bad_point in bad_points:
    data_in[bad_point]  = data_in[bad_point - 1]
    data_out[bad_point] = data_out[bad_point - 1]

fp.close()
