'''
Replaces entries with NaN or inf absorption or scattering coefficients with copies of their
previous entry.
'''

#   SETTINGS

filename        = '../Data/data.h5'
main_group_name = 'main'
group           = 2

import h5py
import numpy as np
import math

fp         = h5py.File(filename, 'r+')
main_group = fp[main_group_name]
group      = main_group[str(group)]
data_in    = group['in']
data_out   = group['out']
bad_points = []

print('max scattering', np.max(data_out[:, 0]))
print('max absorption', np.max(data_out[:, 1]))

for j in range(2):
    for i in range(len(data_out)):
        if math.isnan(data_out[i, j]) or math.isinf(data_out[i, j]):
            print('bad', i)
            bad_points.append([i, j])

#   Untested, but potentially faster.
#       indexes = np.arange(len(data_out))
#       indexes = indexes[math.isnan(data_out[indexes]) or math.isinf(data_out[indexes])]

print('bad points\n', bad_points)

for point in bad_points:
    i = point[0]
    data_in[i] = data_in[i-1]
    data_out[i] = data_out[i-1]

fp.close()
