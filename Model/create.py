'''
Create an HDF5 file to store model performance in.
'''

#   SETTINGS

#   Models file.
filename   = '../Data/models.h5'
group_name = 'models'

import h5py

fp    = h5py.File(filename, 'w')
group = fp.create_group(group_name)

group.attrs['info'] = '''Contained are chronological groups containing training and test data 
                         performance criteria. The name of the subgroup is the name of the model. 
                         It has a timestamp attribute and min max attributes required to undo 
                         normalization.'''

fp.close()
