'''
Create new file to store data in.
'''

#   Data file.
filename = '../Data/data.h5'
group_name = 'main'

import h5py

fp = h5py.File(filename, 'w')
group = fp.create_group(group_name)

group.attrs['info']           = '''Contained are chronological groups containing input and output data. 
                                   Input data creation settings and timestamp are stored as attributes to the group. 
                                   Rows in input correspond to rows in output.'''
group.attrs['columns input']  = ['wavelength (µm)', 'dust size (µm)', 'dust size standard deviation (µm)',
                                 '''mixtures ratios (TiO2, Mg2SiO4, SiO, SiO2, Fe, Al2O3, CaTiO3, FeO, FeS, Fe2O3, 
                                    MgO, MgSiO3, CaSiO3, Fe2SiO4, C, KCl) (UNITLESS)''']
group.attrs['columns output'] = ['scattering coefficient to cross section ratio (UNITLESS)',
                                 'absorption/extinction coefficient to cross section ratio (UNITLESS)']

fp.close()
