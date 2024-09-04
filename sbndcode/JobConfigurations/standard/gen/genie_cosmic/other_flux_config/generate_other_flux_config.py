'''
This script must be run from the sbndcode/JobConfiguration/standard/gen/genie_cosmic folder
'''

config = 'configf'
flux_table = 'sbnd_flux_bnb_nu_Fv1'

import os
from glob import glob

fcl_files = glob('../*.fcl')

directory = 'flux_' + config

if os.path.isdir(directory):
	print(f'Directory {directory} already exists. Exiting.')
	exit()

print(f'Creating directory {directory}.')
os.mkdir(directory)

for fcl_file in fcl_files:

	new_name = directory + '/' + fcl_file[3:].split('_sbnd.fcl')[0] + '_' + config + '_sbnd.fcl'

	print('Creating fcl file:', new_name)

	with open(new_name, "w") as new_fcl_file:

		new_fcl_file.write(f'#include "{fcl_file[3:]}" \n')
		new_fcl_file.write('\n')
		new_fcl_file.write('physics.producers.generator: {\n')
		new_fcl_file.write('    @table::physics.producers.generator\n')
		new_fcl_file.write(f'    @table::{flux_table}\n')
		new_fcl_file.write('}\n')

