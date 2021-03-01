import sys
import pandas as pd
from cdwave import data
sys.path.append('.')

print('Reading parameter file ...')
parameter_df = pd.read_csv('cdwave/parameters.csv', index_col=0)
name_dict = parameter_df['Name'].to_dict()
anno_dict = parameter_df['Annotation'].to_dict()
print('There are {} parameters putting to docs/support_parameters.rst'.format(len(parameter_df)))
strings = ""
for var in anno_dict:
    name = name_dict[var]
    anno = anno_dict[var]
    line = "*\t{} (`{}`): {}\n".format(name, var, anno)
    strings += line
strings += ''
with open('docs/support_parameters.rst', 'wb') as fp:
    fp.write(strings.encode("utf-8"))
