import json
import sys
from cdwave import data
sys.path.append('.')


with open('cdwave/param_annot.json') as fp:
    anno_dict = json.load(fp)
name_dict = data.WaveformFull.parameter_annotations

strings = ""
for var in anno_dict:
    name = name_dict[var]
    anno = anno_dict[var]
    line = "*\t{} (`{}`): {}\n".format(name, var, anno)
    strings += line
strings += ''
with open('docs/support_parameters.rst', 'wb') as fp:
    fp.write(strings.encode("utf-8"))
