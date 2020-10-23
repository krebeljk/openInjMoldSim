import numpy as np
from io import StringIO
import matplotlib.pyplot as plt
import re
import sys
'''
example usage:
python plot_mm.py postProcessing/minMaxp/0/fieldMinMax.dat p_rgh
'''
# get data
cas = ''
fieldMin= ''
fieldMax = ''
fieldName = str(sys.argv[2])
with open(str(sys.argv[1])) as origin_file:
    data = origin_file.read()
    for match in re.findall(r'(?m).*'+fieldName+".*", data):
        cas = cas + match.split()[0] + '\n'
        fieldMin = fieldMin  + match.split()[2] + '\n'
        fieldMax = fieldMax  + match.split()[7] + '\n'

cas = np.loadtxt(StringIO(cas))
fieldMin = np.loadtxt(StringIO(fieldMin))
fieldMax = np.loadtxt(StringIO(fieldMax))

# plot
fig, ax = plt.subplots()

ax.set(xlabel='t [s]'
      ,ylabel='Min/max'
      ,title=fieldName)
ax.grid()

ax.plot(cas, fieldMin, 'b', label='min')
ax.plot(cas, fieldMax, 'r', label='max')

# legend = ax.legend()

plt.show()
# fig.savefig("kappa.png")
