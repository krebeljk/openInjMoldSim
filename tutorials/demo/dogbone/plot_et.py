import numpy as np
from io import StringIO
import matplotlib.pyplot as plt
import re
import sys
'''
example usage:
python plot_dt.py log.openInjMoldSimFimaaIbaa
'''
# get data
cas = ''
dcas = ''
with open(str(sys.argv[1])) as origin_file:
    data = origin_file.read()
    for match in re.findall(r'(?m)^Time\s=.*', data):
        cas = cas + match.split('=')[1] + '\n'
    for match in re.findall(r'(?m)^ExecutionTime\s=.*', data):
        dcas = dcas + match.split('=')[1] + '\n'

dcas = dcas.replace("s","")
cas = np.loadtxt(StringIO(cas))
dcas = np.loadtxt(StringIO(dcas))

if cas.size == dcas.size+1:
    cas = cas[1:]

# plot
fig, ax = plt.subplots()

ax.set(xlabel='t [s]'
      ,ylabel='Execution time [h]'
      ,title='Execution time'
)
ax.grid()

ax.plot(cas, dcas/3600.0, 'k', label='dt')

# legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')

plt.show()
# fig.savefig("kappa.png")
