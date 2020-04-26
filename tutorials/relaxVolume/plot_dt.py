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
    for match in re.findall(r'(?m)^deltaT\s=.*', data):
        dcas = dcas + match.split('=')[1] + '\n'

cas = np.loadtxt(StringIO(cas))
dcas = np.loadtxt(StringIO(dcas))

# plot
fig, ax = plt.subplots()

ax.set(xlabel='t [s]'
      ,ylabel='dt [s]'
      ,title='Casovni korak')
ax.grid()

ax.semilogy(cas, dcas, 'k', label='dt')

# legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')

plt.show()
# fig.savefig("kappa.png")
