import numpy as np
from io import StringIO
import matplotlib.pyplot as plt
import re
import sys
'''
example usage:
python plot_iters.py log.openInjMoldSim
'''
# get data
cas = ''
iters = ''
with open(str(sys.argv[1])) as origin_file:
    data = origin_file.read()
    for match in re.findall(r'(?m)^Time\s=.*', data):
        cas = cas + match.split('=')[1] + '\n'
    for match in re.findall(r'(?m)^PIMPLE:.*in.*iterations', data):
        iters = iters + match.split(' ')[-2] + '\n'

cas = np.loadtxt(StringIO(cas))
iters = np.loadtxt(StringIO(iters))

if cas.size == iters.size+1:
    cas = cas[1:]

# plot
fig, ax = plt.subplots()

ax.set(xlabel='t [s]'
      ,ylabel='iters [1]'
      ,title='Pimple iterations')
ax.grid()

ax.plot(cas, iters, 'k', label='iters')

# legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')

plt.show()
# fig.savefig("iters.png")
