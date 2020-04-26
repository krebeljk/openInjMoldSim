import numpy as np
from io import StringIO
import matplotlib.pyplot as plt
import re
import sys
'''
example usage:
python plot_dt.py log.openInjMoldSimFimaaIbac
'''
# get data
cas = ''
coMean = ''
coMax = ''
with open(str(sys.argv[1])) as origin_file:
    data = origin_file.read()
    for match in re.findall(r'(?m)^Time\s=.*', data):
        cas = cas + match.split('=')[1] + '\n'
    for match in re.findall(r'(?m)^Courant\sNumber\smean:.*', data):
        coMean = coMean + match.split(' ')[3] + '\n'
        coMax = coMax + match.split(' ')[5] + '\n'

cas = np.loadtxt(StringIO(cas))
coMean = np.loadtxt(StringIO(coMean))[1:]#skip first because calculated at time 0
coMax = np.loadtxt(StringIO(coMax))[1:]

if cas.size == coMean.size+1:
    cas = cas[1:]

# plot
fig, ax = plt.subplots()

ax.set(xlabel='t [s]'
      ,ylabel='Co [1]'
      ,title='Courant number')
ax.grid()

ax.plot(cas, coMean, 'b', label='mean')
ax.plot(cas, coMax, 'r', label='max')

# legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')

plt.show()
# fig.savefig("kappa.png")
