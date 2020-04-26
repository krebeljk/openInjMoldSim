import numpy as np
from io import StringIO
import matplotlib.pyplot as plt

# get data
f = open('constant/cpTable')
bes = f.read()
f.close()
pod = bes.replace("("," ").replace(")", " ")
tabela = np.loadtxt(StringIO(pod), usecols=(0,2))
temp = np.add(tabela[:,0], -273.15)
cp = tabela[:,1]

# plot
fig, ax = plt.subplots()
ax.plot(temp, cp)
ax.set(xlabel='T (degC)', ylabel='c_p (J/kgK)',
       title='Specific heat at constant pressure')
ax.grid()

fig.savefig("cp.png")
