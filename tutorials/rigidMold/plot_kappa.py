import numpy as np
from io import StringIO
import matplotlib.pyplot as plt

# get data
f = open('constant/kappaTable')
bes = f.read()
f.close()
pod = bes.replace("("," ").replace(")", " ")
tabela = np.loadtxt(StringIO(pod), usecols=(0, 2, 4, 6))
temp = np.add(tabela[:,0], -273.15)
kappa1 = tabela[:,1]
kappa2 = tabela[:,2]
kappa3 = tabela[:,3]

# plot
fig, ax = plt.subplots()

ax.set(xlabel='T (degC)'
    ,ylabel='kappa (W/mK)'
    ,title='Thermal conductivity')
ax.grid()

ax.plot(temp, kappa1, 'k', label='20 MPa')
ax.plot(temp, kappa2, 'k--', label='80 MPa')
ax.plot(temp, kappa3, 'k:', label='120 MPa')

legend = ax.legend(loc='upper right', shadow=True, fontsize='x-large')

# plt.show()
fig.savefig("kappa.png")
