#!/usr/bin/python

#import scipy
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import re
# from math import pi
# from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
# import sys

# rc('text', usetex=True)

def crossWLF(T, shr, p):
    n    = 0.252 #1
    Tau  = 30800.0 #Pa
    D1   = 4.76e10 #Pa s
    D2   = 373.15 #K
    D3   = 0.51e-6 #K
    A1   = 25.7 #1
    A2   = 61.06 #K
    Tz   = D2
    eta0 = D1 * np.exp( -A1 * (T-Tz-D3*p) /(A2 + (T-Tz)) )
    eta  = eta0 / ( 1 + np.power( eta0 * shr / Tau , 1-n) )
    return eta


def main():

    Tc = 200 # degC
    p = 60e6 # Pa
    L = 0.1e-3 # m
    W = 1e-4 # m
    Thc = 1e-3 # m
    U = 0.1 # m/s
    T = Tc + 273.15 # K
    A0 = L*W
    shr = U/Thc

    forcesFileName = "./postProcessing/forces/0/forces.dat"

    forcesFileContent = open(forcesFileName).read()
    forcesFileContent = re.sub('[()]', '', forcesFileContent)

    open(forcesFileName, 'w').write(forcesFileContent)

    forces = np.loadtxt(forcesFileName, skiprows=3)

    time = forces[:, 0]
    force = forces[:, 5]

    fig, ax = plt.subplots(1)

    labels = ["$F_{visc}$, final [N]: %.4E" % force[-1] ]

    plt.xlabel("Time [s]")
    plt.ylabel("Force [N]")
    plt.plot(time, force, color='red')
    plt.legend(labels, loc="center right")
    plt.savefig('forceInTime.png')

    print('Viscosity in the model:', force[-1]/(A0 * shr), 'Pa s')
    viscP = crossWLF(T, shr, p)
    print('Expected viscosity:    ', viscP, 'Pa s')

if __name__ == "__main__":
    main()
