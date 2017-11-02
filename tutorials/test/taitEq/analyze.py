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

def tait(T, p):
    b1m    =    9.76e-4
    b2m    =    5.8e-7
    b3m    =    1.67e8
    b4m    =    3.6e-3
    b1s    =    9.76e-4
    b2s    =    2.3e-7
    b3s    =    2.6e8
    b4s    =    3.0e-3
    b5     =    373.0
    b6     =    5.1e-7
    b7     =    0.0
    b8     =    0.0
    b9     =    0.0
    Cuniv  =    0.0894
    if T > b5 + b6*p:
        v0 = b1m + b2m * (T-b5)
        B = b3m * np.exp(-b4m * (T-b5))
        v1 = 0.0
    else:
        v0 = b1s + b2s * (T-b5)
        B = b3s * np.exp(-b4s * (T-b5))
        v1 = b7 * np.exp(b8 * (T-b5) - b9*p)
    specV = v0 * ( 1.0 - Cuniv* np.log(1.0 + p/B)) + v1
    return specV


def main():

    pressFileName = "./postProcessing/probes/0/p"
    tmprFileName = "./postProcessing/probes/0/Tc"
    rhoFileName = "./postProcessing/probes/0/rho"

    pressFileContent = open(pressFileName).read()
    tmprFileContent = open(tmprFileName).read()
    rhoFileContent = open(rhoFileName).read()

    press = np.loadtxt(pressFileName, skiprows=3)
    tmpr = np.loadtxt(tmprFileName, skiprows=3)
    rho = np.loadtxt(rhoFileName, skiprows=3)
    n_vals = len(tmpr[:,0])
    dens = np.zeros((n_vals,2))

    for i in range(n_vals):
        dens[i,1] = 1/tait(tmpr[i,1]+273.15, press[i,1])
        dens[i,0] = tmpr[i,0]

    fig, ax = plt.subplots(1)

    labels = ["$\\rho_{py}(T_{OF}(t),p_{OF}(t))$", "$\\rho_{OF}(t)$"]
    axes = plt.gca()
    plt.xlabel("Time [s]")
    plt.ylabel("Density [kg/m3]")
    plt.plot(dens[:,0], dens[:,1], color='red')
    plt.plot(rho[:,0], rho[:,1], color='blue', linestyle='dashed')
    plt.legend(labels, loc="center right")
    plt.savefig('densityInTime.png')

if __name__ == "__main__":
    main()
