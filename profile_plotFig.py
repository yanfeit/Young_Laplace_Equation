#!/bin/python
# Yanfei Tang

import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py

hdfFile = h5py.File(sys.argv[1], 'r')
rho = hdfFile.get('rho')

interface = np.array(hdfFile.get('interface'))
rdata = np.array(hdfFile.get('rdata'))

theory_curve = np.loadtxt("theoretical_x_y.txt")


# data = np.loadtxt(sys.argv[1])
# rho = data[1:, 1:]

binsize = 1.0
span = (50.0, 100.0) # (r, z)

z, r = np.mgrid[slice(0, span[1] + binsize, binsize),
                slice(0, span[0] + binsize, binsize)]

plt.pcolormesh(r, z, rho, cmap = "Blues")

plt.plot(rdata, interface, marker = "o", fillstyle = "none")

plt.plot(theory_curve[:,0], theory_curve[:,1], 'r')
plt.axis("equal")

plt.savefig("test.png")