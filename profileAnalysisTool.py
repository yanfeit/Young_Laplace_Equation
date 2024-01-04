#    -------------------------------------------------------------------------- 
#    Zhejiang Lab, Zhejiang, China
#    Yanfei Tang, tangyf@zhejianglab.com

#    Copyright (2023) Zhejiang Lab. This work is supported under the terms of 
#    contract 2021PB0AC02 with Zhejiang Lab. 
#    This software is distributed under the GNU General Public License.
#    -------------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
import dump as dp
import sys
from scipy.optimize import curve_fit
import h5py
import argparse

def func(r, r0, ds):
    rho_l = 0.927
    rho_g = 0.000
    rho   = 0.5 * (rho_l + rho_g) - 0.5 * (rho_l - rho_g) * np.tanh(2*(r - r0)/ds)
    return rho

def dis(x, y, x0, y0):
    dx = x - x0
    dy = y - y0
    return np.sqrt(dx*dx + dy*dy)

def ringvol(r0, r1, h):
    """
    inner radius : r0
    outer radius : r1
    height : h
    """
    dS = np.pi * (r1 * r1 - r0 * r0)
    return dS * h

parser = argparse.ArgumentParser(description=\
    "Collect and analyze data from LAMMPS dump to get the density profiles.")
parser.add_argument("-idump", type=str, default="./pos60/dump_janus_measure.1", \
    help = "dump file", required=False)
parser.add_argument("-ofile", type= str, default="simulation_profile.h5", \
    help = "Output to a HDF5 data set") 




d = dp.dump(sys.argv[1])
binsize = 1.0

# Dimension calculation
axis = (50.0, 50.0)
span = (50.0, 100.0)
Nr, Nz = int(span[0]/binsize), int(span[1]/binsize)
rho = np.zeros((Nz, Nr), dtype= 'float64')
vol = np.zeros(Nr, dtype = "float64")

for i in range(len(vol)):
    vol[i] = ringvol(binsize * i, binsize * (i+1), binsize)

for i in range(d.nsnaps):

    time, box, atoms, bonds, tris, c = d.viz(i)
    print(f"Processing Snap {i}")

    for atom in atoms:
        
        if int(atom[1]) == 1: 
            r = dis(atom[2], atom[3], axis[0], axis[1])
            z = atom[4]

            nr = int(r/binsize)
            nz = int(z/binsize)
            rho[nz, nr] += 1

rho = rho/vol/d.nsnaps

# Data Output

for i in range(Nr):
    rticks = (i + 0.5) * binsize

for i in range(Nz):
    zticks = (i + 0.5) * binsize

rdata  = np.arange(binsize/2.0, span[0] + binsize/2.0, binsize)
zdata  = np.arange(binsize/2.0, span[1] + binsize/2.0, binsize)
interface = []

for i in range(Nr):
    popt, pcov = curve_fit(func, zdata[4:-4], rho[4:-4, i])
    interface.append(popt[0])


# rho, 二维矩阵, nz=100, nr=50, 
# rdata, r方向上的坐标，长度=50
# interface, 气液界面的z坐标， 
# rticks, r方向的？
# zticks，z方向的？

hdfFile = h5py.File(sys.argv[2], 'w')
hdfFile.create_dataset("rho", data = rho)
hdfFile.create_dataset("rdata", data = rdata)
hdfFile.create_dataset("interface", data = interface)
hdfFile.create_dataset("rticks", data = rticks)
hdfFile.create_dataset("zticks", data = zticks)

hdfFile.close()