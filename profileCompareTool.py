#    -------------------------------------------------------------------------- 
#    Zhejiang Lab, Zhejiang, China
#    Yanfei Tang, tangyf@zhejianglab.com

#    Copyright (2023) Zhejiang Lab. This work is supported under the terms of 
#    contract 2021PB0AC02 with Zhejiang Lab. 
#    This software is distributed under the GNU General Public License.
#    -------------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py
import argparse

parser = argparse.ArgumentParser(description=\
    "Compare the density profile from simulatoin results and theoretical results")
parser.add_argument("-sln", type=str, default="simulation_profile.h5", \
    help = "Simulational results of HDF5 datasets", required=False)
parser.add_argument("-thy", type= str, default="theory_data.h5", \
    help = "Theoretical results of HDF5 datasets")
parser.add_argument("-ofile", type= str, default="profile_compare.pdf", \
    help = "Theoretical results of HDF5 datasets")

args  = parser.parse_args()
sln = args.sln
thy = args.thy
ofile = args.ofile

print(f"Input Simulational HDF5 datasets : {sln}")
print(f"Input Theoretical HDF5 datasets : {thy}")
print(f"Output density profile comparison figure : {ofile}")

# Fetch the Data from HDF5 datasets
with h5py.File(sln, 'r') as f:
    #
    rho       = f.get('rho')[:]
    interface = f.get('interface')[:]
    rdata     = f.get('rdata')[:]
    zdata     = f.get('zdata')[:]

with h5py.File(thy, 'r') as f:
    groupt = f["target"]
    x = groupt["x"][:]
    y = groupt["y"][:]


# Draw the Plots

z, r = np.mgrid[zdata,
                rdata]

plt.pcolormesh(r, z, rho, cmap = "Blues")

plt.plot(rdata, interface, marker = "o", fillstyle = "none")

plt.plot(x, y, 'r')
plt.axis("equal")

plt.savefig(ofile)


