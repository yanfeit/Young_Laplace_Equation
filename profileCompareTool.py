#    -------------------------------------------------------------------------- 
#    Zhejiang Lab, Zhejiang, China
#    Yanfei Tang, tangyf@zhejianglab.com

#    Copyright (2023) Zhejiang Lab. This work is supported under the terms of 
#    contract 2021PB0AC02 with Zhejiang Lab. 
#    This software is distributed under the GNU General Public License.
#    -------------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.patches import Wedge
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
parser.add_argument("-ylimn", type= int, default=3, \
    help = "Y axis limit in the neg direction")

args  = parser.parse_args()
sln = args.sln
thy = args.thy
ofile = args.ofile
ylim_neg = args.ylimn

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
    r         = f.get('r')[:]
    z         = f.get('z')[:]

with h5py.File(thy, 'r') as f:
    groupt = f["target"]
    x = groupt["x"][:]
    y = groupt["y"][:]
    R = groupt["R"][()]
    D = groupt["D"][()]
    d = groupt["d"][()]
    psi = groupt["psi"][()]


# Draw the Plots
fig = plt.figure(figsize=(6, 6))
ax  = fig.add_subplot()

ax.pcolormesh(r/R, z/R, rho, cmap = "Blues"\
    ,linewidth=0,rasterized=True)
ax.plot(rdata[0:-1:2]/R, interface[0:-1:2]/R, marker = "o", fillstyle = "none", color = 'k'\
    ,ls = '')
ax.plot(x, y, 'r', lw = 2)

e1 = Wedge((0, d + 1), 1, theta1 = 0.0, theta2 = 180.0, \
fill = True, fc = "orange", ec = "orange", lw = 0.0)
e2 = Wedge((0, d + 1), 1, theta1 = 180.0, theta2 = 360.0, \
fill = True, fc = "cyan", ec = "cyan", lw = 0.0)

ax.add_patch(e1)
ax.add_patch(e2)

# ylim_neg = y[-1]-1
xlim_pos = r[-1]/R 
ylim_pos = ylim_neg + xlim_pos 

psi_deg = psi /np.pi * 180.0

ax.text(xlim_pos-2, ylim_pos-1,  \
 r"$\psi = {0:.1f} \degree$".format(psi_deg), fontsize = 24)
ax.text(xlim_pos-2, ylim_pos-1.5,  \
 r"$D/R={0:.2f}$".format(d), fontsize = 24)


ax.set_aspect('equal',adjustable='box')
ax.set_xlim([0, xlim_pos])
ax.set_ylim([ylim_neg, ylim_pos])


plt.xlabel(r"r/R", fontsize = 32)
plt.ylabel(r"z/R", fontsize = 32)

plt.xticks(fontsize=32)
plt.yticks(fontsize=32)

fig.tight_layout()
plt.savefig(ofile)
plt.close()