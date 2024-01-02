#    -------------------------------------------------------------------------- 
#    Zhejiang Lab, Zhejiang, China
#    Yanfei Tang, tangyf@zhejianglab.com

#    Copyright (2023) Zhejiang Lab. This work is supported under the terms of 
#    contract 2021PB0AC02 with Zhejiang Lab. 
#    This software is distributed under the GNU General Public License.
#    -------------------------------------------------------------------------- 

import os
import younglaplace as yl
import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.patches import Wedge
from shutil import copyfile
import h5py
import argparse

parser = argparse.ArgumentParser(description=\
    "Plot force and menisci profile figures from theoretical data")
parser.add_argument("-ifile", type=str, default="theory_data.h5", \
    help = "input HDF5 file", required=False)
parser.add_argument("-oforcefig", type= str, default="force.pdf", \
    help = "Output to a force displacement PDF file")
parser.add_argument("-omfig", type= str, default="menisci.pdf", \
    help = "Output to a menisci density profile PDF file")

args = parser.parse_args()
ifile = args.ifile
oforcefig = args.oforcefig
omfig = args.omfig
print(f"Input HDF5 file: {ifile}")

# Step(1): Process HDF5 data
f = h5py.File(ifile, 'r')

# Get the menisci density profile plots
group = f["target"]
content = dict()
d = group["d"][()]
x = group["x"][:]
y = group["y"][:]

force = []
deltaz = []

for i in range(len(f.keys()) - 1):
    group = f[f"{i}"]
    force.append(group["force"][()])
    deltaz.append(group["deltaz"][()])

f.close()

fig = plt.figure(figsize=(8, 6))

ax1 = fig.add_subplot(111, aspect = 'equal', ylim = (1.0, 9.0), xlim = (0, 5))
e1 = Wedge((0, d + 1), 1, theta1 = 0.0, theta2 = 180.0, \
fill = True, fc = "g", ec = "g", lw = 0.0)
e2 = Wedge((0, d + 1), 1, theta1 = 180.0, theta2 = 360.0, \
fill = True, fc = "r", ec = "r", lw = 0.0)
ax1.add_patch(e1)
ax1.add_patch(e2)
ax1.plot( x, y, 'C0', linewidth = 4.0)
ax1.set_xlabel(r"$r/R$", fontsize = 24)
ax1.set_ylabel(r"$z/R$", fontsize = 24)
# ax1.set_title(r"$\psi = {0:.1f}, \theta = {1:.2f},d = {2:.2f},V = {3:.1f}$"\
# .format(model[i]["psi"]/np.pi*180.0, model[i]["theta1"]/np.pi*180.0, model[i]["d"], model[i]["V"]))

plt.savefig(omfig)
plt.close()

fig = plt.figure(figsize=(8, 6))
ax2 = fig.add_subplot(111, ylim = (-1, 1.25), xlim = (-2.5, 2.5))

ax2.plot(np.array(deltaz), np.array(force) )

ax2.set_xlabel(r"$\Delta z/R$", fontsize = 24)
ax2.set_ylabel(r"$\frac{F}{2 \pi \gamma R}$", fontsize = 24)
plt.tight_layout()
plt.savefig(oforcefig)

plt.close()