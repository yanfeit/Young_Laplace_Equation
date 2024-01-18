#    -------------------------------------------------------------------------- 
#    Zhejiang Lab, Zhejiang, China
#    Yanfei Tang, tangyf@zhejianglab.com

#    Copyright (2023) Zhejiang Lab. This work is supported under the terms of 
#    contract 2021PB0AC02 with Zhejiang Lab. 
#    This software is distributed under the GNU General Public License.
#    -------------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py

docstr = """
Compare simualtional and theoretical data.
"""

parser = argparse.ArgumentParser(description=\
    "Compare simualtional and theoretical data.")
parser.add_argument("-ifile1", type=str, default="theory_data.h5", \
    help = "Theoretical HDF5 data", required=False)
parser.add_argument("-ifile2", type= str, default="simulation_force.h5", \
    help = "Simulational force displacement HDF5 data", required=False)
parser.add_argument("-ofile", type = str,default="force.pdf", \
    help = "Output to a figure file" )

args = parser.parse_args()
thydataname = args.ifile1
sltdataname = args.ifile2
ofilename   = args.ofile

print(f"Input theoretical dataset : {thydataname}")
print(f"Input simulational dataset : {sltdataname}")

# Load Data from simulation
with h5py.File(sltdataname, 'r') as f:
    sltdata = f["force"][:]

# Load Data from theory
force, deltaz, psi, deltazp = [], [], [], []
with h5py.File(thydataname, 'r') as f:
    groupt = f["target"]
    D = groupt["D"][()]
    L = groupt["L"][()]
    R = groupt["R"][()]
    h0 = groupt["h0"][()]

    for i in range(len(f.keys()) - 1):
        group = f[f"{i}"]
        force.append (group["force"][()])
        deltaz.append(group["deltaz"][()])
        deltazp.append(group["deltazp"][()])
        psi.append   (group["psi"][()])


# surface tension, Plot the Simulational data
gamma = 1.018 
exp_x = sltdata[:, 0]/R - h0
exp_F = -sltdata[:, 1]/2/np.pi/gamma/R
exp_err = sltdata[:, 2]/2/np.pi/gamma/R

plt.errorbar(exp_x, exp_F, exp_err, fmt = "o", capsize = 2.0, fillstyle = 'none', color = 'k', label = "Simulation")

# Categorize the data, plot the theoretical data
force  = np.array(force)
deltaz = np.array(deltaz)
psi    = np.array(psi)

displacement1, displacement2, displacement3, force1, force2, force3 = [], [], [], [], [], []
for i in range(len(force)):
    if psi[i] >= np.pi/2.0 + 0.01:
        displacement1.append(deltazp[i])
        force1.append(force[i])
    elif psi[i] > (np.pi/2.0 - 0.01) and psi[i] < (np.pi/2.0 + 0.01):
        displacement2.append(deltazp[i])
        force2.append(force[i])
    elif psi[i] <= (np.pi/2.0 - 0.01):
        displacement3.append(deltazp[i])
        force3.append(force[i])

plt.plot(displacement1, force1, '-', color = "g")
plt.plot(displacement2, force2, '-', color = "b", label = "Theory")
plt.plot(displacement3, force3, '-', color = "r")

# plt.axhline(y = 0.0)
# plt.axvline(x = 0.0)

plt.xlabel(r"$\Delta z'/R$", fontsize = 24)
plt.ylabel(r"$\frac{F}{2\pi \sigma R}$", fontsize  =24)

plt.xticks(fontsize = 16)
plt.yticks(fontsize = 16)

plt.legend(loc='upper left', frameon = False)

plt.tight_layout()
plt.savefig(ofilename)
# plt.show()
