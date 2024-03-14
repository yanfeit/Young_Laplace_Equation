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

parser = argparse.ArgumentParser(description=
"Plot the Process of Janus Particle Detachment")
parser.add_argument("-ifile", type=str, default="theory_data.h5", \
    help = "input HDF5 file", required=False)
parser.add_argument("-ofile", type= str, default="sequential_fig.pdf", \
    help = "Output to a PDF file", required=False)

args = parser.parse_args()
ifile = args.ifile
ofile = args.ofile
print(f"Input HDF5 file: {ifile}")

# Step(1): Process HDF5 data
model = []
f = h5py.File(ifile, 'r')
for i in range(len(f.keys()) - 1):
    group = f[f"{i}"]
    content = dict()
    content["R"] = group["R"][()]
    content["D"] = group["D"][()]
    content["L"] = group["L"][()]
    content["d"] = group["d"][()]
    content["l"] = group["l"][()]
    content["theta1"] = group["theta1"][()]
    content["psi"] = group["psi"][()]
    content["x"] = group["x"][:]
    content["y"] = group["y"][:]
    content["force"] = group["force"][()]
    content["V"] = group["V"][()]
    content["deltaz"] = group["deltaz"][()]
    content["regime"] = group["regime"][()]
    content["deltazp"] = group["deltazp"][()]
    model.append(content)

f.close()


displacement1, displacement2, displacement3 = [], [], []
force1, force2, force3 = [], [], []


for i in range(0, len(model), 2):
    if model[i]["psi"]/np.pi*180.0 > 90:
        displacement1.append( model[i]["deltazp"] )
        force1.append( model[i]["force"] )
    elif int(model[i]["psi"]/np.pi*180.0) == 90:
        displacement2.append( model[i]["deltazp"] )
        force2.append( model[i]["force"] )
    elif model[i]["psi"]/np.pi*180.0 < 90 :
        displacement3.append( model[i]["deltazp"] )
        force3.append( model[i]["force"] )


fig = plt.figure(figsize=(10, 4.5))

ax1 = fig.add_subplot(121, aspect = 'equal', adjustable = 'box', ylim = (2.0, 7.5), xlim = (0, model[i]["x"][-1]))

# ax1.set_facecolor("#222")

ct = len(model)

step = 60
for i in range(0, len(model), step):
    col = ((ct-i)/1.25/ct, (ct-i)/ct/1.25, (ct-i)/ct/1.25)
    # col = (0.0, 0.0, 0.0) # pure black
    print(col)
    ax1.plot(model[i]["x"], model[i]["y"],color = col, lw=2.0, label = 'State {}'.format(int(i/step)))

    e1 = Wedge((0, model[i]["d"] + 1), 1, theta1 = 0.0, theta2 = 180.0, \
    fill = False, fc = "None", ec = col, lw = 2.0)
    e2 = Wedge((0, model[i]["d"] + 1), 1, theta1 = 180.0, theta2 = 360.0, \
    fill = False, fc = "None", ec = col, lw = 2.0)

    ax1.add_patch(e1)
    ax1.add_patch(e2)


ax1.set_xlabel(r"$r/R$", fontsize = 24)
ax1.set_ylabel(r"$z/R$", fontsize = 24)

ax1.tick_params('both', labelsize=16)

ax1.legend(loc='lower right', frameon=False)

ax1.annotate('(a)', xy=(0.8, 0.9), xycoords='axes fraction', fontsize=24)

ax2 = fig.add_subplot(122, ylim = (-1, 1.25), xlim = (-2.2, 2.2))

ax2.plot(displacement1, force1, '-', color = "g", fillstyle = "none", lw=2)
ax2.plot(displacement2, force2, '-', color = "b", fillstyle = "none", lw=2)
ax2.plot(displacement3, force3, '-', color = "r", fillstyle = "none", lw=2)

ax2.tick_params('both', labelsize=16)
ax2.set_xlabel(r"$\Delta z'/R$", fontsize = 24)
ax2.set_ylabel(r"$\frac{F}{2 \pi \gamma R}$", fontsize = 24)

ax2.annotate('(b)', xy=(0.8, 0.9), xycoords='axes fraction', fontsize=24)
# if model[i]["regime"] == 1:
#     ax2.text(-2, 0.5, 'sliding', fontsize = 36)
# elif model[i]["regime"] == 2:
#     ax2.text(-2, 0.5, 'pinned', fontsize = 36)


ax2.annotate('Sliding', xy=(0.2, 0.09),  xycoords='axes fraction', fontsize=16)
ax2.annotate('Pinned',  xy=(0.52, 0.4),  xycoords='axes fraction', fontsize=16)
ax2.annotate('Sliding', xy=(0.75, 0.69), xycoords='axes fraction', fontsize=16)


plt.tight_layout()
plt.savefig(ofile)
plt.clf()
plt.close()