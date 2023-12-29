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
    model.append(content)

f.close()

if isFull:
    # Step (2) : Initialize some memory for figures
    displacement1, displacement2, displacement3 = [], [], []
    force1, force2, force3 = [], [], []


    # Step (3) : Plot figures
    for i in range(0, len(model), 2):
        if model[i]["psi"]/np.pi*180.0 > 90:
            displacement1.append( model[i]["deltaz"] )
            force1.append( model[i]["force"] )
        elif int(model[i]["psi"]/np.pi*180.0) == 90:
            displacement2.append( model[i]["deltaz"] )
            force2.append( model[i]["force"] )
        elif model[i]["psi"]/np.pi*180.0 < 90 :
            displacement3.append( model[i]["deltaz"] )
            force3.append( model[i]["force"] )

        fig = plt.figure(figsize=(12, 4.5))

        ax1 = fig.add_subplot(121, aspect = 'equal', ylim = (1.0, 9.0), xlim = (-5, 5))
        e1 = Wedge((0, model[i]["d"] + 1), 1, theta1 = 0.0, theta2 = 180.0, \
        fill = True, fc = "g", ec = "g", lw = 0.0)
        e2 = Wedge((0, model[i]["d"] + 1), 1, theta1 = 180.0, theta2 = 360.0, \
        fill = True, fc = "r", ec = "r", lw = 0.0)
        ax1.add_patch(e1)
        ax1.add_patch(e2)
        ax1.plot( model[i]["x"], model[i]["y"], 'C0', linewidth = 4.0)
        ax1.plot(-model[i]["x"], model[i]["y"], 'C0', linewidth = 4.0)
        ax1.set_xlabel(r"$r/R$", fontsize = 24)
        ax1.set_ylabel(r"$z/R$", fontsize = 24)
        ax1.set_title(r"$\psi = {0:.1f}, \theta = {1:.2f},d = {2:.2f},V = {3:.1f}$"\
        .format(model[i]["psi"]/np.pi*180.0, model[i]["theta1"]/np.pi*180.0, model[i]["d"], model[i]["V"]))

        ax2 = fig.add_subplot(122, ylim = (-1, 1.25), xlim = (-2.5, 2.5))

        ax2.plot(displacement1, force1, 'o-', color = "g", fillstyle = "none")
        ax2.plot(displacement2, force2, 'o-', color = "b", fillstyle = "none")
        ax2.plot(displacement3, force3, 'o-', color = "r", fillstyle = "none")

        ax2.set_xlabel(r"$\Delta z/R$", fontsize = 24)
        ax2.set_ylabel(r"$\frac{F}{2 \pi \gamma R}$", fontsize = 24)
        if model[i]["regime"] == 1:
            ax2.text(-2, 0.5, 'sliding', fontsize = 36)
        elif model[i]["regime"] == 2:
            ax2.text(-2, 0.5, 'pinned', fontsize = 36)

        plt.tight_layout()
        plt.savefig("theta1{0:.1f}psi{1:.1f}.png".format(model[i]["theta1"]/np.pi*180.0, \
            model[i]["psi"]/np.pi*180.0), transparent = True)
        plt.clf()
        plt.close()

# Step (4) ： Generate movies
for i in range(0, len(model), 2):

    ifilename = "theta1{0:.1f}psi{1:.1f}.png".format(model[i]["theta1"]/np.pi*180.0, \
        model[i]["psi"]/np.pi*180.0)
    framename = "frame%05d.png" % int(i/2)
    copyfile(ifilename, framename)
    
os.system("ffmpeg -framerate 10 -i frame%05d.png -c:v libx264  {}".format(ofile))