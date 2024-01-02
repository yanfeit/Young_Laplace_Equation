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
"Young-Laplace Equation of Process of Particle Detachment")
parser.add_argument("-ifile", type=str, default="theory_homo_data.h5", \
    help = "input HDF5 file", required=False)
parser.add_argument("-isFull", type=bool, default=True, \
    help = "Need both plot figures and generate data?", required=False)
parser.add_argument("-ofile", type= str, default="janus.mp4", \
    help = "Output to a mp4 file")

args = parser.parse_args()
ifile = args.ifile
isFull = args.isFull
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
    model.append(content)

f.close()

if isFull:
    # Step (2) : Initialize some memory for figures
    displacement = []
    force        =  []

    # Step (3) : Plot figures
    for i in range(0, len(model), 2):
        displacement.append( model[i]["deltaz"] )
        force.append( model[i]["force"] )
        
        fig = plt.figure(figsize=(12, 4.5))

        ax1 = fig.add_subplot(121, aspect = 'equal', ylim = (1.0, 9.0), xlim = (-5, 5))
        e1 = Wedge((0, model[i]["d"] + 1), 1, theta1 = 0.0, theta2 = 180.0, \
        fill = True, fc = "b", ec = "b", lw = 0.0)
        e2 = Wedge((0, model[i]["d"] + 1), 1, theta1 = 180.0, theta2 = 360.0, \
        fill = True, fc = "b", ec = "b", lw = 0.0)
        ax1.add_patch(e1)
        ax1.add_patch(e2)
        ax1.plot( model[i]["x"], model[i]["y"], 'C0', linewidth = 4.0)
        ax1.plot(-model[i]["x"], model[i]["y"], 'C0', linewidth = 4.0)
        ax1.set_xlabel(r"$r/R$", fontsize = 24)
        ax1.set_ylabel(r"$z/R$", fontsize = 24)
        ax1.set_title(r"$\psi = {0:.1f}, \theta = {1:.2f},d = {2:.2f},V = {3:.1f}$"\
        .format(model[i]["psi"]/np.pi*180.0, model[i]["theta1"]/np.pi*180.0, model[i]["d"], model[i]["V"]))

        ax2 = fig.add_subplot(122, ylim = (-0.5, 1.25), xlim = (-1.0, 3.0))

        ax2.plot(displacement, force, 'o-', color = "b", fillstyle = "none")

        ax2.set_xlabel(r"$\Delta z/R$", fontsize = 24)
        ax2.set_ylabel(r"$\frac{F}{2 \pi \gamma R}$", fontsize = 24)

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