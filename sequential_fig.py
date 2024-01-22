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
    model.append(content)

f.close()


displacement1, displacement2, displacement3 = [], [], []
force1, force2, force3 = [], [], []


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


