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
Collect and analyze data from LAMMPS simulator.
"""

def blocking(x, n):
    size = len(x)
    b = int(size/n)
    data = np.zeros(n)
    for i in range(n):
        data[i] = sum(x[b*i: b*(i+1)])
        data[i] /= b

    return data.mean(), data.std()


parser = argparse.ArgumentParser(description=\
    "Collect and anlyaze force data from LAMMPS simulator")
parser.add_argument("-idir", type=str, default="./", \
    help = "directory folder", required=False)
parser.add_argument("-ofile", type= str, default="simulation_force.h5", \
    help = "Output to a force displacement HDF5 data")
parser.add_argument("-height", type = int, default=25, \
    help= "Initial center of the particle")
parser.add_argument("-endh", type = int, default=62, \
    help= "Lift height for the particle, add discprancy")

args = parser.parse_args()
idir = args.idir
ofile = args.ofile
height = args.height
endh  = args.endh

print(f"Input Directory: {idir}")

# height   = 25
aveForce, stdForce = [], []

# Collect data from ferg.txt files in each directory.
for i in range(2, endh, 2):
    # I use pos{i} for storing different measurements.
    dirname   = idir +  "pos{0}/".format(i)
    ifilename = dirname + "ferg.txt"
    ifile     = np.loadtxt(ifilename, skiprows = 2)
    force, force_std = blocking(ifile[:, 4], 10)
    aveForce.append(force)
    stdForce.append(force_std)

h = height + np.array(range(2, endh, 2))
aveForce = np.array(aveForce)
stdForce = np.array(stdForce)

dataset = np.vstack((h, aveForce, stdForce)).transpose()

with h5py.File(ofile, 'w') as f:
    f.create_dataset("force", data = dataset)
