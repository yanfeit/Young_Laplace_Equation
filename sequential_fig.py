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

