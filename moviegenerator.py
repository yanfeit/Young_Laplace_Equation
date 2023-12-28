#    -------------------------------------------------------------------------- 
#    Zhejiang Lab, Zhejiang, China
#    Yanfei Tang, tangyf@zhejianglab.com

#    Copyright (2023) Zhejiang Lab.  Under the terms of Contract
#    2021PB0AC02 with Zhejiang Lab, the Zhejiang Government retains
#    certain rights in this software.  This software is distributed under
#    the GNU General Public License.
#    -------------------------------------------------------------------------- 


import younglaplace as yl
import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.patches import Wedge
from shutil import copyfile
import h5py
import argparse

parser = argparse.ArgumentParser(description=
"Young-Laplace Equation of Process of Particle Detachment")
parser.add_argument("ifile", type=str, default="theory_data.h5", 
help = "input HDF5 file")

args = parser.parse_args()
ifile = args.ifile
print(f"Input HDF5 file: {ifile}")

# Process HDF5 data
with h5py.File(ifile, 'r') as f:
    

data = np.loadtxt("info.txt", skiprows = 1)

R = data[:, 0]
L = data[:, 1]
D = data[:, 2]
theta1 = data[:, 3]
psi = data[:, 4]

l = L[0]/R[0]
d = D[0]/R[0]
height = 50.0

interface = 0.5 * 4.0/3.0 * np.pi/l/l/np.pi*R[0] + height

displacement1 = []
displacement2 = []
displacement3 = []
force1 = []
force2 = []
force3 = []
sign = 1

for i in range(0, len(R), 2):

    model = yl.YL(R = R[i], L = L[i], D = D[i], theta1 = theta1[i], psi = psi[i])
    if psi[i] > 90:
        displacement1.append( (model.D + R[0] - interface)/R[0] )
        force1.append( model.force )
        sign = 1
    elif int(psi[i]) ==  90:
        displacement2.append( (model.D + R[0] - interface)/R[0] )
        force2.append( model.force )
        sign = 2
    elif psi[i] < 90:
        displacement3.append( (model.D + R[0] - interface)/R[0] )
        force3.append( model.force )
        sign = 1

    fig = plt.figure(figsize=(12, 4.5))

    ax1 = fig.add_subplot(121, aspect = 'equal', ylim = (1.0, 9.0), xlim = (-5, 5))
    e1 = Wedge((0, model.d + 1), 1, theta1 = 0.0, theta2 = 180.0, fill = True, fc = "g", ec = "g", lw = 0.0)
    e2 = Wedge((0, model.d + 1), 1, theta1 = 180.0, theta2 = 360.0, fill = True, fc = "r", ec = "r", lw = 0.0)
    ax1.add_patch(e1)
    ax1.add_patch(e2)
    ax1.plot( model.x, model.y, 'C0', linewidth = 4.0)
    ax1.plot(-model.x, model.y, 'C0', linewidth = 4.0)
    ax1.set_xlabel(r"$r/R$", fontsize = 24)
    ax1.set_ylabel(r"$z/R$", fontsize = 24)
    ax1.set_title(r"$\psi = {0:.1f}, \theta = {1:.2f},d = {2:.2f},V = {3:.1f}$".format(model.psi/np.pi*180.0, model.theta1/np.pi*180.0, model.d, model.V))

    ax2 = fig.add_subplot(122, ylim = (-1, 1.25), xlim = (-2.5, 2.5))

    ax2.plot(displacement1, force1, 'o-', color = "g", fillstyle = "none")
    ax2.plot(displacement2, force2, 'o-', color = "b", fillstyle = "none")
    ax2.plot(displacement3, force3, 'o-', color = "r", fillstyle = "none")

    ax2.set_xlabel(r"$\Delta z/R$", fontsize = 24)
    ax2.set_ylabel(r"$\frac{F}{2 \pi \gamma R}$", fontsize = 24)
    if sign == 1:
        ax2.text(-2, 0.5, 'sliding', fontsize = 36)
    elif sign == 2:
        ax2.text(-2, 0.5, 'pinned', fontsize = 36)
    


    plt.tight_layout()
    plt.savefig("theta1{0:.1f}psi{1:.1f}.png".format(model.theta1/np.pi*180.0, model.psi/np.pi*180.0), transparent = True)
    plt.clf()


