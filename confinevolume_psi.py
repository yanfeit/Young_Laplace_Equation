#!/bin/python
# Yanfei Tang 
# Date: 2018, 2023 @ Virginia Tech, Zhejiang Lab

import younglaplace as yl
import numpy as np
import matplotlib.pyplot as plt
# from   matplotlib.patches import Circle
# import scipy.special as spe
# import sys

docstr = """
In the molecular dynamics simulation, the volume of the liquids is a constant. During the 
pulling process, the distance D between the bottom of the sphere and cylindrical container
increases. However, the script younglaplace.py computes the volumn given the distance. To 
enforce the given volumn, it is required to pre-known the distance D. 

Given a confined volume and a contact angle, compute the force v.s. displacement relationship.
This is what I called slinding region.
"""

# INPUT 
# the location of interface, 0.76 is the effective wall location
interface = 50.9  
L         = 49.7
R         = 10.35
D         = 50.0
d         = D/R
l         = L/R
theta1    =  83.655540
# reduced volume of liquid
V         = np.pi*l*l*interface/R
# End input

# The first step is we need a plot relation
# filling angle v.s. volume given a contant particle position!

psis = np.linspace(1.0, 179.0, 100)
Vs   = []
for psi in psis:
    model = yl.YL(R = R, L = L, D = D, theta1 = theta1, psi = psi)
    Vs.append(model.V)

Vs = np.array(Vs)
indexmin = np.where(Vs == Vs.min())[0][0]
indexmax = np.where(Vs == Vs.max())[0][0]

# Notice that the line of filling angle v.s volume is not 
# monotonic, here needs some stable analysis. We think that
# large filling anlge is more stable than small filing angle. Why???
# And we only consider filling in this region
# psi in [psi[indexmin], psi[indexmax]]

psi_min = psis[indexmin]
psi_max = psis[indexmax]

psis = np.linspace(psi_min, psi_max, 100)
Ds = []

for psi in psis:
    model  = yl.YL(R = R, L = L, D = D, theta1 = theta1, psi = psi)
    deltah = (V - model.V)/np.pi/l/l
    Ds.append(D + deltah*R)

Ds = np.array(Ds)
forces = []
ofile = open("YL_Ans100.txt", 'w')
ofile.write("relative distance  force\n")
for i in range(len(psis)):

    model = yl.YL(R = R, L = L, D = Ds[i], theta1 = theta1, psi = psis[i])
    forces.append(model.force)
    ofile.write("{0} {1}\n".format((Ds[i]+R-interface)/R, model.force))
    print((Ds[i]+R - interface)/R, model.force, model.V)

plt.plot((Ds+R-interface)/R, forces)
plt.xlabel(r"$z/R$")
plt.ylabel(r"$\frac{F}{2 \pi \sigma R}$")

plt.tight_layout()
plt.savefig("Figure_L{0}.pdf".format(L))
plt.show()