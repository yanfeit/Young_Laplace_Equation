import younglaplace as yl
import numpy as np
import matplotlib.pyplot as plt
import sys
from shutil import copyfile
import os

data = np.loadtxt("info.txt", skiprows = 1)

R = data[:, 0]
L = data[:, 1]
D = data[:, 2]
theta1 = data[:, 3]
psi = data[:, 4]

for i in range(0, len(R), 2):

    ifilename = "theta1{0:.1f}psi{1:.1f}.png".format(theta1[i], psi[i])
    framename = "frame%05d.png" % int(i/2)
    
    copyfile(ifilename, framename)
    

os.system("ffmpeg -framerate 10 -i frame%05d.png -c:v libx264  janus.mp4")
#os.system("rm frame*.png")
