import numpy as np
import matplotlib.pyplot as plt
import os

"""
to do list: 1.autocorrlation!
2. F vs h
3. air/liquid interface
4. small angle!
5. 
"""

def blocking(x, n):
    size = len(x)
    b = int(size/n)
    data = np.zeros(n)
    for i in range(n):
        data[i] = sum(x[b*i: b*(i+1)])
        data[i] /= b

    return data.mean(), data.std()


height   = 25
aveForce = []
stdForce = []
for i in range(2, 62, 2):
    dirname   = "./pos{0}/".format(i)
    ifilename = dirname + "ferg.txt"
    ifile     = np.loadtxt(ifilename, skiprows = 2)
    force, force_std = blocking(ifile[:, 4], 10)
    aveForce.append(force)
    stdForce.append(force_std)

h = np.linspace(27, 85, 30)
aveForce = np.array(aveForce)

with open("Fvsh_Ans100.txt", 'w') as ofile:
    ofile.write("H \t Force \t std\n")
    for i in range(len(h)):
        ofile.write(str(h[i]) + "\t" + str(aveForce[i]) + "\t" + str(stdForce[i]) + "\n")





# plt.plot(h, aveForce,'o-')
# plt.xlabel(r"$z/\sigma$", fontsize= 16)
# plt.ylabel(r"$F_{colloid-particles}^{z}/(\epsilon/\sigma)$", fontsize = 16)
# plt.savefig("aveForce.jpg")
# plt.show()




