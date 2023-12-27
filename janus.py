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

# from   matplotlib.patches import Circle
# import scipy.special as spe
# import sys

docstr = """
The script generates a quasi-static pulling process of a Janus particle from 
the liquid. The Janus particle has a hydrophobic part of the hemisphere at the top and 
a hydrophilic part of the hemisphere at the bottom. 

This is a physical decription of the process. In the beginning of sliding region, the Janus particle is 
merged under the liquid with menisci formed. The filling angle decreases with the undergoing pulling process. 
It then pins at the equator of the sphere with a right filling angle. At some point, the process 
is back to the sliding region with a rising menisci.

input: 
hydrophobic, hydrophilic : contact angle of the two part janus particle
height :                   liquid height without considering the particle
L :                        container radius      
R :                        radius of the janus particle
D :                        distance between the bottom of the particle and bottom of the containder.

output:
HDF5
"""

class Janus(object):
    """
    hydrophobic: contact angle of hydrophobic part of the sphere
    hydrophilic: contact angle of hydrophilic part of the sphere
    """
    def __init__(self, hydrophobic = 150.0, hydrophilic = 30.0, height = 50.0, L = 50.0, R = 10.0, D = 50.0):

        self.hydrophobic = hydrophobic
        self.hydrophilic = hydrophilic
        self.height      = height
        self.L           = L
        self.R           = R
        self.D           = D

        self.d           = D/R
        self.l           = L/R

        # reduced volume of liquid
        self.V           = np.pi*self.l*self.l*height/R
        # definition of the interface
        # In the equilibrium state, half of the particle is merged under the water
        # we can computed the location of the interface.
        self.interface   = 0.5 *4.0/3.0*np.pi/self.l/self.l/np.pi*self.R + self.height
        
        # Pulling up from a finite filling angle with a constant hydrophobic contact angle
        # to filling angle 90 degree with a constant hydrophobic contact angle.

        self.psis = np.linspace(179, 90.0, 100)
        self.Vs   = []

        for psi in self.psis:
            model = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = self.hydrophobic, psi = psi)
            self.Vs.append(model.V)

        self.Vs = np.array(self.Vs)

        #plt.plot(self.psis, self.Vs)
        #plt.show()
        
        indexmin = -1
        indexmax = np.where(self.Vs == self.Vs.max())[0][0]

        psi_min = self.psis[indexmin]
        psi_max = self.psis[indexmax]

        psis = np.linspace(psi_max, psi_min, 100)
        Ds   = []

        for psi in psis:

            model  = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = self.hydrophobic, psi = psi)
            deltah = (self.V - model.V)/np.pi/self.l/self.l
            Ds.append(self.D + deltah*self.R)

        Ds = np.array(Ds)
        self.models = []

        for i in range(len(psis)):

            model = yl.YL(R = self.R, L = self.L, D = Ds[i], theta1 = self.hydrophobic, psi = psis[i])
            self.models.append(model)

        # End of sliding contact angle at the hydropphobic part
        # Transition to a pinned contact angle condtion at the equator of the Janus
        # particle.
        # Filling angle is maintained at a right angle, while contact angle at the equator
        # decrease from hydrophobic angle to hydrophilic angle.

        theta1s = np.linspace(self.hydrophobic, self.hydrophilic, 100)
        Vs = []

        for theta1 in theta1s:
            model = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = theta1, psi = 90.0)
            Vs.append(model.V)

        Vs = np.array(Vs)

        #plt.plot(theta1s, Vs)
        #plt.show()
        
        indexmin = np.where(Vs == Vs.max())[0][0] # should be 0
        indexmax = np.where(Vs == Vs.min())[0][0] # should be -1, (99). 

        theta1_max = theta1s[indexmax]
        theta1_min = theta1s[indexmin]

        theat1s = np.linspace(theta1_min, theta1_max, 100)
        Ds = []

        for theta1 in theta1s:
            model  = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = theta1, psi = 90.0)
            deltah = (self.V - model.V)/np.pi/self.l/self.l
            Ds.append(self.D + deltah*self.R)

        Ds = np.array(Ds)

        for i in range(len(theta1s)):

            model = yl.YL(R = self.R, L = self.L, D = Ds[i], theta1 = theta1s[i], psi = 90.0)
            self.models.append(model)

        # Keep pulling the particles to sliding contact angle condition at the
        # hydrophilic part of the Janus particle till the rupture of the meniscus

        self.psis = np.linspace(90.0, 1.0, 100)
        self.Vs   = []

        for psi in self.psis:

            model = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = self.hydrophilic, psi = psi)
            self.Vs.append(model.V)

        self.Vs = np.array(self.Vs)

        #plt.plot(self.psis, self.Vs)
        #plt.show()
        
        indexmin = 0
        indexmax = np.where(self.Vs == self.Vs.min())[0][0]

        psi_min = self.psis[indexmax]
        psi_max = self.psis[indexmin]

        psis = np.linspace(psi_max, psi_min, 100)
        Ds   = []

        for psi in psis:

            model  = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = self.hydrophilic, psi = psi)
            deltah = (self.V - model.V)/np.pi/self.l/self.l
            Ds.append(self.D + deltah*self.R)

        Ds = np.array(Ds)

        for i in range(len(psis)):

            model = yl.YL(R = self.R, L = self.L, D = Ds[i], theta1 = self.hydrophilic, psi = psis[i])
            self.models.append(model)


if __name__ == "__main__":

    # model = Janus( hydrophobic = 123.25, hydrophilic = 52.68, height = 50.9, L = 49.3, R = 10.9, D = 50.0)

    # model = Janus( hydrophobic = 112.06, hydrophilic = 51.832688, height = 50.9, L = 49.3, R = 10.5, D = 50.0)
    model = Janus( hydrophobic = 112.06, hydrophilic = 112.06, height = 50.9, L = 49.3, R = 10.5, D = 50.0)

    displacement, force, L, R, D, theta1, psi, regime = [], [], [], [], [], [], [], []
        
    for i in range(len(model.models)):
        displacement.append((model.models[i].D + model.R  - model.interface)/model.R)
        force.append(model.models[i].force)
        L.append(model.models[i].L)
        R.append(model.models[i].R)
        D.append(model.models[i].D)
        theta1.append(model.models[i].theta1)
        psi.append(model.models[i].psi)
        if model.models[i].psi > np.pi/2.0 - 0.01 and model.models[i].psi < np.pi/2.0 + 0.01:
            regime.append("pinned")
        else:
            regime.append("sliding")
                            
    plt.plot(displacement, force, '-')
    plt.show()
    
    ofile = open("sample.txt", 'w')
    ofile.write("relative_distance  force  psi\n")
    for i in range(len(force)):

        ofile.write("{0} {1} {2}\n".format(displacement[i], force[i], psi[i]))

    ofile.close()

    ofile = open("info.txt", 'w')
    ofile.write("R L D theta1 psi regime\n")
    for i in range(len(force)):

        ofile.write("{0} {1} {2} {3} {4}\n".format(R[i], L[i], D[i], theta1[i]/np.pi*180.0, psi[i]/np.pi*180.0))

    ofile.close()