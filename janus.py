#    -------------------------------------------------------------------------- 
#    Zhejiang Lab, Zhejiang, China
#    Yanfei Tang, tangyf@zhejianglab.com

#    Copyright (2023) Zhejiang Lab. This work is supported under the terms of 
#    contract 2021PB0AC02 with Zhejiang Lab. 
#    This software is distributed under the GNU General Public License.
#    -------------------------------------------------------------------------- 

import sys
import younglaplace as yl
import numpy as np
import matplotlib.pyplot as plt
import h5py

docstr = """
The script generates a quasi-static pulling process of a janus particle from 
the liquid. The Janus particle has a hydrophobic part of the hemisphere at the top and 
a hydrophilic part of the hemisphere at the bottom. 

This is a physical decription of the dynamic process:
For a Janus particle with a hydropholic part large than 90 degree and a hydrophilic part
smaller than 90 degree. The pulling process of the particle is composed by three steps:

(1) Sliding regime: The filling angle decreases from phi_1 to 90 degree while the particle  maintains 
a constant contact angle at the upper part of the surface of the janus particle. 
phi_1 is the critical angle that forms the menisci.
(2) Pinned regime: When the particle is pulled in this regime, the filling angle maintains a constant.
In the meantime, the menisci pinned at the equator of the sphere. The contact angle decreases from the 
hydrophobic degree to hydrophilic degree.
(3) Sliding regime: the pinned menisci leaves the equator of the sphere. The filling angle decreases
from 90 degree to phi_2 while the particles maintains a constant contact angle at the lower part of 
the surface of the janus particle. phi_2 is the critical angle that the menisci breaks.

In the beginning of sliding region, the Janus particle is merged under the liquid with menisci formed. 
The filling angle decreases with the undergoing pulling process. 
It then pins at the equator of the sphere with a right filling angle. At some point, the process 
is back to the sliding region with a rising menisci.
"""

class Janus(object):
    """   
    input: 
        hydrophobic, hydrophilic : contact angle of the two part of the janus particle
        height :                   liquid height without considering the particle
        L :                        container radius      
        R :                        radius of the janus particle
        D :                        distance between the bottom of the particle and bottom of the containder.
        interface :                location of the interface if the particle is in equilibrium.
        regime :                   Sliding(1), pinned(2)

    methods :
        pullproc : pulling process of the particle and store to HDF5 files.
        frmproc  : given the volume of liquids, solve the menisci and force, and store to HDF5 files.
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
        # if self.hydrophilic != self.hydrophobic:
            # this is the case for janus particle,
            # we consider hydrophobic > 90.0 and hydrophilic < 90.0
            # other cases will be considered in future.
        self.interface = 0.5 *4.0/3.0*np.pi/self.l/self.l/np.pi*self.R + self.height
        # else:
            # the case should seamlessly reduces to homogenous spherical particle case.
        #    theta_angle = self.hydrophilic / 180.0 * np.pi
        #    costheta    = np.cos(theta_angle)
        #    cos3theta   = costheta*costheta*costheta
        #    self.interface = 1/3.0 * (2 + 3*costheta + cos3theta)/self.l/self.l + self.height
    
    def pullproc(self):

        self.regime = []
        # Step (1) : 
        # Pulling up from a finite filling angle with a constant hydrophobic contact angle
        # to filling angle 90 degree with a constant hydrophobic contact angle

        # Step (1.1) Solve the critical phi_1, approximately.
        self.psis = np.linspace(179, 90.0, 100)
        self.Vs   = []

        for psi in self.psis:
            model = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = self.hydrophobic, psi = psi)
            self.Vs.append(model.V)

        self.Vs = np.array(self.Vs)
        
        indexmin = -1
        indexmax = np.where(self.Vs == self.Vs.max())[0][0]
        
        # psi_min = 90, psi_max is the critical filling angle
        # that we consider the menisci forms.
        psi_min = self.psis[indexmin]
        psi_max = self.psis[indexmax] 

        # Step (1.2) Given the constrained volume of liquids, solve D paratmers as an input.
        psis = np.linspace(psi_max, psi_min, 100)
        Ds   = []

        for psi in psis:

            model  = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = self.hydrophobic, psi = psi)
            deltah = (self.V - model.V)/np.pi/self.l/self.l
            # compensate the D
            Ds.append(self.D + deltah*self.R)

        Ds = np.array(Ds)
        self.models = []

        # Step (1.3) Recompute the model by YL equation.
        for i in range(len(psis)):

            model = yl.YL(R = self.R, L = self.L, D = Ds[i], theta1 = self.hydrophobic, psi = psis[i])
            self.models.append(model)
            self.regime.append(1)

        # End of sliding contact angle at the hydropphobic part
        # Transition to a pinned contact angle condtion at the equator of the Janus
        # particle.

        # Step (2)
        # Filling angle is maintained at a right angle, while contact angle at the equator
        # decrease from hydrophobic angle to hydrophilic angle.

        # Step (2.1)
        theta1s = np.linspace(self.hydrophobic, self.hydrophilic, 100)
        # Vs = []

        # for theta1 in theta1s:
        #     model = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = theta1, psi = 90.0)
        #     Vs.append(model.V)

        # Vs = np.array(Vs)
        
        # indexmin = np.where(Vs == Vs.max())[0][0] # should be 0
        # indexmax = np.where(Vs == Vs.min())[0][0] # should be -1, (99). 

        # theta1_max = theta1s[indexmax]
        # theta1_min = theta1s[indexmin]

        # theta1s = np.linspace(theta1_min, theta1_max, 100)
        Ds = []

        for theta1 in theta1s:
            model  = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = theta1, psi = 90.0)
            deltah = (self.V - model.V)/np.pi/self.l/self.l
            Ds.append(self.D + deltah*self.R)

        # Step (2.2)
        Ds = np.array(Ds)

        for i in range(len(theta1s)):

            model = yl.YL(R = self.R, L = self.L, D = Ds[i], theta1 = theta1s[i], psi = 90.0)
            self.models.append(model)
            self.regime.append(2)

        # Step (3)
        # Keep pulling the particles to sliding contact angle condition at the
        # hydrophilic part of the Janus particle till the rupture of the meniscus
        
        # Step (3.1) : Solve the critical filling angle
        self.psis = np.linspace(90.0, 1.0, 100)
        self.Vs   = []

        for psi in self.psis:

            model = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = self.hydrophilic, psi = psi)
            self.Vs.append(model.V)

        self.Vs = np.array(self.Vs)
        
        indexmin = 0
        indexmax = np.where(self.Vs == self.Vs.min())[0][0]

        psi_min = self.psis[indexmax]
        psi_max = self.psis[indexmin]  # indexmin = 0, pis_max=90

        # Step (3.2) : Given the constrained volume of liquids, solve D paratmers as an input.
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
            self.regime.append(1)


    def storeHDF5(self, ofilename="theory_data.h5"):
        # Store Data to HDF5 format
        with h5py.File(ofilename, 'w') as f:
            
            for i in range(len(self.models)):

                group = f.create_group("{}".format(i))
                group.create_dataset("R", data = self.models[i].R)
                group.create_dataset("D", data = self.models[i].D)
                group.create_dataset("L", data = self.models[i].L)
                group.create_dataset("d", data = self.models[i].d)
                group.create_dataset("l", data = self.models[i].l)
                group.create_dataset("theta1", data = self.models[i].theta1)
                group.create_dataset("psi", data = self.models[i].psi)
                group.create_dataset("x", data = self.models[i].x)
                group.create_dataset("y", data = self.models[i].y)
                group.create_dataset("force", data = self.models[i].force)
                group.create_dataset("V", data = self.models[i].V)
                group.create_dataset("deltaz", data = self.models[i].deltaz)
                group.create_dataset("h0", data=self.models[i].h0)
                group.create_dataset("regime", data = self.regime[i])
                group.create_dataset("deltazp", data = self.models[i].deltazp)
                # group.create_dataset("d", data = self.models[i].d)

            groupt = f.create_group("target")
            groupt.create_dataset("R", data = self.targetmodel.R)
            groupt.create_dataset("D", data = self.targetmodel.D)
            groupt.create_dataset("L", data = self.targetmodel.L)
            groupt.create_dataset("d", data = self.targetmodel.d)
            groupt.create_dataset("l", data = self.targetmodel.l)
            groupt.create_dataset("theta1", data = self.targetmodel.theta1)
            groupt.create_dataset("psi", data = self.targetmodel.psi)
            groupt.create_dataset("x", data = self.targetmodel.x)
            groupt.create_dataset("y", data = self.targetmodel.y)
            groupt.create_dataset("force", data = self.targetmodel.force)
            groupt.create_dataset("V", data = self.targetmodel.V)
            groupt.create_dataset("deltaz", data = self.targetmodel.deltaz)
            groupt.create_dataset("h0", data=self.targetmodel.h0)
            groupt.create_dataset("deltazp", data=self.targetmodel.deltazp)



    def frmproc(self):
        """
        To solve the profile of the menisci given a reasonable D,
        we should first find out a whole pulling process.
        Then we match D by using the data from self.models.d.
        """
        self.pullproc()

        # Determine D is in the range of the pulling process which the menisci can form.
        if self.d < self.models[0].d or self.d > self.models[-1].d:
            sys.exit("Can't form a menisci")
        
        idx = 0
        for i in range(len(self.models)):
            if self.d < self.models[i].d:
                idx = i
                break

        targetD = ((self.models[idx].d + self.models[idx-1].d)/2.0  * self.R ) # convert to nonscaled unit
        targettheta1 = ((self.models[idx].theta1 + self.models[idx-1].theta1)/2.0  /np.pi * 180.0 )
        targetpsi = ((self.models[idx].psi + self.models[idx-1].psi)/2.0  / np.pi * 180.0)

        self.targetmodel = yl.YL(R = self.R, L = self.L, D = targetD, theta1 = targettheta1, psi = targetpsi)


        

if __name__ == "__main__":

    model = Janus( hydrophobic = 112.06, hydrophilic = 51.83, height = 50.9, L = 49.3, R = 10.5, D = 50.0)
    model.frmproc()
    model.storeHDF5()

