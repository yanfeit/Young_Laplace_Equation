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
The script generates a quasi-static pulling process of a homogenous particle from 
the liquid.
"""

class Homogeneous(object):
    """   
    input: 
        angle  :                   contact angle of the two part of the  particle
        height :                   liquid height without considering the particle
        L :                        container radius      
        R :                        radius of the particle
        D :                        distance between the bottom of the particle and bottom of the containder.
        interface :                location of the interface if the particle is in equilibrium.

    methods :
        pullproc : pulling process of the particle and store to HDF5 files.
        frmproc  : given the volume of liquids, solve the menisci and force, and store to HDF5 files.
    """
    def __init__(self, angle = 30.0, height = 50.0, L = 50.0, R = 10.0, D = 50.0):

        self.angle       = angle
        self.height      = height
        self.L           = L
        self.R           = R
        self.D           = D

        self.d           = D/R
        self.l           = L/R

        # reduced volume of liquid
        self.V           = np.pi*self.l*self.l*height/R
        # definition of the interface
        # In the equilibrium state, particle is merged under the water

        theta_angle = self.angle / 180.0 * np.pi
        costheta    = np.cos(theta_angle)
        cos3theta   = costheta * costheta * costheta
        self.interface = self.height + 1/3.0 * (2 + 3*costheta - cos3theta)/self.l/self.l * self.R
    

    def pullproc(self):

        # Step (1) : 
        # We need a relation between the filling angle v.s volume
        # give a contant particle location
        self.psis = np.linspace(179.0, 1.0, 100)
        self.Vs   = []

        for psi in self.psis:
            model = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = self.angle, psi = psi)
            self.Vs.append(model.V)

        self.Vs = np.array(self.Vs)
        
        indexmin = np.where(self.Vs == self.Vs.min())[0][0]
        indexmax = np.where(self.Vs == self.Vs.max())[0][0]
        
        # Notice that the line of filling angle v.s volume is not 
        # monotonic, here needs some stable analysis. We think that
        # large filling anlge is more stable than small filing angle. 
        # (An stability analysis is required, but I don't have the time)
        # And we only consider filling in this region
        # psi in [psi[indexmin], psi[indexmax]]
        psi_min = self.psis[indexmin]
        psi_max = self.psis[indexmax] 

        # Step 2: Given the constrained volume of liquids, solve D paratmers as an input.
        psis = np.linspace(psi_max, psi_min, 100)
        Ds   = []

        for psi in psis:

            model  = yl.YL(R = self.R, L = self.L, D = self.D, theta1 = self.angle, psi = psi)
            deltah = (self.V - model.V)/np.pi/self.l/self.l
            # compensate the D
            Ds.append(self.D + deltah*self.R)

        Ds = np.array(Ds)
        self.models = []

        # Step (1.3) Recompute the model by YL equation.
        for i in range(len(psis)):

            model = yl.YL(R = self.R, L = self.L, D = Ds[i], theta1 = self.angle, psi = psis[i])
            self.models.append(model)


    def storeHDF5(self, ofilename="theory_homo_data.h5"):
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
        targetpsi = ((self.models[idx].psi + self.models[idx-1].psi)/2.0  / np.pi * 180.0)

        self.targetmodel = yl.YL(R = self.R, L = self.L, D = targetD, theta1 = self.angle, psi = targetpsi)


        

if __name__ == "__main__":

    model = Homogeneous( angle = 51.83, height = 50.9, L = 49.3, R = 10.5, D = 40.0)
    model.frmproc()
    model.storeHDF5()
