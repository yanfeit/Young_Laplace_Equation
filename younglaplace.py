#    -------------------------------------------------------------------------- 
#    Zhejiang Lab, Zhejiang, China
#    Yanfei Tang, tangyf@zhejianglab.com

#    Copyright (2023) Zhejiang Lab.  Under the terms of Contract
#    2021PB0AC02 with Zhejiang Lab, the Zhejiang Government retains
#    certain rights in this software.  This software is distributed under
#    the GNU General Public License.
#    -------------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.patches import Wedge
import scipy.special as spe
import sys

docstr = """
Solve Young-Laplace Equation of nanoparticle in a confining cylinder wall.

Reference: 
1. Orr, Scriven and Rivas; Pendular rings between solids: meniscus properties and capillary force, JFM 67 723 (1975)
2. Rubinstein and Fel; Theory of Pendular Rings Revisited, ArXiv 1207.7096v1 (2012)
"""

class YL(object):
    """
    R      : radius of nanoparticle
    L      : radius of cylinder
    D      : distance between of bottom of nanoparticle and coordinate
    d      : scaled D, D/R
    l      : scaled radius of cylinder, L/R

    theta1 : contact angle at the sphere
    psi    : filling angle
    t      : angle of the normal to meniscus with the vertical axis
    t1     : angle of the normal to meniscus at sphere with the vertical axis

    H      : dimensionless mean curvature
    c      : parameter in the formula

    k2     : k square, used in elliptic integral!
    x      : scaled horizontal coordination
    y      : scaled vertical coordination

    y1     : scaled vertical boundary condition at sphere
    force  : scaled force on nanoparticle / (2*pi*gamma*R)
    h      : vertical distance between center of the nanoparicle 
            and boundary condition at the wall!
    Dy     : vertical scaled interface distortion length
    """

    def __init__(self, R = 10.35, L = 50.0, D = 50.0, theta1 = 30.0, psi = 130.0):
        
        self.R = R
        self.D = D
        self.L = L
        self.d = self.D/self.R
        self.l = self.L/self.R

        self.theta1 = theta1/180.0*np.pi
        self.psi    = psi/180.0*np.pi

        # meriodinal slope, need explanation here!
        # meridional   经线的
        # azimuthal    方位角的

        if self.theta1 + self.psi < np.pi:
            # For rising meniscus
            self.t1 = self.theta1 + self.psi
            self.t  = np.linspace(self.t1, np.pi, 1000)
        elif self.theta1 + self.psi > np.pi:
            # For depressing meniscus
            self.t1 = self.theta1 + self.psi - np.pi
            self.t  = np.linspace(self.t1, 0, 1000) 
        else:
            sys.exit("t1 = np.pi, Flat surface.\n")

        self.H     = np.sin(self.psi)*np.sin(self.t1)/(np.sin(self.psi)*np.sin(self.psi) - self.l*self.l)
        self.c     = 4 * self.H * np.sin(self.psi) * (self.H * np.sin(self.psi) - np.sin(self.t1))

        self.k2    = -1.0/self.c
        self.s     = -1 

        self.x     = 1/2.0/self.H * (np.sin(self.t) \
                                     + self.s * np.sqrt(np.sin(self.t) * np.sin(self.t) + self.c))

        self.y1    = 1 + self.d - np.cos(self.psi)
        self.part1 = 1.0/2.0/self.H * (np.cos(self.t1) - np.cos(self.t))
        self.part2 = self.s * np.sqrt(self.c)/2.0/self.H * \
                     (spe.ellipeinc(self.t,self.k2) - spe.ellipeinc(self.t1,self.k2) - spe.ellipkinc(self.t,self.k2) + spe.ellipkinc(self.t1, self.k2))

        self.y     = self.y1 + self.part1 + self.part2
        self.force = np.sin(self.psi) * np.sin(self.theta1 + self.psi) - self.H * np.sin(self.psi) * np.sin(self.psi)
        self.h     = 1 + self.d - self.y[-1]
        self.Dy    = self.y[0] - self.y[-1]

        self.t2    = self.t[-1]
        self.y2    = self.y[-1]

        cost1      = np.cos(self.t1)
        cost2      = np.cos(self.t2)
        cos3t1     = cost1*cost1*cost1
        cos3t2     = cost2*cost2*cost2
        self.Jpart1 = (4.0 + self.c) * (-cost1 + cost2)
        self.Jpart2 = 4.0/3.0 * (cos3t1 - cos3t2)
        self.Jpart3 = -np.sqrt(self.c) * ((8 + self.c)/3.0 * (spe.ellipeinc(self.t1, self.k2) - spe.ellipeinc(self.t2, self.k2)) -\
                                          (4 + self.c)/3.0 * (spe.ellipkinc(self.t1, self.k2) - spe.ellipkinc(self.t2, self.k2)))
        self.Jpart4 = 2.0/3.0 * (np.sin(2*self.t1) * np.sqrt(np.sin(self.t1)*np.sin(self.t1) + self.c) - \
                                 np.sin(2*self.t2) * np.sqrt(np.sin(self.t2)*np.sin(self.t2) + self.c))
        self.Js     = self.Jpart1 + self.Jpart2 + self.Jpart3 + self.Jpart4
        self.Vc     = np.pi/8.0/self.H/self.H/self.H*self.Js
        cos3psi     = np.cos(self.psi) * np.cos(self.psi) * np.cos(self.psi) 
        self.V      = self.Vc + np.pi * self.l * self.l * self.y2 - 1/3.0 * np.pi * (2 - 3*np.cos(self.psi) + cos3psi)


    def __str__(self):
        """
        printalbe info        
        print "h   force    H   psi\n"
        """
        self.entry = "{0:5f} \t {1:5f} \t {2:5f} \t {3:5f}\n".format(self.h, self.force, self.H, self.psi)
        return self.entry

    def draw(self):
        """
        Draw pendular ring and nanoparticle.
        """
        fig = plt.figure(0)
        ax = fig.add_subplot(111, aspect = 'equal', ylim = (1.0, 9.0), xlim = (-5, 5))
        e1 = Wedge((0, self.d + 1), 1, theta1 = 0.0,   theta2 = 180.0, fill = True, fc = "g", ec = "g", lw = 0.0)
        e2 = Wedge((0, self.d + 1), 1, theta1 = 180.0, theta2 = 360.0, fill = True, fc = "r", ec = "r", lw = 0.0)
        ax.add_patch(e1)
        ax.add_patch(e2)
        ax.plot( self.x, self.y, 'C0', linewidth = 4.0)
        ax.plot(-self.x, self.y, 'C0', linewidth = 4.0)
        ax.set_xlabel(r"$r/R$")
        ax.set_ylabel(r"$z/R$")
        ax.set_title(r"$\psi = {0:.1f}, \theta = {1:.2f},d = {2:.2f},V = {3:.1f}$".format(self.psi/np.pi*180.0, self.theta1/np.pi*180.0, self.d, self.V))
        
        plt.tight_layout()
        plt.savefig("theta1{0:.1f}psi{1:.1f}.png".format(self.theta1/np.pi*180.0, self.psi/np.pi*180.0), transparent = True)
        plt.clf()
        
        
    def meniscus(self):
        """
        local curvature detail, azimuthal 方位角曲率 meridional 经角的曲率
        """
        info = "# Local information of Meniscus\n" +\
               "t        x        y        meridional azimuthal \n"
        for i in range(len(self.t)):
            azimuthal = np.sin(self.t[i])/self.x[i]
            meridional = 2*self.H - azimuthal
            info += "{0:3f} {1:3f} {2:3f} {3:3f} {4:3f}\n".format(self.t[i], self.x[i], self.y[i], meridional, azimuthal)

        return info


if __name__ == "__main__":
    # L:      length of your cylinder
    # theta1: contact angle at sphere
    # psi   : filling angle

    model = YL(R = 10.6, L = 49.24, D = 50.0, theta1 = 180.0, psi = 90.0)
    print(model)
    model.draw()
    print(model.meniscus())
        