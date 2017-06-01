
import math
import numpy as np
import eos

class Particle:

    def __init__(self, i, data, params):

        self.i = i

        (dum, self.x, self.y, self.z, self.m, self.h, self.rho, \
        self.vx, self.vy, self.vz, self.vxdot, self.vydot, self.vzdot, \
        self.u, self.udot, self.grpot, self.mmw, self.cc, self.divv, dum) = data

        self.r = np.array([self.x, self.y, self.z])
        self.v = np.array([self.vx, self.vy, self.vz])

        if params:
            #self.T = eos.get_temperature(self.rho, self.u, self.mmw, params.K, params.gam, params.beta, params.xmix, self.cc, params.ul0, params.u0)
            #self.p = eos.get_pressure(self.rho, self.u, self.mmw, params.K, params.gam, params.beta, params.xmix, self.cc, params.ul0, params.uh0, params.u0)
            #self.t_cool = eos.get_cooling_time(self.m, self.u, self.rho, self.T)
            # copy over all the attributes from params
            for k, v in params.__dict__.items():
                self.__dict__[k] = v
        else:
            #self.T = None
            #self.t_cool = None
            self.core = 0
            self.core2 = 0





    def get_p(self):

        self.p = eos.get_pressure(self.rho, self.u, self.mmw, self.K, self.gam, self.beta, self.xmix, self.cc, self.ul0, self.uh0, self.u0)




    def get_T(self):

        self.T = eos.get_temperature(self.rho, self.u, self.mmw, self.K, self.gam, self.beta, self.xmix, self.cc, self.ul0, self.uh0, self.u0)




    def calculate_angular_momentum(self, v, r, j):

        g   = 6.67390*math.pow(10,-8)                               # gravitational constant [cgs]
        rs  = 6.9599*math.pow(10,10)                                # unit of length [cm]
        ms  = 1.9891*math.pow(10,33)                                # unit of mass [g]
        ts  = math.sqrt(math.pow(rs,3)/(g*ms))                      # unit of time [s]

        # calculate the angular momentum by projecting onto the spin axis
        d = np.cross(list(self.r - (r + j)), list(self.r - r))/np.linalg.norm(j)

        self.j = self.m*np.linalg.norm(np.cross(d, list(self.v - v)))*ms*rs*rs/ts
        self.spec_j = self.j/self.m




class ParticleParameters:

    def __init__(self, i, data, u0, core, core2, K, gam):

        gravconst = 6.67390*10**-8.
        rs  = 6.9599*math.pow(10,10)                                # unit of length [cm]
        ms  = 1.9891*math.pow(10,33)                                # unit of mass [g]

        self.i = i

        (dum, self.rho0, self.xmix, self.beta, self.T0, self.umix, \
        self.uh0, self.ul0, self.dmix, self.rhoh, self.rhol, dum) = data

        self.u0 = u0*gravconst*ms/rs        # initial internal energy
        self.core = core    # flag for mantle/iron particles
        self.core2 = core2  # flag for iron particles
        self.K = K
        self.gam = gam




















