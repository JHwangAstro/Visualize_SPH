import numpy as np
import math
import matplotlib.pyplot as plt

import mercury_routines
import orbit_calc


class Ref_system:

    def __init__(self, name):

        G = 1.0

        self.m, self.d, self.a, self.e, self.star_mass = self.get_system_attributes(name)
        self.exists = 1
        self.pratio = 1
        self.rho_ratio = 1

        if self.m[0]:
            self.mu = G*(self.star_mass+self.m)
            self.E = -self.mu/(2.*self.a)
            self.h = (self.a*(1.-self.e**2.)*self.mu)**0.5
            self.specific_e = -0.5*(self.star_mass*self.m)/self.a
            self.specific_h = self.m*self.star_mass/(self.m+self.star_mass)*(self.a*(1.-self.e**2.)*(self.star_mass+self.m))**.5
            self.rho_ratio = self.d[0]/self.d[1]
            self.pratio = (self.a[1]/self.a[0])**1.5*((self.star_mass+self.m[1])/(self.star_mass+self.m[0]))**-.5
            self.AMD = orbit_calc.calculate_AMD(self.star_mass, self.m, [self.a[0], self.e[0]], [self.a[1], self.e[1]])
            if self.pratio < 1.: self.pratio = 1./self.pratio
        else:
            self.exists = 0




    # returns attributes of the system
    # m - array of planet masses in MS
    # a - array of semi-major axes in AU
    # e - array of eccentricities
    # star_mass - mass of host star in MS
    def get_system_attributes(self, name):

        ms = 1.98855*math.pow(10.,33.)                          #Solar mass in g
        me = 5.972*math.pow(10.,27.)                            #Earth mass in g

        if name == 'Kepler-36':
            m = [4.45, 8.08]
            a = [0.1154, 0.1283]
            e = [0.04, 0.04]
            d = [0.86, 6.8]
            star_mass = 1.071

        else:
            print 'no system with this name'
            m = [None]
            a = [None]
            e = [None]
            d = [None]
            star_mass = None

        if m[0]: return np.array(m)*me/ms, np.array(d), np.array(a), np.array(e), star_mass
        if not m[0]: return m, d, a, e, star_mass


























