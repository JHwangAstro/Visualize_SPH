import numpy as np
import math


class Planet:

    def __init__(self, particles, star):

        G = 1.0
        ms = 1.98855*math.pow(10.,33.)                          #Solar mass in g
        me = 5.972*math.pow(10.,27.)                            #Earth mass in g

        self.m = np.sum([particle.m for particle in particles])
        self.n = len(particles)
        self.r = self.calculate_center_of_mass(particles)
        self.v = self.calculate_velocity(particles)
        self.mu = G*(star.m+self.m)
        self.aei = self.calculate_orbital_elements(self.r, self.v)
        self.E = -self.mu/(2.*self.aei['a'])
        self.h = math.sqrt(self.aei['a']*(1.-self.aei['e']**2.)*self.mu)
        self.particles = [particle.i for particle in particles]
        self.specific_e = -0.5*(star.m*self.m)/self.aei['a']
        self.specific_h = self.m*star.m/(self.m+star.m)*math.sqrt(self.aei['a']*(1.-self.aei['e']**2.)*(star.m+self.m))

        # Convert mass to Solar masses
        self.m = self.m*ms/me

        print self.r, self.v, self.aei

#        self.profile = PlanetProfile()



    def calculate_velocity(self, particles):
        if self.m == 0: return -1
        vx = np.sum([particle.m*particle.vx for particle in particles])/self.m
        vy = np.sum([particle.m*particle.vy for particle in particles])/self.m
        vz = np.sum([particle.m*particle.vz for particle in particles])/self.m
        return np.array([vx,vy,vz])



    def calculate_center_of_mass(self, particles):
        if self.m == 0: return -1
        x = np.sum([particle.m*particle.x for particle in particles])/self.m
        y = np.sum([particle.m*particle.y for particle in particles])/self.m
        z = np.sum([particle.m*particle.z for particle in particles])/self.m
        return np.array([x,y,z])



    #in:  r,v
    #out: a,e,i,g,n,M
    def calculate_orbital_elements(self, r, v):
        if self.m == 0: return -1

        eps = math.pow(10.,-6)

        h=np.cross(r,v)
        norm=np.cross([0,0,1],h)

        evec = ((np.linalg.norm(v)**2.-self.mu/np.linalg.norm(r))*r-np.dot(r,v)*v)/self.mu
        e = np.linalg.norm(evec)

        energy = np.linalg.norm(v)**2./2.-self.mu/np.linalg.norm(r)

        if math.fabs(e-1.0)>eps:
            a = -self.mu/(2.*energy)
            p = a*(1.-e**2.)
        else:
            a = math.inf
            p = np.linalg.norm(h)**2./self.mu

        i = math.acos(h[2]/np.linalg.norm(h))

        n = math.acos(norm[0]/np.linalg.norm(norm))
        if norm[1]<0:
            n = 2.*math.pi-n

        g = math.acos(np.dot(norm,evec)/(np.linalg.norm(norm)*e))
        if evec[2]<0.0:
            g = 2.*math.pi-g

        M = math.acos(np.dot(evec,r)/(e*np.linalg.norm(r)))
        if np.dot(r,v)<0.0:
            M = 2.*math.pi-M

        # Convert semi-major axis to AU
        a = a*0.00465

        return {'a': a, 'e': e, 'i': i, 'g': g, 'n': n, 'M': M}



    def add_core(self, core):

        self.core = core
        self.core_mass_fraction = core.m/self.m
        self.gas_mass_fraction = 1.-self.core_mass_fraction







class PlanetProf:

  def __init__(self,Planet):
    nbins = 100
    self.r,self.rho,self.m,self.T,self.P,self.u,self.gamm,self.cc = Calcs.Generate_Profile(nbins,Planet.com,Planet.Particles)




























