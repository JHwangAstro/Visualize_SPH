import numpy as np

class Particle:
    def __init__(self,i,f):

        #index
        self.i = i

        dum        = np.fromfile(f,dtype=np.int32,count=1)           #dummy variable
        self.x     = np.fromfile(f,dtype=float,count=1)              #x position
        self.y     = np.fromfile(f,dtype=float,count=1)              #y position
        self.z     = np.fromfile(f,dtype=float,count=1)              #z position
        self.m     = np.fromfile(f,dtype=float,count=1)              #mass
        self.h     = np.fromfile(f,dtype=float,count=1)              #smoothing length
        self.rho   = np.fromfile(f,dtype=float,count=1)              #density
        self.vx    = np.fromfile(f,dtype=float,count=1)              #x velocity
        self.vy    = np.fromfile(f,dtype=float,count=1)              #y velocity
        self.vz    = np.fromfile(f,dtype=float,count=1)              #z velocity
        self.vxdot = np.fromfile(f,dtype=float,count=1)              #x acceleration
        self.vydot = np.fromfile(f,dtype=float,count=1)              #y acceleration
        self.vzdot = np.fromfile(f,dtype=float,count=1)              #z acceleration
        self.u     = np.fromfile(f,dtype=float,count=1)              #internal energy
        self.udot  = np.fromfile(f,dtype=float,count=1)              #change in internal energy
        self.grpot = np.fromfile(f,dtype=float,count=1)              #gravitational potential energy
        self.mmw   = np.fromfile(f,dtype=float,count=1)              #mean molecular weight
        self.cc    = np.fromfile(f,dtype=np.int32,count=1)           #composition flag
        self.divv  = np.fromfile(f,dtype=float,count=1)              #mean molecular weight
        dum        = np.fromfile(f,dtype=np.int32,count=1)           #dummy variable



    def determine_bound(self, planets, star):

        # constants
        G = 1.0
        RS = 1.0

        # calculate the distances, velocities and energies from each particle to the respective cores
        dx = np.array([((self.x -planet.r[0])**2 + (self.y -planet.r[1])**2 + (self.z -planet.r[2])**2)**0.5 for planet in planets])
        dv = np.array([((self.vx-planet.v[0])**2 + (self.vy-planet.v[1])**2 + (self.vz-planet.v[2])**2)**0.5 for planet in planets])
        E  = np.array([self.m*(dv[i]**2./2.0-G*planets[i].m/dx[i]) for i in range(len(planets))])

        # calculate the distances, velocities and energies of each particle
        r = (self.x**2 + self.y**2 + self.z**2)**0.5
        dv2 = (self.vx**2 + self.vy**2 + self.vz**2)**0.5
        Ecom = self.m*(dv2**2/2.0 + self.grpot)

        # assign mass bound numbers
        # 0 - Bound to planet 1
        # 1 - Bound to planet 2
        # 2 - Bound to the system
        # 3 - Accreted
        # 4 - Ejecta
        self.bound = 4
        if Ecom < 0.0: self.bound = 2
        if len(planets) == 2:
            if (E[1] < 0.0) & (E[1] < E[0]): self.bound = 1
        if (E[0] < 0.0) & (E[0] < E[1]): self.bound = 0
        if r < RS: self.bound = 3




















