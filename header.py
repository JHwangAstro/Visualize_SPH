import numpy as np
import math

class Header:
    def __init__(self,f):

        g   = 6.67390*math.pow(10,-8)                                 #gravitational constant [cgs]
        rs  = 6.9599*math.pow(10,10)                                  #unit of length [cm]
        ms  = 1.9891*math.pow(10,33)                                  #unit of mass [g]
        ts  = math.sqrt(math.pow(rs,3)/(g*ms))/31556926.              #unit of time [yr]

        dum         = np.fromfile(f,dtype=np.int32,count=1)           #dummy variable
        self.ntot   = np.fromfile(f,dtype=np.int32,count=1)           #total number of particles
        self.nnopt  = np.fromfile(f,dtype=np.int32,count=1)           #optimal number of neighbors
        self.hmin   = np.fromfile(f,dtype=float,count=1)              #smoothing length for black holes
        self.hmax   = np.fromfile(f,dtype=float,count=1)              #smoothing length for black holes
        self.sep0   = np.fromfile(f,dtype=float,count=1)              #initial separation of two stars
        self.tf     = np.fromfile(f,dtype=float,count=1)              #time the code stops in code units
        self.dtout  = np.fromfile(f,dtype=float,count=1)              #time interval between writing
        self.nout   = np.fromfile(f,dtype=np.int32,count=1)           #snapshot number
        self.nit    = np.fromfile(f,dtype=np.int32,count=1)           #number of interations completed
        self.time   = np.fromfile(f,dtype=float,count=1)*ts           #current time in code units
        self.nav    = np.fromfile(f,dtype=np.int32,count=1)           #chooses artificial viscosity scheme
        self.alpha  = np.fromfile(f,dtype=float,count=1)              #artificial viscosity parameter
        self.beta   = np.fromfile(f,dtype=float,count=1)              #artificial viscosity parameter
        self.tskip  = np.fromfile(f,dtype=float,count=1)              #??
        self.ngr    = np.fromfile(f,dtype=np.int32,count=1)           #gravity boundary conditions
        self.nrelax = np.fromfile(f,dtype=np.int32,count=1)           #relaxation flag
        self.trelax = np.fromfile(f,dtype=float,count=1)              #timescale for artificial drag force
        self.dt     = np.fromfile(f,dtype=float,count=1)              #current timestep
        self.omega2 = np.fromfile(f,dtype=float,count=1)              #square of the orbital velocity
        dum         = np.fromfile(f,dtype=np.int32,count=1)           #dummy variable
























