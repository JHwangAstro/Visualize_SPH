import numpy as np
import math

class Header:
    def __init__(self,f):

#        g   = 6.67390*math.pow(10,-8)                               # gravitational constant [cgs]
#        rs  = 6.9599*math.pow(10,10)                                # unit of length [cm]
#        ms  = 1.9891*math.pow(10,33)                                # unit of mass [g]
#        ts  = math.sqrt(math.pow(rs,3)/(g*ms))/31556926.            # unit of time [yr]
        ts  = 5.05000318388e-05                                     # Convert from code units to years

        dtype = np.dtype("i4,i4,i4,f8,f8,f8,f8,f8,i4,i4,f8,i4, \
                f8,f8,f8,i4,i4,f8,f8,f8,i4")

        [(dum,self.ntot,self.nnopt,self.hmin,self.hmax,self.sep0, \
                self.tf,self.dtout,self.nout,self.nit,self.time, \
                self.nav,self.alpha,self.beta,self.tskip,self.ngr, \
                self.nrelax,self.trelax,self.dt,self.omega2,dum)] = \
                np.fromfile(f,dtype=dtype,count=1)

        self.time = self.time*ts

























