import numpy as np
import struct
import math
import ReadData
import pickle
from Header import Header
from Particle import Particle


#Read_Snapshot imports the data from each snapshot and outputs an array of particle objects and a header object
#in:  sdir, ss, ts
#out: exts, t, trel, omega2, Particles
def Read_Snapshot(sdir,ss):

  #_______________________________________________________________SETS FILENAME FOR SNAPSHOT
  fname = sdir + 'out000' + repr(ss) + '.sph'
  if ss >= 10:   fname = sdir + 'out00' + repr(ss) + '.sph'
  if ss >= 100:  fname = sdir + 'out0'  + repr(ss) + '.sph'
  if ss >= 1000: fname = sdir + 'out'   + repr(ss) + '.sph'

  #_______________________________________________________________READS DATA FILE(S)
  f = open(fname, 'rb')
  Header = ReadData.Read_Header(f)
  Particles = [ReadData.get_particle(i,f) for i in range(Header.ntot)]
  f.close()

  return Header, Particles




def Read_Header(f):
  g   = 6.67390*math.pow(10,-8)                                   #gravitational constant [cgs]
  rs  = 6.9599*math.pow(10,10)                                    #unit of length [cm]
  ms  = 1.9891*math.pow(10,33)                                    #unit of mass [g]
  ts  = math.sqrt(math.pow(rs,3)/(g*ms))/31556926.                #unit of time [yr]
  print 'ts: ', ts

  dum       = struct.unpack('i',f.read(4))[0]                     #dummy variable
  ntot      = struct.unpack('i',f.read(4))[0]                     #total number of particles
  nnopt     = struct.unpack('i',f.read(4))[0]                     #optimal number of neighbors
  hmin      = struct.unpack('d',f.read(8))[0]                     #smoothing length for black holes
  hmax      = struct.unpack('d',f.read(8))[0]                     #smoothing length for black holes
  sep0      = struct.unpack('d',f.read(8))[0]                     #initial separation of two stars
  tf        = struct.unpack('d',f.read(8))[0]                     #time the code stops in code units
  dtout     = struct.unpack('d',f.read(8))[0]                     #time interval between writing outxxx.sph files
  nout      = struct.unpack('i',f.read(4))[0]                     #xxx value of outputxxx.sph
  nit       = struct.unpack('i',f.read(4))[0]                     #number of iterations completed
  time      = struct.unpack('d',f.read(8))[0]*ts                  #current time in code units
  nav       = struct.unpack('i',f.read(4))[0]                     #chooses artificial viscosity scheme
  alpha     = struct.unpack('d',f.read(8))[0]                     #artificial viscosity parameter
  beta      = struct.unpack('d',f.read(8))[0]                     #artificial viscosity parameter
  tskip     = struct.unpack('d',f.read(8))[0]                     #??
  ngr       = struct.unpack('i',f.read(4))[0]                     #gravity boundary conditions
  nrelax    = struct.unpack('i',f.read(4))[0]                     #relaxation flag
  trelax    = struct.unpack('d',f.read(8))[0]                     #timescale for artificial drag force
  dt        = struct.unpack('d',f.read(8))[0]                     #current timestep
  omega2    = struct.unpack('d',f.read(8))[0]                     #square of the orbital velocity
  dum       = struct.unpack('i',f.read(4))[0]                     #dummy variable

  return Header(ntot,nnopt,hmin,hmax,sep0,tf,dtout,nout,nit,time,nav,alpha,beta,tskip,ngr,nrelax,trelax,dt,omega2)




def Read_Energy_Relax(f):
  t,Eg,Ek,Ei,E,S,A = np.loadtxt(f+'energy0.sph',skiprows=0,unpack=True,usecols=(0,1,2,3,4,5,6))
  print t
  return [t,E,Eg,Ek,Ei]




def Get_Core_Constant(f):
  #rho_e,rho_c,rho_m,rho_i,K_e,gamma_e,u_factor = np.loadtxt(f,skiprows=0,unpack=True,usecols=(0,1,2,3,4,5,6))
  rho_e,rho_c,rho_m,rho_i,K_e,gamma_e = np.loadtxt(f,skiprows=0,unpack=True,usecols=(0,1,2,3,4,5))
  return rho_c




def get_particle(i,f):
  dum      = struct.unpack('i',f.read(4))[0]                    #dummy variable
  x        = struct.unpack('d',f.read(8))[0]                    #x position
  y        = struct.unpack('d',f.read(8))[0]                    #y position
  z        = struct.unpack('d',f.read(8))[0]                    #z position
  m        = struct.unpack('d',f.read(8))[0]                    #mass of the particle
  h        = struct.unpack('d',f.read(8))[0]                    #smoothing length
  rho      = struct.unpack('d',f.read(8))[0]                    #density
  vx       = struct.unpack('d',f.read(8))[0]                    #x-velocity
  vy       = struct.unpack('d',f.read(8))[0]                    #y-velocity
  vz       = struct.unpack('d',f.read(8))[0]                    #z-velocity
  vxdot    = struct.unpack('d',f.read(8))[0]                    #x-acceleration
  vydot    = struct.unpack('d',f.read(8))[0]                    #y-acceleration
  vzdot    = struct.unpack('d',f.read(8))[0]                    #z-acceleration
  u        = struct.unpack('d',f.read(8))[0]                    #potential
  udot     = struct.unpack('d',f.read(8))[0]                    #change in potential
  grpot    = struct.unpack('d',f.read(8))[0]                    #gravitational potential
  mmw      = struct.unpack('d',f.read(8))[0]                    #mean molecular weight
  cc       = struct.unpack('i',f.read(4))[0]                    #??
  divv     = struct.unpack('d',f.read(8))[0]                    #divergence of velocities
  dum      = struct.unpack('i',f.read(4))[0]                    #dummy variable

  return Particle(i,x,y,z,m,h,rho,vx,vy,vz,vxdot,vydot,vzdot,u,udot,grpot,mmw,cc,divv)





def ReadAEI(fname):

  n_elements=14
  x = [[] for i in range(n_elements)]

  f = open(fname,'rb')
  dstr = f.readline()                               #reads in dummy lines
  dstr = f.readline()
  dstr = f.readline()
  dstr = f.readline()

  for line in f:
    flag = 0
    dummy = line.split()                            #reads in final orbital parameters
    for i in range(len(dummy)):
      if dummy[i] == '***************': flag = 1
      if dummy[i] == '********': flag = 1
    if flag == 0:
      for i in range(len(dummy)):
        x[i].append(float(dummy[i]))
    if flag == 1:
      x[0].append(float(dummy[0]))
      for i in range(1,len(dummy)):
        x[i].append(0.0)
      x[7][len(x[7])-1] = dummy[7]
  f.close()

  return np.array(x[0]),np.array(x[1]),np.array(x[2])





def restart(sdir,numss):

  time       = np.zeros(numss)
  Energy     = np.zeros(numss)
  Specific_E = np.zeros([2,numss])
  Specific_h = np.zeros([2,numss])
  Mgas       = np.zeros([2,numss])
  Mass       = np.zeros([4,numss])
  aei        = np.zeros([2,6,numss])
  aei2       = np.zeros([2,6,numss])

  with open(sdir+'sph.restart_analysis', 'rb') as input:
    time0       = pickle.load(input)
    Energy0     = pickle.load(input)
    Specific_E0 = pickle.load(input)
    Specific_h0 = pickle.load(input)
    Mgas0       = pickle.load(input)
    Mass0       = pickle.load(input)
    aei0        = pickle.load(input)
    aei20       = pickle.load(input)
    Planets     = pickle.load(input)
    Planets0    = pickle.load(input)
    Star        = pickle.load(input)
    ninit       = pickle.load(input)

  for i in range(0,len(time)):
    time[i]         = time0[i]
    Energy[i]       = Energy0[i]
    Specific_E[0,i] = Specific_E0[0,i]
    Specific_E[1,i] = Specific_E0[1,i]
    Specific_h[0,i] = Specific_h0[0,i]
    Specific_h[1,i] = Specific_h0[1,i]
    Mgas[0,i]       = Mgas0[0,i]
    Mgas[1,i]       = Mgas0[1,i]
    Mass[0,i]       = Mass0[0,i]
    Mass[1,i]       = Mass0[1,i]
    Mass[2,i]       = Mass0[2,i]
    Mass[3,i]       = Mass0[3,i]
    for j in range(0,6):
      aei[0,j,i]    = aei0[0,j,i]
      aei[1,j,i]    = aei0[1,j,i]
      aei2[0,j,i]   = aei20[0,j,i]
      aei2[1,j,i]   = aei20[1,j,i]

  return time,Energy,Specific_E,Specific_h,Mgas,Mass,aei,aei2,Planets,Planets0,Star,ninit












