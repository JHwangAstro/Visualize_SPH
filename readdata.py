import numpy as np
import struct
import math
import readdata
import pickle
from header import Header
from particle import Particle



def Read_Energy_Relax(f):
  t,Eg,Ek,Ei,E,S,A = np.loadtxt(f+'energy0.sph',skiprows=0,unpack=True,usecols=(0,1,2,3,4,5,6))
  return [t,E,Eg,Ek,Ei]




def Get_Core_Constant(f):
  #rho_e,rho_c,rho_m,rho_i,K_e,gamma_e,u_factor = np.loadtxt(f,skiprows=0,unpack=True,usecols=(0,1,2,3,4,5,6))
  rho_e,rho_c,rho_m,rho_i,K_e,gamma_e = np.loadtxt(f,skiprows=0,unpack=True,usecols=(0,1,2,3,4,5))
  return rho_c




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












