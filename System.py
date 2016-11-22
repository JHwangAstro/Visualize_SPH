import numpy as np
import math
import Calcs


class System:

  def __init__(self,Planets,Star,):
    self.Particles = Particles

    self.M = np.sum([Particles[i].m for i in range(len(Particles))])
    self.N = len(Particles)

    self.vx = 0
    self.vy = 0
    self.vz = 0
    self.x = 0
    self.y = 0
    self.z = 0

    if self.M != 0:
      self.vx = np.sum([Particles[i].m*Particles[i].vx for i in range(len(Particles))])/self.M
      self.vy = np.sum([Particles[i].m*Particles[i].vy for i in range(len(Particles))])/self.M
      self.vz = np.sum([Particles[i].m*Particles[i].vz for i in range(len(Particles))])/self.M
      self.v = math.pow(self.vx**2+self.vy**2+self.vz**2,0.5)

      #Generate Profiles
      print 'Check Calculation Control: ', self.M
      self.x = Calcs.Find_COM(Particles) #Finds center of mass

      #Create Core
      self.core_com = Calcs.Find_COM([Particles[i] for i in range(len(Particles)) if Particles[i].rho > 0.4])

      print 'CHECK!: ', self.M, [self.x,self.y,self.z], self.core_com




class PlanetProf:

  def __init__(self,Planet):
    nbins = 100
    self.r,self.rho,self.m,self.T,self.P,self.u,self.gamm,self.cc = Calcs.Generate_Profile(nbins,[Planet.x,Planet.y,Planet.z],Planet.Particles)




























