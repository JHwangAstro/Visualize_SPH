import numpy as np
import math
import Calcs


class Planet:

  def __init__(self,Header,Particles,Star):
    self.Particles = Particles

    self.M = np.sum([Particles[i].m for i in range(len(Particles))])
    self.N = len(Particles)

    self.vx = 0
    self.vy = 0
    self.vz = 0
    self.x = 0
    self.y = 0
    self.z = 0
    self.a = 0
    self.e = 0
    self.i = 0
    self.g = 0
    self.n = 0
    self.nu = 0

    if self.M != 0:
      self.vx = np.sum([Particles[i].m*Particles[i].vx for i in range(len(Particles))])/self.M
      self.vy = np.sum([Particles[i].m*Particles[i].vy for i in range(len(Particles))])/self.M
      self.vz = np.sum([Particles[i].m*Particles[i].vz for i in range(len(Particles))])/self.M
      self.v = math.pow(self.vx**2+self.vy**2+self.vz**2,0.5)

      #Generate Profiles
      self.com = Calcs.Find_COM(Particles) #Finds center of mass
      self.vcom= Calcs.Find_VCOM(Particles)#Finds center of mass velocity

      if Star != 0.0:
        G = 1.0
        self.mu = G*(Star.M+self.M)
        self.a,self.e,self.i,self.g,self.n,self.nu = Calcs.Convert_To_Orbital_Elements(self.mu,self.com,[self.vx,self.vy,self.vz])
        self.E = -self.mu/(2.*self.a)
        self.h = math.sqrt(self.a*(1.-self.e**2.)*self.mu)
    self.vel= [self.vx,self.vy,self.vz]





class PlanetProf:

  def __init__(self,Planet):
    nbins = 100
    self.r,self.rho,self.m,self.T,self.P,self.u,self.gamm,self.cc = Calcs.Generate_Profile(nbins,Planet.com,Planet.Particles)




























