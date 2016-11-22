import numpy as np
import math



class Header:
  def __init__(self,ntot,nnopt,hmin,hmax,sep0,tf,dtout,nout,nit,time,nav,alpha,beta,tskip,ngr,nrelax,trelax,dt,omega2):
    self.ntot = ntot                                              #total number of particles
    self.nnopt = nnopt                                            #optimal number of neighbors
    self.hmin = hmin                                              #smoothing length for black holes
    self.hmax = hmax                                              #smoothing length for black holes
    self.sep0 = sep0                                              #initial separation of two stars
    self.tf = tf                                                  #time the code stops in code units
    self.dtout = dtout                                            #time interval between writing snapshots
    self.nout = nout                                              #numerical value of snapshot file
    self.nit = nit                                                #number of iterations completed
    self.time = time                                              #current time in years
    self.nav = nav                                                #chooses artificial viscosity scheme
    self.alpha = alpha                                            #artificial viscosity parameter
    self.beta = beta                                              #artificial viscosity parameter
    self.tskip = tskip                                            #??
    self.ngr = ngr                                                #gravity boundary conditions
    self.nrelax = nrelax                                          #relaxation flag
    self.trelax = trelax                                          #timescale for artificial drag force
    self.dt = dt                                                  #current timestep
    self.omega2 = omega2                                          #square of the orbital velocity
    























