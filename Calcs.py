import numpy as np
import struct
import math
import Calcs

#comc takes in a binary star configuration and calculates the center of masses and separation
#in:  star, ss, x, y, z, mass, M1, M2, rsep
#out: rsep, com, vcom
def Find_COM(P):
  #P - particles
  M  = sum([P[i].m        for i in range(len(P))])
  px = sum([P[i].x*P[i].m for i in range(len(P))])
  py = sum([P[i].y*P[i].m for i in range(len(P))])
  pz = sum([P[i].z*P[i].m for i in range(len(P))])

  return [px/M,py/M,pz/M]




#Converts cartesian coordinates to orbital elements
#comc takes in a binary star configuration and calculates the center of masses and separation
#in:  star, ss, x, y, z, mass, M1, M2, rsep
#out: rsep, com, vcom
def Find_VCOM(P):
  #P - particles
  M  = sum([P[i].m        for i in range(len(P))])
  px = sum([P[i].vx*P[i].m for i in range(len(P))])
  py = sum([P[i].vy*P[i].m for i in range(len(P))])
  pz = sum([P[i].vz*P[i].m for i in range(len(P))])

  return [px/M,py/M,pz/M]




#in:  r,v
#out: a,e,i,g,n,M
def Convert_To_Orbital_Elements(mu,r,v):

  r = np.array(r)
  v = np.array(v)
  eps = math.pow(10.,-6)

  h=np.cross(r,v)
  norm=np.cross([0,0,1],h)

  evec = ((np.linalg.norm(v)**2.-mu/np.linalg.norm(r))*r-np.dot(r,v)*v)/mu
  e = np.linalg.norm(evec)

  energy = np.linalg.norm(v)**2./2.-mu/np.linalg.norm(r)

  if math.fabs(e-1.0)>eps:
    a = -mu/(2.*energy)
    p = a*(1.-e**2.)
  else:
    a = math.inf
    p = np.linalg.norm(h)**2./mu

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

  return a,e,i,g,n,M








#Generate_Profile calculates values for each shell of a star
#in:  nbins,Particles
#out: Profile arrays
def Generate_Profile(nbins,COM,Particles):

#_________________________________________________________________CONSTANTS

  N     = len(Particles)
  kb    = 1.380658*math.pow(10,-16)                               #boltzmann constant in cgs
  arad  = 7.5657*math.pow(10,-15)                                 #boltzmann constant in cgs
  g     = 6.67390*math.pow(10,-8)                                 #gravitational constant [cgs]
  rs    = 6.9599*math.pow(10,10)                                  #unit of length [cm]
  ms    = 1.9891*math.pow(10,33)                                  #unit of mass [g]
  q     = 1.5*kb/arad

#_________________________________________________________________CALCULATES VARIABLES

  x  = np.array([Particles[i].x for i in range(N)])
  y  = np.array([Particles[i].y for i in range(N)])
  z  = np.array([Particles[i].z for i in range(N)])
  cc = np.array([Particles[i].cc for i in range(N)])
  mmw   = np.array([Particles[i].mmw for i in range(N)])
  ucgs  = np.array([Particles[i].u*g*ms/rs for i in range(N)])
  rhocgs= np.array([Particles[i].rho*ms/math.pow(rs,3) for i in range(N)])
  T = np.array([Calcs.get_T(q*rhocgs[i]/mmw[i],-ucgs[i]*rhocgs[i]/arad) for i in range(N)])

  pgas  = rhocgs*kb*T/mmw                                         #calculates the gas pressure
  prad  = arad*T**4/3.0                                           #calculates the radiation pressure
  ptot  = pgas + prad																							#calculates the total pressure
  beta  = pgas/ptot                                               #calculates the ratio of gas pressure to total pressure
  gamm	= (32.0-24.0*beta-3.0*beta**2)/(24.0-21.0*beta)           #calculates gamma as a function of the pressure ratio

  r = ((x-COM[0])**2+(y-COM[1])**2+(z-COM[2])**2)**0.5            #calculates radius of each particle

#_________________________________________________________________FINDS Z BIN CHARACTERISTICS
  rbin = np.array([(i)*max(r)/nbins for i in range(nbins)])
  index = np.array([np.where(r > rbin[i-1],1,0)*np.where(r <= rbin[i],1,0) for i in range(1,nbins)])
  snum  = np.array([len(np.compress(index[i-1],index[i-1])) for i in range(1,nbins)])
  rho_R = np.array([sum(np.compress(index[i-1],rhocgs))/snum[i-1] for i in range(1,nbins)])
  #m_R   = np.array([sum(np.compress(index[i],np.array([Particles[j].m for j in range(N)])))/snum[i] for i in range(1,nbins)])
  T_R   = np.array([sum(np.compress(index[i-1],T))/snum[i-1] for i in range(1,nbins)])
  u_R   = np.array([sum(np.compress(index[i-1],np.array([Particles[j].u for j in range(N)])))/snum[i-1] for i in range(1,nbins)])
  P_R   = np.array([sum(np.compress(index[i-1],ptot))/snum[i-1] for i in range(1,nbins)])
  gamm_R= np.array([sum(np.compress(index[i-1],gamm))/snum[i-1] for i in range(1,nbins)])
  cc_R  = np.array([sum(np.compress(index[i-1],cc))/snum[i-1] for i in range(1,nbins)])

  index = np.array([np.where(r <= rbin[i],1,0) for i in range(nbins)])
  m = np.array([Particles[j].m for j in range(N)])
  m_R = np.array([sum(np.compress(index[i],m))/sum(m) for i in range(1,nbins)]) #calculates the mass within each shell as a percentage
  return rbin[1:len(rbin)],rho_R,m_R,T_R,P_R,u_R,gamm_R,cc_R






#massbnd calculates which system each particle is bound to
#in:  ntot, x, y, z, vx, vy, vz, u, m, cor1i, cor2i, com, rcor
#out: mt, ss, MB
#def massbnd(ntot, x, y, z, m, h, vx, vy, vz, u, grpot, omega, cor1i, cor2i, com1, com2, rcor, mt, ss, MB, Ecom):
def Mass_Bound(Particles,Planets,Star):

  RS = 1.0

  N = len(Particles)

  m  = np.array([Particles[i].m  for i in range(N)])
  x  = np.array([Particles[i].x  for i in range(N)])
  y  = np.array([Particles[i].y  for i in range(N)])
  z  = np.array([Particles[i].z  for i in range(N)])
  vx = np.array([Particles[i].vx for i in range(N)])
  vy = np.array([Particles[i].vy for i in range(N)])
  vz = np.array([Particles[i].vz for i in range(N)])
  u  = np.array([Particles[i].u  for i in range(N)])
  grpot = np.array([Particles[i].grpot for i in range(N)])

  #_________________________________________________________________CONSTANTS
  G = 1.0

  #calculate the distances, velocities and energies from each particle to the respective cores
  dx = np.array([((x-Planets[i].com[0])**2+(y-Planets[i].com[1])**2+(z-Planets[i].com[2])**2)**0.5 for i in range(len(Planets))])
  dv = np.array([((vx-Planets[i].vx)**2 +(vy-Planets[i].vy)**2 +(vz-Planets[i].vz)**2)**0.5  for i in range(len(Planets))])
  E  = np.array([m*(dv[i]**2./2.0-G*Planets[i].M/dx[i]) for i in range(len(Planets))])

  #calculate the distances, velocities and energies of each particle
  r  = ( x**2+ y**2+ z**2)**0.5
  dv2= (vx**2+vy**2+vz**2)**0.5
  Ecom = m*(dv2**2/2.0+grpot)

  #assign mass bound numbers
  MB = np.zeros(N)
  MB[:] = 4
  MB[np.compress(np.where(Ecom < 0.0,1,0),np.arange(N))] = 2        #Bound to the system
  if len(Planets) == 2:
    MB[np.compress(np.where((E[1] < 0.0) & (E[1] < E[0]),1,0),np.arange(N))] = 1 #Bound to planet 2
    MB[np.compress(np.where((E[0] < 0.0) & (E[0] < E[1]),1,0),np.arange(N))] = 0 #Bound to planet 1
  if len(Planets) == 1:
    MB[np.compress(np.where(E[0] < 0.0,1,0),np.arange(N))] = 0 #Bound to planet 1
  MB[N-1] = 5                                                       #Index of star
  #MB[np.compress(np.where(r < RS,1,0),np.arange(N))] = 3            #Accreted by the star


  return MB
  #0 - Bound to planet 1
  #1 - Bound to planet 2
  #2 - Bound to the system
  #3 - Accreted
  #4 - Ejecta





def get_T(q, r):
#Subroutine to solve 4th order equations to determine the temperature x3
#for an equation of state with both ideal gas and radiation pressure.
#Written by Scott Fleming 10/04/02 and James Lombardi 2002-2003

#The fourth order equation comes from u_gas+ u_rad = u, with
#u_gas proportional to T and u_rad proportional to T^4

#In general, we can transform a 4th order equation to x^4+px^2+qx+r=0
#(see pages 57-58 of Stillwell's "Mathematics and its history" text)
#but we fortunately don't even have an x^2 term (that is, p=0).

#Follow Stillwell, we can transform this into a cubic equation:
#First solve for y by using the fact that B^2-4AC=0
#equation is then:  y^3=ry+q^2/8
#using the solution of cubic equations found in Stillwell page 55:

  k   = 0.125*math.pow(q,2)
  kh  = 0.5*k

  p1 = kh + math.sqrt(math.pow(kh,2) - math.pow((r/3.0),3))
  p2 = math.pow((r/3.0),3)/p1

  y1 = math.pow(p1,1.0/3.0)

#Fortran can't handle cube roots of neg. #'s
  y2 = -math.pow(math.fabs(p2),1.0/3.0)
  yy = y1 + y2

#Equation to solve: (x^2+p+y)^2=Ax^2+Bx+C
#Now take square root of both sides with:

  AA = 2.0*yy
  B  = -q

#Re-writing Ax^2+Bx+C as a square then solving the equation we
#obtain 2 results:
#x^2 + (-(A^(1/2)))x + (-B/(2(A)^(1/2))+p+y) = 0 (1)
#or
#x^2 + (A^(1/2))x + (B/(2(A)^(1/2))+p+y) = 0     (2)

#Our solution we're interested in:
  if AA < 0:
    print 'error in temperature subroutine'
    return 0.0
  b2 = math.sqrt(AA)
  c2 = 0.50*B/b2 + yy

#Therefore, we once again have x^2+bx+c=0, and our answer we want is:
  x3 = -2.0*c2/(b2 + math.sqrt(math.pow(b2,2)-4.0*c2))
  if b2 == 0.0: x3 = math.pow(-r - q*math.pow(-r - q*math.pow(-r,0.25),0.25),0.25)
  if p2 >= 0.0: x3 = -(r+math.pow(r/q,4))/q

  return x3










