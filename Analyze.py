import numpy as np
import math
import ReadData
import WriteData
import Plots
import Analyze
import Calcs
import Planet
import os.path
from Star   import Star
from Planet import Planet
from Planet import PlanetProf
import matplotlib.pyplot as plt

#Planet_Planet_Collision loops over the snapshots and calls the modules to read and plot the data for a planet planet collision simualtion
#in: sdir,bsdir,idir,fdir,numss
def Planet_Planet_Collision(sdir,idir,fdir,mdir,nstart,numss,nskip,restart,np2):
  #_________________________________________________________________VARIABLES

  MS = 1.98855*math.pow(10.,33.)                                    #Solar Mass in g
  ME = 5.972*math.pow(10.,27.)
  G  = 1.0

  #Variables tracked per snapshot
  time      = np.zeros(numss)
  Energy    = np.zeros(numss)
  Specific_E= np.zeros([2,numss])
  Specific_h= np.zeros([2,numss])
  Mgas      = np.zeros([2,numss])
  Mass      = np.zeros([4,numss])
  aei       = np.zeros([2,6,numss])
  aei_b     = np.zeros([2,6,numss])
  ninit     = nstart

  if os.path.isfile(sdir+'sph.restart_analysis') and restart == 1:
    time,Energy,Specific_E,Specific_h,Mgas,Mass,aei,aei_b,Planets,Planets0,Star,ninit = ReadData.restart(sdir,numss)

  for i in range(ninit,numss*nskip,nskip):
    print 'Snapshot: ', i,nskip
    Header,Particles = ReadData.Read_Snapshot(sdir,i)
    #Create objects for first snapshot
    if i == nstart:
      Planets,Star = Analyze.Gen_System_Initial(Header,Particles,(len(Particles)-np2-1)/2,(len(Particles)-np2-1)/2)
      Planets0 = Planets
      Eaten    = [Planets[0].Particles[0]]
      Envelope = [Planets[0].Particles[0]]
      Ejecta   = [Planets[0].Particles[0]]
    #create objects if not first snapshot
    if i !=nstart: Planets,Envelope,Eaten,Ejecta = Analyze.Gen_System(Header,Particles,Planets,Star)

    for j in range(2):
      aei[j,0,i/nskip] = Planets[j].a*0.00465
      aei[j,1,i/nskip] = Planets[j].e
      aei[j,2,i/nskip] = Planets[j].i
      aei[j,3,i/nskip] = Planets[j].g
      aei[j,4,i/nskip] = Planets[j].n
      aei[j,5,i/nskip] = Planets[j].nu

    #tracks the mass of each type of categorization in ME
    Mass[0,i/nskip] = Planets[0].M*MS/ME
    Mass[1,i/nskip] = Planets[1].M*MS/ME
    Mass[2,i/nskip] = sum([Envelope[j].m for j in range(len(Envelope))])*MS/ME
    Mass[3,i/nskip] = sum([Ejecta[j].m for j in range(len(Ejecta))])*MS/ME

    #If there are two planets, track the energy
    if np.min([Planets[0].M,Planets[1].M]) > 0:
      v_cm  = np.array([(Planets[0].M*Planets[0].vel[j]+Planets[1].M*Planets[1].vel[j])/(Planets[0].M+Planets[1].M) for j in range(3)])
      v1_cm = np.array([Planets[0].vel[j]-v_cm[j] for j in range(3)])
      v2_cm = np.array([Planets[1].vel[j]-v_cm[j] for j in range(3)])
      v1_cm2= np.sum(np.array([v1_cm[j]**2. for j in range(3)]))
      v2_cm2= np.sum(np.array([v2_cm[j]**2. for j in range(3)]))
      r  = np.sqrt(np.sum(np.array([(Planets[1].com[j]-Planets[0].com[j])**2. for j in range(3)])))
      K_cm  = 0.5*(Planets[0].M*v1_cm2+Planets[1].M*v2_cm2)
      P_cm  = -G*Planets[0].M*Planets[1].M/r
      Energy[i/nskip] = K_cm+P_cm
    Specific_E[0,i/nskip] = -0.5*((Star.M*Planets[0].M)/Planets[0].a)
    Specific_E[1,i/nskip] = -0.5*((Star.M*Planets[1].M)/Planets[1].a)
    Specific_h[0,i/nskip] = Planets[0].M*Star.M/(Planets[0].M+Star.M)*math.sqrt(Planets[0].a*(1.-Planets[0].e**2.)*(Star.M+Planets[0].M))
    Specific_h[1,i/nskip] = Planets[1].M*Star.M/(Planets[1].M+Star.M)*math.sqrt(Planets[1].a*(1.-Planets[1].e**2.)*(Star.M+Planets[1].M))

    Core0 = Planet(Header,[Planets[0].Particles[j] for j in range(len(Planets[0].Particles)) if Planets[0].Particles[j].rho > ReadData.Get_Core_Constant(sdir+'sph.rhoparams')], Star)
    Core1 = Planet(Header,[Planets[1].Particles[j] for j in range(len(Planets[1].Particles)) if Planets[1].Particles[j].rho > ReadData.Get_Core_Constant(sdir+'sph.rhoparams2')],Star)

    #Calculate gas mass of each planet
    #Mgas[0,i] = 1.-(Core0.M/Planets[0].M)
    #Mgas[1,i] = 1.-(Core1.M/Planets[1].M)
    Mgas[0,i/nskip] = 1.-((Mass[0,0]-Mgas[0,0])/Planets[0].M)
    Mgas[1,i/nskip] = 1.-((Mass[1,0]-Mgas[1,0])/Planets[1].M)

    #Plot snapshot data
    Plots.Plot_Snapshot_Data(idir,i,Planets,Envelope,Ejecta,Star)

    #if option is enabled, then make two sets of snapshots where the origin is at a planet
    if 1:
      if Core0.M > 0.0: WriteData.Origin_Shift(Core0.com,[0.,0.,0.],sdir,i,Header,Particles,'Inner_')
      if Core1.M > 0.0: WriteData.Origin_Shift(Core1.com,[0.,0.,0.],sdir,i,Header,Particles,'Outer_')
      if Energy[i/nskip] < 0:
        #Do calculations for bound planets
        COM = [(Core0.com[k] *Core0.M+Core1.com[k] *Core1.M)/(Core0.M+Core1.M) for k in range(3)]
        VCOM= [(Core0.vcom[k]*Core0.M+Core1.vcom[k]*Core1.M)/(Core0.M+Core1.M) for k in range(3)]
        a_temp,e_temp,i_temp,g_temp,n_temp,nu_temp = Calcs.Convert_To_Orbital_Elements(Planets[0].M+Planets[1].M,Planets[0].com,[Core0.vcom[k]-VCOM[k] for k in range(3)])
        aei_b[0,0,i/nskip] = a_temp
        aei_b[0,1,i/nskip] = e_temp
        a_temp,e_temp,i_temp,g_temp,n_temp,nu_temp = Calcs.Convert_To_Orbital_Elements(Planets[0].M+Planets[1].M,Planets[1].com,[Core1.vcom[k]-VCOM[k] for k in range(3)])
        aei_b[1,0,i/nskip] = a_temp
        aei_b[1,1,i/nskip] = e_temp
        rotation = [0.,i_temp,0.]
        WriteData.Origin_Shift(COM,rotation,sdir,i,Header,Particles,'Bound_')
    P0 = aei[0,0,0]**1.5
    time[i/nskip] = Header.time/P0                                           #convert years to orbital periods of innermost orbit

  Plots.Summary_Plots(sdir,fdir,time,aei,Mass,Energy,Star.M,Specific_E,Specific_h,Mgas,aei_b)
  for i in range(len(Planets)):
    if Planets[i/nskip].M > 0.0: Plots.Profile(PlanetProf(Planets[i/nskip]),PlanetProf(Planets0[i/nskip]),i,fdir)

  WriteData.Write_Summary(fdir,Planets0,Planets,Star,Mass,Energy,Specific_E,Specific_h,Mgas)
  WriteData.Write_Mercury(sdir,mdir,Planets)
  WriteData.Write_Restart(sdir,time,Energy,Specific_E,Specific_h,Mgas,Mass,aei,aei_b,Planets,Planets0,Star,numss)




def Planet_Profile(sdir,idir,fdir,numss):
  for i in range(0,numss):
    #read snapshot and generate planet, core, and planet profile objects
    Header,Particles = ReadData.Read_Snapshot(sdir,i)
    Star = 0.
    Planet1 = Planet(Header,Particles,Star)
    Core = Planet(Header,[Particles[j] for j in range(len(Particles)) if Particles[j].cc == 2],Star)
    PlanetProf1 = PlanetProf(Planet1)

    print 'Snapshot ', i
    print 'Particles in core: ', len(Core.Particles)
    print 'Coordinates of Core:', Core.x, Core.y, Core.z

    if i == 0: PlanetProf0 = PlanetProf(Planet1)
    Plots.Profile(PlanetProf0,PlanetProf1,i,idir)
  Plots.Energy(ReadData.Read_Energy_Relax(sdir),fdir)




def Gen_System_Initial(Header,Particles,N1,N2):
  #assuming that the particles are assigned in index order to the planets
  HostStar = Star(Particles[len(Particles)-1])
  return [Planet(Header,Particles[0:N1-1],HostStar),Planet(Header,Particles[N1:N1+N2-1],HostStar)],HostStar




def Gen_System(Header,Particles,Planets,Star):
  #Update Jupiters with new particle properties
  Planets = [Planet(Header,[Particles[i] for i in [Planets[j].Particles[k].i for k in range(Planets[j].N)]],Star) for j in range(len(Planets))]
  Index = Calcs.Mass_Bound(Particles,[Planets[i] for i in range(len(Planets)) if Planets[i].M > 0.0],Star)
  print 'particles belonging to planet1: ', len([Particles[i] for i in range(len(Particles)) if Index[i]==0])
  print 'particles belonging to planet2: ', len([Particles[i] for i in range(len(Particles)) if Index[i]==1])
  print 'particles in the envelope:      ', len([Particles[i] for i in range(len(Particles)) if Index[i]==2])
  print 'particles accreted by star      ', len([Particles[i] for i in range(len(Particles)) if Index[i]==3])
  print 'particles ejected from system:  ', len([Particles[i] for i in range(len(Particles)) if Index[i]==4])
  print 'just to double check: ', len(Index), len(Particles)
  Planet1 = Planet(Header,np.array([Particles[i] for i in range(len(Particles)) if Index[i] == 0]),Star)
  Planet2 = Planet(Header,np.array([Particles[i] for i in range(len(Particles)) if Index[i] == 1]),Star)
  Envelope= [Particles[i] for i in range(len(Particles)) if Index[i] == 2]
  Eaten   = [Particles[i] for i in range(len(Particles)) if Index[i] == 3]
  Ejecta  = [Particles[i] for i in range(len(Particles)) if Index[i] == 4]
  return [[Planet1,Planet2],Envelope,Eaten,Ejecta]


















