import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import GenPlots
import ReadData
import os

def Profile(PlanetProf0,PlanetProf1,i,fdir):

  mask0 = np.where(PlanetProf0.T != 0,1,0)
  mask1 = np.where(PlanetProf1.T != 0,1,0)
  GenPlots.Double_Plot(PlanetProf0.r,PlanetProf1.r,PlanetProf0.T,  PlanetProf1.T,  mask0,mask1,r'$r/R$',r'$T\ \mathrm{[10^6K]}$',fdir+'Profile_T'+str(i),0)
  GenPlots.Double_Plot(PlanetProf0.r,PlanetProf1.r,PlanetProf0.rho,PlanetProf1.rho,mask0,mask1,r'$r/R$',r'$\rho\ \mathrm{[}M_\odot R_\odot^{-3}\mathrm{]}$',fdir+'Profile_rho'+str(i),1)
  GenPlots.Double_Plot(PlanetProf0.r,PlanetProf1.r,PlanetProf0.P,  PlanetProf1.P,  mask0,mask1,r'$r/R$',r'$P\ \mathrm{[Ba]}$',fdir+'Profile_P'+str(i),1)
  GenPlots.Double_Plot(PlanetProf0.r,PlanetProf1.r,PlanetProf0.m,  PlanetProf1.m,  mask0,mask1,r'$r/R$',r'$m(r)/M$',fdir+'Profile_m'+str(i),0)
  #GenPlots.Double_Plot(PlanetProf0.r,PlanetProf1.r,PlanetProf0.cc, PlanetProf1.cc, mask0,mask1,r'$r/R$',r'$x_\mathrm{n}$',fdir+'Profile'+str(i),0)




def Plot_Snapshot_Data(idir,ss,Jupiters,Envelope,Ejecta,Star):

  x = np.array([np.array([Jupiters[i].Particles[j].x for j in range(Jupiters[i].N)]) for i in range(len(Jupiters))])
  y = np.array([np.array([Jupiters[i].Particles[j].y for j in range(Jupiters[i].N)]) for i in range(len(Jupiters))])
  x1 = np.array([Envelope[i].x for i in range(len(Envelope))])
  x2 = np.array([Ejecta[i].x   for i in range(len(Ejecta))])
  y1 = np.array([Envelope[i].y for i in range(len(Envelope))])
  y2 = np.array([Ejecta[i].y   for i in range(len(Ejecta))])
  #z = np.array([np.array([Jupiters[i].Particles[j].z for j in range(Jupiters[i].N)]) for i in range(len(Jupiters))])
  #print x.shape
  #GenPlots.Particle_Plot(x,y,x1,y1,x2,y2,r'$x/R_\mathrm{*}$',r'$y/R_\mathrm{*}$',idir+'XY'+str(ss))
  #GenPlots.Particle_Plot(x,z,r'$x/R_\mathrm{*}$',r'$z/R_\mathrm{*}$',idir+'XZ'+str(ss))
  #GenPlots.Particle_Plot(y,z,r'$y/R_\mathrm{*}$',r'$z/R_\mathrm{*}$',idir+'YZ'+str(ss))
  #GenPlots.Zoomed_Particle_Plot(x,y,x1,y1,x2,y2,r'$x/R_\mathrm{*}$',r'$y/R_\mathrm{*}$',idir+'Zoom_XY'+str(ss))
  #GenPlots.Zoomed_Particle_Plot(x,z,r'$x/R_\mathrm{*}$',r'$z/R_\mathrm{*}$',idir+'Zoom_XZ'+str(ss))
  #GenPlots.Zoomed_Particle_Plot(y,z,r'$y/R_\mathrm{*}$',r'$z/R_\mathrm{*}$',idir+'Zoom_YZ'+str(ss))




def Summary_Plots(sdir,fdir,time,orbital_elements,Mass,Energy,MStar,E,h,Mgas,bound_orbital_elements):

  c = ['b','r']

  N = len(time)
  #CREATE ORBITS PLOT (a,e)
  for i in range(2):
    t = np.array([time[j] for j in range(N) if orbital_elements[i,0,j] != 0.0])
    a = np.array([orbital_elements[i,0,j] for j in range(N) if orbital_elements[i,0,j] != 0.0]) #conversion to AU
    e = np.array([orbital_elements[i,1,j] for j in range(N) if orbital_elements[i,0,j] != 0.0])
    plt.plot(t, a, color = c[i], linestyle = '-')  #semi-major axes
    plt.plot(t, (1.-e)*a, color = c[i], linestyle = ':')  #periapsis
    plt.plot(t, (1.+e)*a, color = c[i], linestyle = ':')  #apoapsis
  plt.xlim([0,max(time)])
  plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
  plt.ylabel(r'$a\ [\mathrm{AU}]$', fontsize=15)

  fname = fdir + 'orbits.png'
  plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()

  #CREATE BOUND ORBITS PLOT (a,e)
  print 'length: ', len(Energy[Energy<0.])
  if len(Energy[Energy<0.]) > 1:
    for i in range(2):
      N = len(time)
      t = np.array([time[j] for j in range(N) if Energy[j] < 0.])
      a = np.array([bound_orbital_elements[i,0,j] for j in range(N) if Energy[j] < 0.])
      e = np.array([bound_orbital_elements[i,1,j] for j in range(N) if Energy[j] < 0.])
      plt.plot(t, a, color = c[i], linestyle = '-')  #semi-major axes
      plt.plot(t, (1.-e)*a, color = c[i], linestyle = ':')  #periapsis
      plt.plot(t, (1.+e)*a, color = c[i], linestyle = ':')  #apoapsis
    plt.xlim([0,max(time)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$a\ [\mathrm{AU}]$', fontsize=15)

    fname = fdir + 'bound_orbits.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

  if os.path.isfile(sdir+'Kepler0.aei') and os.path.isfile(sdir+'Kepler1.aei'):
    tM0,aM0,eM0 = ReadData.ReadAEI(sdir+'Kepler0.aei')
    tM1,aM1,eM1 = ReadData.ReadAEI(sdir+'Kepler1.aei')

    #for i in range(2):
    #  N = len(time)
    #  t = np.array([time[j] for j in range(N) if orbital_elements[i,0,j] != 0.0])
    #  a = np.array([orbital_elements[i,0,j] for j in range(N) if orbital_elements[i,0,j] != 0.0])
    #  e = np.array([orbital_elements[i,1,j] for j in range(N) if orbital_elements[i,0,j] != 0.0])
    #  plt.plot(t, a, color = c[i], linestyle = '-')  #semi-major axes
    #  plt.plot(t, (1.-e)*a, color = c[i], linestyle = ':')  #periapsis
    #  plt.plot(t, (1.+e)*a, color = c[i], linestyle = ':')  #apoapsis
    plt.plot(tM0, aM0, color = c[0], linestyle = '-')  #semi-major axes
    plt.plot(tM0, (1.-eM0)*aM0, color = c[0], linestyle = ':')  #periapsis
    plt.plot(tM0, (1.+eM0)*aM0, color = c[0], linestyle = ':')  #apoapsis
    plt.plot(tM1, aM1, color = c[1], linestyle = '-')  #semi-major axes
    plt.plot(tM1, (1.-eM1)*aM1, color = c[1], linestyle = ':')  #periapsis
    plt.plot(tM1, (1.+eM1)*aM1, color = c[1], linestyle = ':')  #apoapsis
    plt.xlim([0,max([max(t),max(tM0),max(tM1)])])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$a\ [\mathrm{AU}]$', fontsize=15)

    fname = fdir + 'orbits_extended.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()

  t0 = np.array([time[j]   for j in range(N) if Mass[0,j] != 0.0])
  t1 = np.array([time[j]   for j in range(N) if Mass[1,j] != 0.0])
  M0 = np.array([Mass[0,j] for j in range(N) if Mass[0,j] != 0.0])
  M1 = np.array([Mass[1,j] for j in range(N) if Mass[1,j] != 0.0])
  plt.plot(t0,M0-M0[0], color = 'k', linestyle = '-')       #mass in planet 1
  plt.plot(t1,M1-M1[0], color = 'k', linestyle = '--')      #mass in planet 2
  plt.plot(time, Mass[2,:], color = 'k', linestyle = ':')   #mass bound to system
  plt.plot(time, Mass[3,:], color = 'k', linestyle = '-.')  #mass in ejecta
  #for i in range(0,len(M0)-100,10):
  #  print time[i],M0[i]-M0[0], M1[i]-M1[0], M0[0], M1[0]
  #print min(M0-M0[0]), min(M1-M1[0]),max(M0-M0[0]),max(M1-M1[0])
  plt.axis([0,max(time),min(min(M0-M0[0]),min(M1-M1[0])),max(max(M0-M0[0]),max(M1-M1[0]))])
  #plt.xlim([0,max(t)])
  plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
  plt.ylabel(r'$\dot M/\mathrm{M}_\oplus$', fontsize=15)

  fname = fdir + 'mass_transfer.png'
  plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()

  plt.plot(time,Energy, color = 'k', linestyle = '-')  #mass in planet 1
  plt.plot([0,max(t)],[0.,0.],linestyle='--')
  plt.xlim([0,max(t)])
  plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
  plt.ylabel(r'$E$', fontsize=15)

  fname = fdir + 'planet_energy.png'
  plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()

  plt.plot(time,E[0,:], color = 'k', linestyle = '--')  #mass in planet 1
  plt.plot(time,E[1,:], color = 'k', linestyle = '-.')  #mass in planet 1
  plt.plot(time,E[0,:]+E[1,:], color = 'k', linestyle = '-')  #mass in planet 1
  plt.xlim([0,max(time)])
  plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
  plt.ylabel(r'$E$', fontsize=15)

  fname = fdir + 'specific_orbital_energy.png'
  plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()

  plt.plot(time,Mgas[0,:], color = 'k', linestyle = '--')  #mass in planet 1
  plt.plot(time,Mgas[1,:], color = 'k', linestyle = '-')   #mass in planet 2
  plt.xlim([0,max(time)])
  plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
  plt.ylabel(r'$M_\mathrm{gas}$', fontsize=15)

  fname = fdir + 'Gas_Mass.png'
  plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()

  plt.plot(time,h[0,:], color = 'k', linestyle = '--')  #mass in planet 1
  plt.plot(time,h[1,:], color = 'k', linestyle = '-.')  #mass in planet 1
  plt.plot(time,h[0,:]+h[1,:], color = 'k', linestyle = '-')  #mass in planet 1
  plt.xlim([0,max(time)])
  plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
  plt.ylabel(r'$h$', fontsize=15)

  fname = fdir + 'specific_orbital_angularmomentum.png'
  plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.clf()



def Energy(E,f):
  time = E[0]
  plt.plot(time, E[1], color = 'k')
  plt.plot(time, E[2], color = 'k', linestyle = '--') #dashed (potential energy)
  plt.plot(time, E[3], color = 'k', linestyle = ':')  #dotted (kinetic energy)
  plt.plot(time, E[4], color = 'k', linestyle = '-.') #dash-dotted (internal energy)
  plt.xlabel(r'$t\ \mathrm{[yr]}$', fontsize=15)
  plt.ylabel(r'$E$', fontsize=15)
  #plt.axis([0, max(time), 1.2*np.min(E), 1.2*np.max(E)])

  fname = f + 'Energy.ps'
  plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
  fname = f + 'Energy.png'
  plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
  plt.show()
  plt.clf()














