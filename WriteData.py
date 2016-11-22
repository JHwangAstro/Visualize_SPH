import struct
import math
import pickle
import subprocess
import os
import shutil


def Write_Input():
#_________________________________________________________________WRITES OPTIONS FILE

  f = open(findir + 'Input.txt', 'w')
  f.write('Change formatting in readdat.py in function init\n\n')
  f.write('#_________________________________________________________________INITIAL OPTIONS\n\n')
  if run > 2:
    f.write('\nMass1  = ' + str(m1[0]))
    f.write('\nMass2  = {0}' + str(m2[0]))
    f.write('\nCore1  = {0}' + str(cor1))
    f.write('\nCore2  = ' + str(cor2))
  else:
    f.write('\nMass   = ' + str(m1[0]))
    f.write('\nCore   = ' + str(cor1))

  if run == 3: f.write('\nSepI  = ' + str(rsep));

  f.write('\nntot   = ' + str(ntot))
  f.write('\nnnopt  = ' + str(nnopt))
  f.write('\nhmin   = ' + str(hmin))
  f.write('\nhmax   = ' + str(hmax))
  f.write('\nsep0   = ' + str(sep0))
  f.write('\ntf     = ' + str(tf))
  f.write('\ndtout  = ' + str(dtout))
  f.write('\nnout   = ' + str(nout))
  f.write('\nnit    = ' + str(nit))
  f.write('\ntime   = ' + str(time))
  f.write('\nnav    = ' + str(nav))
  f.write('\nalpha  = ' + str(alpha))
  f.write('\nbeta1  = ' + str(beta1))
  f.write('\ntskip  = ' + str(tskip))
  f.write('\nngr    = ' + str(ngr))
  f.write('\nnrelax = ' + str(nrelax))
  f.write('\ntrelax = ' + str(trelax))
  f.write('\ndt     = ' + str(dt))
  f.write('\nomega2 = ' + str(omega2))

  f.close()




def Origin_Shift(Origin,Rot,sdir,i,Header,Particles,affix):

  #_______________________________________________________________SETS FILENAME FOR SNAPSHOT
  fname = sdir + 'new_out000' + repr(i) + '.sph'
  if i >= 10:   fname = sdir + affix + 'out00' + repr(i) + '.sph'
  if i >= 100:  fname = sdir + affix + 'out0'  + repr(i) + '.sph'
  if i >= 1000: fname = sdir + affix + 'out'   + repr(i) + '.sph'

  #_______________________________________________________________WRITES DATA FILE(S)
  f = open(fname, 'wb')
  #Header = ReadData.Read_Header(f)
  dum = 124
  f.write(struct.pack('i',dum))                       #dummy variable
  f.write(struct.pack('i',Header.ntot))               #total number of particles
  f.write(struct.pack('i',Header.nnopt))              #optimal number of neighbors
  f.write(struct.pack('d',Header.hmin))               #smoothing length for black holes
  f.write(struct.pack('d',Header.hmax))               #smoothing length for black holes
  f.write(struct.pack('d',Header.sep0))               #initial separation of two stars
  f.write(struct.pack('d',Header.tf))                 #time the code stops in code units
  f.write(struct.pack('d',Header.dtout))              #time interval between writing outxxx.sph files
  f.write(struct.pack('i',Header.nout))               #xxx value of outputxxx.sph
  f.write(struct.pack('i',Header.nit))                #number of iterations completed
  f.write(struct.pack('d',Header.time))               #current time in code units
  f.write(struct.pack('i',Header.nav))                #chooses artificial viscosity scheme
  f.write(struct.pack('d',Header.alpha))              #artificial viscosity parameter
  f.write(struct.pack('d',Header.beta))               #artificial viscosity parameter
  f.write(struct.pack('d',Header.tskip))              #??
  f.write(struct.pack('i',Header.ngr))                #gravity boundary conditions
  f.write(struct.pack('i',Header.nrelax))             #relaxation flag
  f.write(struct.pack('d',Header.trelax))             #timescale for artificial drag force
  f.write(struct.pack('d',Header.dt))                 #current timestep
  f.write(struct.pack('d',Header.omega2))             #square of the orbital velocity
  f.write(struct.pack('i',dum))                       #dummy variable

  dum = 140
  for i in range(len(Particles)):
    x1 = (Particles[i].x-Origin[0])
    y1 = (Particles[i].y-Origin[1])
    z1 = (Particles[i].z-Origin[2])
    x2 = x1*math.cos(Rot[1])+0.+z1*math.sin(Rot[1])
    y2 = 0.                 +y1+0.
    z2 =-z1*math.sin(Rot[1])+0.+z1*math.cos(Rot[1])
    f.write(struct.pack('i',dum))                     #dummy variable
    f.write(struct.pack('d',x2))                    #x position
    f.write(struct.pack('d',y2))                    #y position
    f.write(struct.pack('d',z2))                    #z position
    f.write(struct.pack('d',Particles[i].m))          #mass of the particle
    f.write(struct.pack('d',Particles[i].h))          #smoothing length
    f.write(struct.pack('d',Particles[i].rho))        #density
    f.write(struct.pack('d',Particles[i].vx))         #x-velocity
    f.write(struct.pack('d',Particles[i].vy))         #y-velocity
    f.write(struct.pack('d',Particles[i].vz))         #z-velocity
    f.write(struct.pack('d',Particles[i].vxdot))      #x-acceleration
    f.write(struct.pack('d',Particles[i].vydot))      #y-acceleration
    f.write(struct.pack('d',Particles[i].vzdot))      #z-acceleration
    f.write(struct.pack('d',Particles[i].u))          #potential
    f.write(struct.pack('d',Particles[i].udot))       #change in potential
    f.write(struct.pack('d',Particles[i].grpot))      #gravitational potential
    f.write(struct.pack('d',Particles[i].mmw))        #mean molecular weight
    f.write(struct.pack('i',Particles[i].cc))         #flag for particle type
    f.write(struct.pack('d',Particles[i].divv))       #divergence of velocities
    f.write(struct.pack('i',dum))                     #dummy variable


  f.close()



#Writes a summary file
#in - fdir,P0,P1
#P0 - Initial planets
#P1 - Final Planets
def Write_Summary(fdir,P0,P1,Star,Mass,Energy,E,h,Mgas):
  f = open(fdir+'Summary.txt', 'w')
  f.write('Initial - Mass'+'\t'+'Semi-major axis'+'\t'+'Eccentricity'+'\t'+'Inclination'+'\t'+'Argument of Periastron'+'\t'+'Ascending Node'+'\t'+'True Anomaly'+'\n')
  for i in range(0,len(P0)):
    f.write(str(P0[i].M)+'\t'+str(P0[i].a*0.00465)+'\t'+str(P0[i].e)+'\t'+str(P0[i].i)+'\t'+str(P0[i].g)+'\t'+str(P0[i].n)+'\t'+str(P0[i].nu)+'\n')
  f.write('\n')
  f.write('Final - Mass'+'\t'+'Semi-major axis'+'\t'+'Eccentricity'+'\t'+'Inclination'+'\t'+'Argument of Periastron'+'\t'+'Ascending Node'+'\t'+'True Anomaly'+'\n')
  for i in range(0,len(P1)):
    if P1[i].M > 0: f.write(str(P1[i].M)+'\t'+str(P1[i].a*0.00465)+'\t'+str(P1[i].e)+'\t'+str(P1[i].i)+'\t'+str(P1[i].g)+'\t'+str(P1[i].n)+'\t'+str(P1[i].nu)+'\n')
  if len([j for j in range(len(P1)) if P1[j].M > 0]) > 1:
    f.write('Planet-Pair'+'\t'+'Dynamical Separation'+'\n')
    da = math.fabs(P1[0].a-P1[1].a)
    RH = math.pow((P1[0].M+P1[1].M)/(3.*Star.M),1./3.)*(P1[0].a+P1[1].a)*0.5  #mutual Hill radius
    f.write('12'+'\t'+str(da/RH))
    print 'Dynamical Separation: ', da/RH
    if da/RH > math.sqrt(3.)*2.: print 'STABLE!'

  N = len(Mass[0,:])-1
  f.write('\n\nChanges from Initial to Final')
  f.write('\ndM1'+'\t'+'dM2'+'\t'+'M_env'+'\t'+'M_eject'+'\t'+'dE1'+'\t'+'dE2'+'\t'+'dE'+'\t'+'dh1'+'\t'+'dh2'+'\t'+'dh')
  f.write('\n'+str(Mass[0,N]-Mass[0,0])+'\t'+str(Mass[1,N]-Mass[1,0])+'\t'+str(Mass[2,N])+'\t'+str(Mass[3,N])+'\t'+str(E[0,N]-E[0,0])+'\t'+str(E[1,N]-E[1,0])+'\t'+str(E[0,N]+E[1,N]-E[0,0]-E[1,0])+'\t'+str(h[0,N]-h[0,0])+'\t'+str(h[1,N]-h[1,0])+'\t'+str(h[0,N]+h[1,N]-h[0,0]-h[1,0]))
  f.write('\n\nChanges from Initial to Final in percentages')
  f.write('\ndM1'+'\t'+'dM2'+'\t'+'M_env'+'\t'+'M_eject'+'\t'+'dE'+'\t'+'dh')
  f.write('\n'+str((Mass[0,N]-Mass[0,0])/Mass[0,0])+'\t'+str((Mass[1,N]-Mass[1,0])/Mass[1,0])+'\t'+str(Mass[2,N]/(Mass[0,0]+Mass[1,0]))+'\t'+str(Mass[3,N]/(Mass[0,0]+Mass[1,0]))+'\t'+str((E[0,N]-E[0,0])/E[0,0])+'\t'+str((E[1,N]-E[1,0])/E[1,0])+'\t'+str((E[0,N]+E[1,N]-E[0,0]-E[1,0])/(E[0,0]+E[1,0]))+'\t'+str((h[0,N]-h[0,0])/h[0,0])+'\t'+str((h[1,N]-h[1,0])/h[1,0])+'\t'+str((h[0,N]+h[1,N]-h[0,0]-h[1,0])/(h[0,0]+h[1,0])))
  f.close()




def Write_Restart(sdir,t,E,Specific_E,Specific_h,Mgas,M,aei,aei2,P,P0,S,n):

  with open(sdir+'sph.restart_analysis', 'wb') as output:
    pickle.dump(t,    output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(E,    output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(Specific_E, output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(Specific_h, output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(Mgas, output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(M,    output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(aei,  output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(aei2, output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(P,    output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(P0,   output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(S,    output, pickle.HIGHEST_PROTOCOL)
    pickle.dump(n,    output, pickle.HIGHEST_PROTOCOL)





def Write_Mercury(dst,sdir,Planets):
  #dst is the location to make Mercury Run
  #sdir is the location of the base Mercury directory
  #provide all in mercury units (AU,degrees for last 3 angles)
  fdir = dst+'/Mercury'
  if not os.path.exists(fdir):
    shutil.copytree(sdir,fdir)
  os.chdir(fdir)

  mass = [Planets[0].M,Planets[1].M]
  density = [10000000.0,10000000.0]

  a = [Planets[0].a*0.00465,Planets[1].a*0.00465]
  e = [Planets[0].e,Planets[1].e]
  I = [Planets[0].i,Planets[1].i]
  g = [Planets[0].g,Planets[1].g]
  n = [Planets[0].n,Planets[1].n]
  M = [Planets[0].nu,Planets[1].nu]

  ang_x = [0.0,0.0]
  ang_y = [0.0,0.0]
  ang_z = [0.0,0.0]

  f = open(fdir+'/big.in', 'w')
  f.write(')O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n')
  f.write(')---------------------------------------------------------------------\n')
  f.write(' style (Cartesian, Asteroidal, Cometary) = Asteroidal\n')
  f.write(' epoch (in days) = 0.0\n')#)2451000.5\n')
  f.write(')---------------------------------------------------------------------\n')

  for j in range(0,len(a)):
    f.write('Kepler'+str(j)+'    m='+str(mass[j])+' r='+str(3)+' d='+str(density[j])+'\n')
    f.write(' '+str(a[j])+' '+str(e[j])+' '+str(I[j])+' '+str(g[j])+' '+str(n[j])+' '+str(M[j])+'\n')
    f.write(str(ang_x[j])+' '+str(ang_y[j])+' '+str(ang_z[j])+'\n')

  f.close()

  #________________________________________________RUN MERCURY

  #write pbs file
  f = open(fdir+'/jobscript.pbs', 'w')
  f.write('#interpret commands (except PBS directives) using the bash shell\n')
  f.write('#!/bin/bash\n')
  f.write('#Request one nodes using 1 cores per node\n')
  f.write('#PBS -l nodes=1:ppn=1\n')
  f.write('#PBS -l walltime=7:00:00:00\n')
  f.write('#PBS -N Mercury'+'\n')
  f.write('#PBS -j oe\n')
  f.write('#Change the current working directory to the directory from which the script was submitted\n')
  f.write('cd $PBS_O_WORKDIR\n')
  f.write('#These Commands are what would be needed to run this job natively\n')
  f.write('module load intel/2011.3\n')
  f.write('./mercury6')
  f.close()
  #submit pbs file
  #subprocess.call("msub -A b1011 -q astro jobscript.pbs", shell=True)












