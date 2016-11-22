#!/usr/bin/env python

# Input:  *.sph files from StarCrash
# Output: visualization of the data of the system

# example command
# python Run_GasPlanet.py -n 1 -r 'RelaxationRuns/Kepler36b_0.111AU/'
# python Run_GasPlanet.py -n 2 -r 'Collisions/Kepler36/Run124/'

#_________________________________________________________________PARAMETERS

import Analyze                      #Loops through each snapshot, reads and plots the data
from time import clock              #Calculate time spent
import os                           #For creating folders
import argparse                     #Reading in arguments
import sys

#Default values for input arguments
base  = "/projects/"                   #sets base directory
rund  = "Collisions/Run17/"                             #sets run directory
numss = -1#201                                              #number of snapshots
num_p = 2                                                       #number of planets
nstart= 0
restart=0                                                       #Read restart file? (1=yes)
num_p2 = 0

parser = argparse.ArgumentParser()
parser.add_argument("-n","--nplanets"   ,dest="num_p" ,default=num_p ,help="Number of Planets"  ,type=int)
parser.add_argument("-m","--nplanets2"  ,dest="num_p2",default=num_p2,help="Number of non-SPH Planets"  ,type=int)
parser.add_argument("-s","--nsnapshots" ,dest="numss" ,default=numss ,help="Number of Snapshots",type=int)
parser.add_argument("-k","--numskip"    ,dest="nskip" ,default=1     ,help="Snapshot Interval"  ,type=int)
parser.add_argument("-z","--startnum"   ,dest="nstart",default=nstart,help="Starting Snapshot",  type=int)
parser.add_argument("-b","--projectdir" ,dest="bdir"  ,default=base  ,help="Project Directory")
parser.add_argument("-r","--rundir"     ,dest="rdir"  ,default=rund  ,help="Run Directory")

args = parser.parse_args()

sdir  = args.bdir + args.rdir# + "SS/"                          #sets folder for snapshots
idir  = args.bdir + args.rdir + 'IMG/'                          #sets folder for images
fdir  = args.bdir + args.rdir + 'Final/'                        #sets folder for final plots
mdir  = '/projects/b1011/jhwang/Mercury/Runs_Kepler36/base'     #directory for base mercury folder

#If numss is default, then automatically find number of snapshots
if args.numss < 0:
  args.numss = 0
  fname = sdir + 'out000' + repr(args.numss) + '.sph'
  if args.numss >= 10:   fname = sdir + 'out00' + repr(args.numss) + '.sph'
  if args.numss >= 100:  fname = sdir + 'out0'  + repr(args.numss) + '.sph'
  if args.numss >= 1000: fname = sdir + 'out'   + repr(args.numss) + '.sph'
  while os.path.exists(fname):
    args.numss += args.nskip
    fname = sdir + 'out000' + repr(args.numss) + '.sph'
    if args.numss >= 10:   fname = sdir + 'out00' + repr(args.numss) + '.sph'
    if args.numss >= 100:  fname = sdir + 'out0'  + repr(args.numss) + '.sph'
    if args.numss >= 1000: fname = sdir + 'out'   + repr(args.numss) + '.sph'
  print 'there are ', args.numss, ' snapshots in this run'

if not os.path.exists(sdir):
  print 'source directory does not exist'
  sys.exit
if not os.path.exists(idir): os.makedirs(idir)
if not os.path.exists(fdir): os.makedirs(fdir)

if args.num_p == 2: Analyze.Planet_Planet_Collision(sdir,idir,fdir,mdir,nstart,args.numss/args.nskip,args.nskip,restart,args.num_p2)
if args.num_p == 1: Analyze.Planet_Profile(sdir,idir,fdir,args.numss)

















