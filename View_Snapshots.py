#!/usr/bin/env python

# Input:  *.sph files from StarCrash
# Output: visualization of the data of the system
# Run_GasPlanet.py (SPH Visualization)

# example command to run analysis on relaxation calculation
# python Run_GasPlanet.py -n 1 -d 'RelaxationRuns/Kepler36b_0.111AU/'
# example command to run analysis on collision calculation
# python Run_GasPlanet.py -n 2 -d 'Collisions/Kepler36a/Run049/'

# This script plots the selected snapshots

import os                           #For creating folders
import argparse                     #Reading in arguments
import sys
from snapshot import Snapshot
import numpy as np


#Main function
def main():

    args = getArgs()
    if args.numss < 0: args.numss = args.nstart+1
    snapshot_keys = np.arange(args.nstart, args.numss, args.nskip)
    snapshot = {}
    for k in snapshot_keys:
        snapshot[k] = Snapshot(k)
        snapshot[k].initialize(args, snapshot[k])



#Creates the input arguments for the analysis
def getArgs():

    # Default values for input arguments
    base    = "/jhwang/StarSmasher/"                 #sets base directory
    rund    = "Collisions/Kepler36a/Run049/"                        #sets run directory
    numss   = -1                                                    #number of snapshots
    num_p   = 2                                                     #number of planets
    nstart  = 0                                                     #starting snapshot
    restart = 0                                                     #Read restart file? (1=yes)
    num_p2  = 0
    defaultDir = base
    mdir    = '/jhwang/Mercury/Runs_Kepler36/base'   #directory for base mercury folder
    rname   = 'sph.restart_analysis'                                #restart file name

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nplanets",  dest="num_p",    default=num_p,  type = int, help="Number of Planets")
    parser.add_argument("-m", "--nplanets2", dest="num_p2",   default=num_p2, type = int, help="Number of non-SPH Planets")
    parser.add_argument("-s", "--nsnapshots",dest="numss",    default=numss,  type = int, help="Number of Snapshots")
    parser.add_argument("-k", "--numskip",   dest="nskip",    default=1,      type = int, help="Snapshot Interval")
    parser.add_argument("-z", "--startnum",  dest="nstart",   default=nstart, type = int, help="Starting Snapshot")
    parser.add_argument("-r", "--restart",   dest="restart",  default=rund,   help="Restart", action="store_true")
    parser.add_argument("-b", "--projectdir",dest="bdir",     default=base,   help="Project Directory")
    parser.add_argument("-d", "--rundir",    dest="rdir",     default=rund,   help="Run Directory")

    args = parser.parse_args()

    # Set default values
    args.sdir = args.bdir + args.rdir# + "SS/"                          #sets folder for snapshots
    args.idir = args.bdir + args.rdir + 'IMG/'                          #sets folder for images
    args.fdir = args.bdir + args.rdir + 'Final/'                        #sets folder for final plots

    return args



# Create directories as needed; return an error if source directory does not exist
def createFolders(args):

    if not os.path.exists(args.sdir):
        print 'source directory does not exist'
        sys.exit
    if not os.path.exists(args.idir): os.makedirs(args.idir)




if __name__ == "__main__":
    main()








