#!/usr/bin/env python

# Input:  *.sph files from StarCrash
# Output: visualization of the data of the system
# Run_GasPlanet.py (SPH Visualization)

# example command to run analysis on relaxation calculation
# python Run_GasPlanet.py -n 1 -d 'RelaxationRuns/Kepler36b_0.111AU/'
# example command to run analysis on collision calculation
# python Run_GasPlanet.py -n 2 -d 'Collisions/Kepler36a/Run049/'

import argparse                     #Reading in arguments
import os                           #For creating folders
import sys
from run import Run


# if caled from another program
def called_from_external(args):

    run = Run(args)




# Main function
def main():

    args = getArgs()
    run = Run(args)




# Creates the input arguments for the analysis
def getArgs():

    # Default values for input arguments
    base    = "/jhwang/StarSmasher/"                 #sets base directory
    rund    = "Collisions/Kepler36a/Run049/"                        #sets run directory
    numss   = -1                                                    #number of snapshots
    num_p   = 2                                                     #number of planets
    nstart  = 0                                                     #starting snapshot
    restart = 0                                                     #Read restart file? (1=yes)
    num_p2  = 0
    mdir    = '/jhwang/Mercury/Runs_Kepler36/post_SPH'   #directory for base mercury folder
    #ndir    = 'mercury/'
    ndir    = 'mercury_xyz/'
    rname   = 'sph.restart_analysis'                                #restart file name
    write   = 0
    compare = None

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--num_p",    dest="num_p",    default=num_p,  type = int, help="Number of Planets")
    parser.add_argument("-m", "--num_p2",   dest="num_p2",   default=num_p2, type = int, help="Number of non-SPH Planets")
    parser.add_argument("-s", "--numss",    dest="numss",    default=numss,  type = int, help="Number of Snapshots")
    parser.add_argument("-k", "--nskip",    dest="nskip",    default=1,      type = int, help="Snapshot Interval")
    parser.add_argument("-z", "--nstart",   dest="nstart",   default=nstart, type = int, help="Starting Snapshot")
    parser.add_argument("-w", "--rewrite",  dest="rewrite",  default=write,  help="Rewrite dump files", action="store_true")
    parser.add_argument("-r", "--restart",  dest="restart",  default=restart,help="Restart", action="store_true")
    parser.add_argument("-b", "--bdir",     dest="bdir",     default=base,   help="Project Directory")
    parser.add_argument("-d", "--rdir",     dest="rdir",     default=rund,   help="Run Directory")
    parser.add_argument("-d2","--ndir",     dest="ndir",     default=ndir,   help="N-Body Run Directory")
    parser.add_argument("-d3","--mdir",     dest="mdir",     default=mdir,   help="N-Body Source Directory")
    parser.add_argument("-rf","--rname",    dest="rname",    default=rname,  help="Restart file name")
    parser.add_argument("-c", "--compare",  dest="compare",  default=compare,help="Compare to a system")

    args = parser.parse_args()

    # Set default values
    args.sdir = args.bdir + args.rdir# + "SS/"                          #sets folder for snapshots
    args.idir = args.bdir + args.rdir + 'IMG/'                          #sets folder for images
    args.fdir = args.bdir + args.rdir + 'Final/'                        #sets folder for final plots
    args.ndir = args.bdir + args.rdir + args.ndir   #directory for base mercury folder

    #If numss is default (-1), then automatically find number of snapshots
    if args.numss < 0: args.numss = findSnapshotCount(args.sdir,args.nskip)
    createFolders(args)

    return args



# Finds and returns the number of appropriate sph files in the source directory
def findSnapshotCount(f, nskip):

    numss = 0
    fname = f + 'out' + format(numss, "04") + '.sph'
    while os.path.exists(fname):
        numss += nskip
        fname = f + 'out' + format(numss, "04") + '.sph'
    print 'there are ', numss, ' snapshots in this run'

    return numss




# Create directories as needed; return an error if source directory does not exist
def createFolders(args):

    if not os.path.exists(args.sdir):
        print 'source directory does not exist'
        sys.exit
    if not os.path.exists(args.idir): os.makedirs(args.idir)
    if not os.path.exists(args.fdir): os.makedirs(args.fdir)




if __name__ == "__main__":
    main()








