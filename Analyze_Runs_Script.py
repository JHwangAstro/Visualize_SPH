
#!/usr/bin/env python

# Input:  *.sph files from StarCrash
# Output: visualization of the data of the system
# Aggregate_Runs.py (SPH Visualization)

# example command to run analysis
# python Analyze_Runs_Script.py -d 'Collisions/V10/Kepler11/' -i ['054', '068']
# python Analyze_Runs_Script.py -d 'Collisions/V10/Grid/4M_4M/85Mc/' -d3 '/jhwang/Mercury/Runs_Grid/post_SPH'

import argparse
import os
import sys
import subprocess


#Main function
def main():

    args = getArgs()

    f = './'
    submit_jobs(f,  args)




#Creates the input arguments for the analysis
def getArgs():

    # Default values for input arguments
    base    = ""                                                    # sets base directory
    rund    = "Collisions/V10/"                                     # sets run directory
    prefix  = "Run"                                                 # prefix for runs
    numss   = -1                                                    #number of snapshots
    nstart  = 0                                                     #starting snapshot
    restart = 0                                                     #Read restart file? (1=yes)
    ids     = []                                                    # affixes for runs; if empty will look for any folder with prefix
    num_p   = 2                                                     # number of planets; assume for now always 2
    num_p2  = 0                                                     # number of other planets
    mdir    = '/Mercury/Runs_Kepler36/post_SPH'                     #directory for base mercury folder
    ndir    = 'mercury_xyz/'
    rname   = 'sph.restart_analysis'                                # restart file name
    write   = 0
    compare = None

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--num_p",    dest="num_p",    default=num_p,     type = int, help="Number of Planets")
    parser.add_argument("-m", "--num_p2",   dest="num_p2",   default=num_p2,    type = int, help="Number of non-SPH Planets")
    parser.add_argument("-s", "--numss",    dest="numss",    default=numss,     type = int, help="Number of Snapshots")
    parser.add_argument("-k", "--nskip",    dest="nskip",    default=1,         type = int, help="Snapshot Interval")
    parser.add_argument("-z", "--nstart",   dest="nstart",   default=nstart,    type = int, help="Starting Snapshot")
    parser.add_argument("-w", "--rewrite",  dest="rewrite",  default=write,     help="Rewrite dump files", action="store_true")
    parser.add_argument("-r", "--restart",  dest="restart",  default=restart,   help="Restart", action="store_true")
    parser.add_argument("-b", "--bdir",     dest="bdir",     default=base,      help="Project Directory")
    parser.add_argument("-d", "--rdir",     dest="rdir",     default=rund,      help="Run Directory")
    parser.add_argument("-d2","--ndir",     dest="ndir",     default=ndir,      help="N-Body Run Directory")
    parser.add_argument("-d3","--mdir",     dest="mdir",     default=mdir,      help="N-Body Source Directory")
    parser.add_argument("-rf","--rname",    dest="rname",    default=rname,     help="Restart file name")
    parser.add_argument("-c", "--compare",  dest="compare",  default=compare,   help="Compare to a system")

    parser.add_argument("-i", "--ids",      dest="ids",      default=ids,       help="run ids to aggregate", nargs='+')
    parser.add_argument("-p", "--prefix",   dest="prefix",   default=prefix,    help="prefix of runs")

    args = parser.parse_args()

    #If ids is an empty list, then automatically find assign each folder with the correct prefix to the list
    if not len(args.ids): args.ids = getRunFolders(args.bdir + args.rdir, args.prefix)

    # Set default values
    args.sdir = [args.rdir + args.prefix + affix + '/' for affix in args.ids]     # sets folders for snapshots

    return args




# Finds and returns the number of appropriate sph files in the source directory
def getRunFolders(base_directory, prefix):

    # get all run folders in base directory
    run_dirs = [x[len(prefix):] for x in next(os.walk(base_directory))[1] if x[0:len(prefix)] == prefix]

    # return run folders that have a started run (existence of out0000.sph)
    return [x for x in run_dirs if os.path.isfile(base_directory + prefix + x + '/out0000.sph')]




def write_pbs_script(f, args, sdir):

    ignored_arguments = ['ids', 'sdir', 'prefix', 'rdir']
    string_args = ['rdir', 'rname', 'compare', 'ndir', 'mdir', 'bdir']

    arg_string = "--rdir '" + sdir + "'"
    for arg in vars(args):
        if str(arg) not in ignored_arguments and getattr(args, arg):
            if str(arg) in string_args:
                arg_string = arg_string + ' --' + str(arg) + " '" + str(getattr(args, arg)) + "'"
            else:
                arg_string = arg_string + ' --' + str(arg) + ' ' + str(getattr(args, arg))

    command = 'python Run_GasPlanet.py ' + arg_string
    print command


    # Write pbs file
    f = open(f + '/jobscript.pbs', 'w')
    f.write('#interpret commands (except PBS directives) using the bash shell\n')
    f.write('#!/bin/bash\n')
    f.write('#Request one nodes using 1 cores per node\n')
    f.write('#PBS -l nodes=1:ppn=1\n')
    f.write('#PBS -l walltime=4:00:00:00\n')
    f.write('#PBS -N Analyze_SPH\n')
    f.write('#PBS -j oe\n')
    f.write('#Change the current working directory to the directory from which the script was submitted\n')
    f.write('cd $PBS_O_WORKDIR\n')
    f.write('#These Commands are what would be needed to run this job natively\n')
    f.write('#module load intel/2011.3\n')

    f.write(command)
    f.close()




def submit_jobs(f, args):


    for sdir in args.sdir:

        write_pbs_script(f, args, sdir)
        subprocess.call("msub " + f + "jobscript.pbs", shell = True)




if __name__ == "__main__":
    main()







