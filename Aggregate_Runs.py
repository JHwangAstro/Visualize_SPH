#!/usr/bin/env python

# Input:  *.sph files from StarCrash
# Output: visualization of the data of the system
# Aggregate_Runs.py (SPH Visualization)

# example command to run analysis
# python Aggregate_Runs.py -d 'Collisions/V10/Grid/4M_4M/85Mc/Real_Runs/' -ms 1.0 -c None -rd '4M4M_85Mc' -ic '4M4M_85Mc_collision_parameters.txt' -sf 'save.4M4M_85Mc_summary' -sn '$q=1; m_\mathrm{c}=0.85$'
# python Aggregate_Runs.py -d 'Collisions/V10/Grid/4M_4M/95Mc/' -ms 1.0 -c None -rd '4M4M_95Mc' -ic '4M4M_95Mc_collision_parameters.txt' -sf 'save.4M4M_95Mc_summary' -sn '$q=1; m_\mathrm{c}=0.95$'
# python Aggregate_Runs.py -d 'Collisions/V10/Grid/4M_12M/85Mc/Real_Runs/' -ms 1.0 -c None -rd '4M12M_85Mc' -ic '4M12M_85Mc_collision_parameters.txt' -sf 'save.4M12M_85Mc_summary' -sn '$q=1/3; m_\mathrm{c}=0.85$'
# python Aggregate_Runs.py -d 'Collisions/V10/Grid/4M_12M/95Mc/' -ms 1.0 -c None -rd '4M12M_95Mc' -ic '4M12M_95Mc_collision_parameters.txt' -sf 'save.4M12M_95Mc_summary' -sn 'q=1/3; m_\mathrm{c}=0.95$'
# python Aggregate_Runs.py -d 'V10/Kepler11/' -i ['054', '068'] -ms 0.95 -c None -rd 'Kepler-11'

import argparse                     #Reading in arguments
import os                           #For creating folders
import sys
from summary import Summary


#Main function
def main():

    args = getArgs()
    summary = Summary(args)




#Creates the input arguments for the analysis
def getArgs():

    # Default values for input arguments
    base    = ""                                                    # sets base directory
    rund    = "Collisions/V10/"                                     # sets run directory
    relax_dir = "RelaxationRuns/V10/"
    relax_run = ["Kepler36b_0.111AU/", "Kepler36c_0.132AU/"]
    sum_dir = "summary/"                                            # directory to write summary files
    prefix  = "Run"                                                 # prefix for runs
    ids     = []                                                    # affixes for runs; if empty will look for any folder with prefix
    num_p   = 2                                                     # number of planets; assume for now always 2
    num_p2  = 0                                                     # number of other planets
    ndir    = "mercury_xyz/"
    sum_name = os.getcwd() + "/save.kepler36_summary"
    rname   = "sph.restart_analysis_new"                            # restart file name
    rname_old = "sph.restart_analysis"                              # old restart file name
    fname_ic = "kepler36_collision_parameters.txt"                  # filename of text file with initial conditions
    star_mass = 1.071                                               # mass of Kepler-36
    write   = 0
    compare = "Kepler-36"
    structure = "Kepler-36"
    rerun = 0
    set_name = 'Kepler-36 Progenitor'

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nplanets",   dest="num_p",    default=num_p,    type = int, help="Number of Planets")
    parser.add_argument("-m", "--nplanets2",  dest="num_p2",   default=num_p2,   type = int, help="Number of non-SPH Planets")
    parser.add_argument("-b", "--projectdir", dest="bdir",     default=base,     help="Project Directory")
    parser.add_argument("-d", "--rundir",     dest="rdir",     default=rund,     help="Run Directory")
    parser.add_argument("-rx","--relaxdir",   dest="relax_dir",default=relax_dir,help="Relax Directory")
    parser.add_argument("-re","--relaxruns",  dest="relax_run",default=relax_run,help="Relax Runs")
    parser.add_argument("-d2","--nbody_dir",  dest="ndir",     default=ndir,     help="N-Body Run Directory")
    parser.add_argument("-rf","--restartf",   dest="rname",    default=rname,    help="Restart file name")
    parser.add_argument("-sf","--runsum_f",   dest="sum_name", default=sum_name, help="Run summary restart file name")
    parser.add_argument("-of","--restart_old",dest="rname_old",default=rname_old,help="Old restart file name")
    parser.add_argument("-c", "--compare",    dest="compare",  default=compare,  help="Compare to a system")
    parser.add_argument("-i", "--run_ids",    dest="ids",      default=ids,      help="run ids to aggregate", nargs='+')
    parser.add_argument("-ic","--fname_ic",   dest="fname_ic", default=fname_ic, help="name of initial conditions text file, relative to rdir")
    parser.add_argument("-p", "--run_prefix", dest="prefix",   default=prefix,   help="prefix of runs")
    parser.add_argument("-s", "--sum_dir",    dest="sum_dir",  default=sum_dir,  help="directory to store summary files")
    parser.add_argument("-ms","--star_mass",  dest="star_mass",default=star_mass,type = float, help="Mass of host star")
    parser.add_argument("-rd","--structure",  dest="structure",default=structure,help="initial structure of planets")
    parser.add_argument("-rr","--rerun",      dest="rerun",    default=rerun,    help="Rewrite dump file", action="store_true")
    parser.add_argument("-sn","--set_name",   dest="set_name", default=set_name, help="Set name")

    args = parser.parse_args()

    # If ids is an empty list, then automatically find assign each folder with the correct prefix to the list
    if not len(args.ids): args.ids = getRunFolders(args.bdir + args.rdir, args.prefix)
    # if radii is a string, look up values
    if not isinstance(args.structure, list):
        args.structure = get_structure(args.structure)

    # Set default values
    args.sdir = [args.bdir + args.rdir + args.prefix + affix + '/' for affix in args.ids]     # sets folders for snapshots
    args.fname_ic = args.bdir + args.rdir + args.fname_ic                               # sets folder for initial conditions
    args.idir = [directory + 'IMG/' for directory in args.sdir]                         # sets folder for images
    args.fdir = [directory + 'Final/' for directory in args.sdir]                       # sets folder for final plots
    args.ndir = [directory + args.ndir for directory in args.sdir]                      # directory for base mercury folder
    args.sum_dir = args.bdir + args.rdir + args.sum_dir                                 # directory to write summary files
    args.relax_run = [args.bdir + args.relax_dir + relax_runs + "sph.profile" for relax_runs in args.relax_run]
    if args.compare == 'None': args.compare = None

    #TODO - create aggregate folders
    createFolders(args)

    return args




# Finds and returns the number of appropriate sph files in the source directory
def getRunFolders(base_directory, prefix):

    return [x[len(prefix):] for x in next(os.walk(base_directory))[1] \
            if x[0:len(prefix)] == prefix and os.path.isfile(base_directory + x + '/restartrad.sph')]




def get_structure(name):

    g = 6.6738480*10.**-8.
    me = 5.9736*10**27.
    rs = 6.957*10.**10


    if name == 'Kepler-36':
        prad1 = 1687802440.83
        prad2 = 2144828233.8
        crad1 = 0.55728878053
        crad2 = 0.493310387376
        crad1_sph = 0.0143316888978*rs/prad1
        crad2_sph = 0.0161931906978*rs/prad2
        mass1 = 2.78790893168*10.**28.
        mass2 = 4.70014381397*10.**28.
        #mass1 = 4.45*me
        #mass2 = 8.08*me
        e_b = g*(mass1**2./prad1 + mass2**2./prad2)
        return {'r': [prad1, prad2], 'c': [crad1, crad2], 'm': [mass1, mass2], 'eb': e_b, 'c_sph': [crad1_sph, crad2_sph]}


    if name == '4M4M_85Mc':
        prad1 = 2680121082.15
        prad2 = 2523307175.48
        crad1 = 876268634.896/prad1
        crad2 = 876268634.867/prad2
        crad1_sph = 0.0138068224768*rs/prad1
        crad2_sph = 0.0138001562654*rs/prad2
        mass1 = 2.39043982304e+28
        mass2 = 2.39043982304e+28
        e_b = g*(mass1**2./prad1 + mass2**2./prad2)


    if name == '4M4M_95Mc':
        prad1 = 1808451511.47
        prad2 = 1734460278.95
        crad1 = 902738238.59/prad1
        crad2 = 902738238.547/prad2
        crad1_sph = 0.0138202931241*rs/prad1
        crad2_sph = 0.0137935141783*rs/prad2
        mass1 = 2.39043982304e+28
        mass2 = 2.39043982304e+28
        e_b = g*(mass1**2./prad1 + mass2**2./prad2)


    if name == '4M12M_85Mc':
        prad1 = 2680121082.15
        prad2 = 2554659954.54
        crad1 = 876268634.896/prad1
        crad2 = 1149568184.6/prad2
        crad1_sph = 0.0138068224562*rs/prad1
        crad2_sph = 0.0177501731554*rs/prad2
        mass1 = 2.39043982304e+28
        mass2 = 7.17131946913e+28
        e_b = g*(mass1**2./prad1 + mass2**2./prad2)


    if name == '4M12M_95Mc':
        prad1 = 1808451511.47
        prad2 = 1872208576.5
        crad1 = 902738238.59/prad1
        crad2 = 1178013833.44/prad2
        crad1_sph = 0.0138202905564*rs/prad1
        crad2_sph = 0.0178185412253*rs/prad2
        mass1 = 2.39043982304e+28
        mass2 = 7.17131946913e+28
        e_b = g*(mass1**2./prad1 + mass2**2./prad2)


    return {'r': [prad1, prad2], 'c': [crad1, crad2], 'm': [mass1, mass2], 'eb': e_b, 'c_sph': [crad1_sph, crad2_sph]}

    return None




# Create directories as needed; return an error if source directory does not exist
def createFolders(args):

    if not os.path.exists(args.bdir + args.rdir):
        print 'source directory does not exist'
        sys.exit
    if not os.path.exists(args.sum_dir): os.makedirs(args.sum_dir)




if __name__ == "__main__":

    main()








