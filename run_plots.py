
import math
import matplotlib.pyplot as plt
import numpy as np
import os


def plot_energy(snapshots, fdir, run_id):

    g   = 6.67390*math.pow(10,-8)
    ms = 1.98855*math.pow(10.,33.)
    rs  = 6.9599*math.pow(10,10)
    ts  = math.sqrt(math.pow(rs,3)/(g*ms))

    orbits_dir = fdir + 'Orbits'
    if not os.path.exists(orbits_dir): os.makedirs(orbits_dir)

    t = [snapshot.unitless_time for snapshot in snapshots]
    e_orbit = [snapshot.orbital_energy*ms*rs*rs/ts/ts for snapshot in snapshots]
    #print e_orbit

    plt.plot(t, e_orbit, color = 'k', linestyle = '-')  #semi-major axes

    plt.xlim([0,max(t)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$E_\mathrm{orbit}$', fontsize=15)

    fname = orbits_dir + '/Run_' + run_id + '_e_orbit.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = orbits_dir + '/Run_' + run_id + '_e_orbit.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()


