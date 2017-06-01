
import math
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

def energy_plots(fdir, energy, P0):

    ts  = 5.05000318388e-05                                     # Convert from code units to years
    keys = ['eg', 'ek', 'ei', 'et']

    # plots energy as a function of time
    t = np.array([k*ts/P0 for k in sorted(energy)]) # convert from code units to year to P0
    eg = np.array([energy[k]['eg'] for k in sorted(energy)])
    ek = np.array([energy[k]['ek'] for k in sorted(energy)])
    ei = np.array([energy[k]['ei'] for k in sorted(energy)])
    et = np.array([energy[k]['et'] for k in sorted(energy)])
    et0 = math.fabs(et[0])

    plt.plot(t, et/et0, color = 'k', linestyle = '-')
    plt.plot(t, eg/et0, color = 'k', linestyle = '--')
    plt.plot(t, ek/et0, color = 'k', linestyle = ':')
    plt.plot(t, ei/et0, color = 'k', linestyle = '-.')

    plt.xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize = 15)
    plt.ylabel(r'$e/e_\mathrm{i}$', fontsize = 15)

    plt.xlim([0, max(t)])

    fname = fdir + 'energy.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'energy.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()




def combined_orbits(fdir, snapshots, ref_system, pre_sph, post_sph):

    c = ['b','r']

    skip1 = 1
    skip = 1000

    # Create a list of used planets attributes
    t = np.array([snapshot.unitless_time for snapshot in snapshots])
    a = [np.array([snapshot.planets[i].aei['a'] for snapshot in snapshots]) for i in range(2)]
    e = [np.array([snapshot.planets[i].aei['e'] for snapshot in snapshots]) for i in range(2)]

    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey = True)

    for i in range(2):

        if pre_sph.exists:
            ax1.plot(pre_sph.t[i][::skip1], pre_sph.a[i][::skip1], color = c[i], linestyle = '-')
            ax1.plot(pre_sph.t[i][::skip1], pre_sph.a[i][::skip1]*(1.-pre_sph.e[i][::skip1]), color = c[i], linestyle = ':')
            ax1.plot(pre_sph.t[i][::skip1], pre_sph.a[i][::skip1]*(1.+pre_sph.e[i][::skip1]), color = c[i], linestyle = ':')

        am = np.ma.masked_where(a[i] == 0., a[i])
        em = np.ma.masked_where(a[i] == 0., e[i])
        ax2.plot(t, am, color = c[i], linestyle = '-')  #semi-major axes
        ax2.plot(t, (1.-em)*am, color = c[i], linestyle = ':')  #periapsis
        ax2.plot(t, (1.+em)*am, color = c[i], linestyle = ':')  #apoapsis

        if post_sph.exists:
            ax3.plot(post_sph.t[i][::skip]/10.**6., post_sph.a[i][::skip], color = c[i], linestyle = '-')
            ax3.plot(post_sph.t[i][::skip]/10.**6., post_sph.a[i][::skip]*(1.-post_sph.e[i][::skip]), color = c[i], linestyle = ':')
            ax3.plot(post_sph.t[i][::skip]/10.**6., post_sph.a[i][::skip]*(1.+post_sph.e[i][::skip]), color = c[i], linestyle = ':')

    if pre_sph.exists:
        ax1.set_xlim([0,max(pre_sph.t[0])])
        ax1.set_xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize = 15)
        ax1.set_ylabel(r'$a\ [\mathrm{AU}]$', fontsize = 15)
        ax1.xaxis.major.locator.set_params(nbins=5)

    ax2.set_xlim([0,max(t)])
    ax2.set_xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize=15)
    ax2.xaxis.major.locator.set_params(nbins=5)

    if post_sph.exists:
        ax3.set_xlim([0,max(post_sph.t[0]/10.**6.)])
        ax3.set_xlabel(r'$t/(10^6\mathrm{P}_\mathrm{inner})$')
        ax3.xaxis.major.locator.set_params(nbins=5)

    fname = fdir + 'orbits_combined.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_combined.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    # Create orbit plots
    for i in range(2):
        if ref_system.exists:
            if pre_sph.exists: ax1.plot([min(pre_sph.t[0]), max(pre_sph.t[0])], [ref_system.a[i], ref_system.a[i]], linestyle = '--', color = c[i])
            ax2.plot([min(t), max(t)], [ref_system.a[i], ref_system.a[i]], linestyle = '--', color = c[i])
            if post_sph.exists: ax3.plot([min(post_sph.t[0]), max(post_sph.t[0])/10.**6.], [ref_system.a[i], ref_system.a[i]], linestyle = '--', color = c[i])

    fname = fdir + 'orbits_combined_compare.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_combined_compare.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()


    # Create a list of used planets attributes

    f, (ax1, ax2) = plt.subplots(1, 2, sharey = True, gridspec_kw = {'width_ratios':[1,2]})

    for i in range(2):

        if pre_sph.exists:
            ax1.plot(pre_sph.t[i][::skip1], pre_sph.a[i][::skip1], color = c[i], linestyle = '-')
            ax1.plot(pre_sph.t[i][::skip1], pre_sph.a[i][::skip1]*(1.-pre_sph.e[i][::skip1]), color = c[i], linestyle = ':')
            ax1.plot(pre_sph.t[i][::skip1], pre_sph.a[i][::skip1]*(1.+pre_sph.e[i][::skip1]), color = c[i], linestyle = ':')

        am = np.ma.masked_where(a[i] == 0., a[i])
        em = np.ma.masked_where(a[i] == 0., e[i])
        ax2.plot(t, am, color = c[i], linestyle = '-')  #semi-major axes
        ax2.plot(t, (1.-em)*am, color = c[i], linestyle = ':')  #periapsis
        ax2.plot(t, (1.+em)*am, color = c[i], linestyle = ':')  #apoapsis

    if pre_sph.exists:
        ax1.set_xlim([0,max(pre_sph.t[0])])
        ax1.set_xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize = 15)
        ax1.set_ylabel(r'$a\ [\mathrm{AU}]$', fontsize = 15)
        ax1.xaxis.major.locator.set_params(nbins=5)

    ax2.set_xlim([0,max(t)])
    ax2.set_xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize=15)
    ax2.xaxis.major.locator.set_params(nbins=5)

    fname = fdir + 'orbits_combined.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_combined.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    # Create orbit plots
    for i in range(2):
        if ref_system.exists:
            if pre_sph.exists: ax1.plot([min(pre_sph.t[0]), max(pre_sph.t[0])], [ref_system.a[i], ref_system.a[i]], linestyle = '--', color = c[i])
            ax2.plot([min(t), max(t)], [ref_system.a[i], ref_system.a[i]], linestyle = '--', color = c[i])

    fname = fdir + 'orbits_combined_compare_nopost.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_combined_compare_nopost.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()


    f, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey = True)

    for i in range(2):

        if pre_sph.exists:
            ax1.plot(pre_sph.t[i][::skip1], pre_sph.e[i][::skip1], color = c[i], linestyle = '-')

        em = np.ma.masked_where(a[i] == 0., e[i])
        ax2.plot(t, em, color = c[i], linestyle = '-')  #semi-major axes

        if post_sph.exists:
            ax3.plot(post_sph.t[i][::skip]/10.**6., post_sph.e[i][::skip], color = c[i], linestyle = '-')

    ax1.set_ylabel(r'$e$', fontsize = 15)
    if pre_sph.exists:
        ax1.set_xlim([0,max(pre_sph.t[0])])
        ax1.set_xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize = 15)
        ax1.xaxis.major.locator.set_params(nbins=5)

    ax2.set_xlim([0,max(t)])
    ax2.set_xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize=15)
    ax2.xaxis.major.locator.set_params(nbins=5)

    if post_sph.exists:
        ax3.set_xlim([0,max(post_sph.t[0]/10.**6.)])
        ax3.set_xlabel(r'$t/(10^6\mathrm{P}_\mathrm{inner})$')
        ax3.xaxis.major.locator.set_params(nbins=5)

    fname = fdir + 'orbits_combined_e.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_combined_e.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    # Create orbit plots
    if ref_system.exists:
        for i in range(2):
            if pre_sph.exists: ax1.plot([min(pre_sph.t[0]), max(pre_sph.t[0])], [ref_system.e[i], ref_system.e[i]], linestyle = '--', color = c[i])
            ax2.plot([min(t), max(t)], [ref_system.e[i], ref_system.e[i]], linestyle = '--', color = c[i])
            if post_sph.exists: ax3.plot([min(post_sph.t[0]), max(post_sph.t[0])/10.**6.], [ref_system.e[i], ref_system.e[i]], linestyle = '--', color = c[i])

        fname = fdir + 'orbits_combined_compare_e.png'
        plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
        fname = fdir + 'orbits_combined_compare_e.ps'
        plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()




def orbit_plots(fdir, snapshots, ref_system):

    c = ['b','r']
    ms = 1.98855*math.pow(10.,33.)                          #Solar mass in g
    me = 5.972*math.pow(10.,27.)                            #Earth mass in g

    # Create an array of used snapshot attributes
    t = np.array([snapshot.unitless_time for snapshot in snapshots])
    #t = np.array([snapshot.t for snapshot in snapshots])
    m_bound_gas = np.array([snapshot.bound_mass for snapshot in snapshots])
    m_unbound_gas = np.array([snapshot.ejected_mass for snapshot in snapshots])
    planet_energy = np.array([snapshot.planet_energy for snapshot in snapshots])

    # Create a list of used planets attributes
    a = [np.array([snapshot.planets[i].aei['a'] for snapshot in snapshots]) for i in range(2)]
    e = [np.array([snapshot.planets[i].aei['e'] for snapshot in snapshots]) for i in range(2)]

    m = [np.array([snapshot.planets[i].m for snapshot in snapshots]) for i in range(2)]
    specific_e = [np.array([snapshot.planets[i].specific_e for snapshot in snapshots]) for i in range(2)]
    specific_h = [np.array([snapshot.planets[i].specific_h for snapshot in snapshots]) for i in range(2)]
    gas_mass_fraction = [np.array([snapshot.planets[i].gas_mass_fraction for snapshot in snapshots]) for i in range(2)]
    gas_mass = [np.array([snapshot.planets[i].gas_mass_fraction*snapshot.planets[i].m*ms/me for snapshot in snapshots]) for i in range(2)]
    utot = [np.array([snapshot.planets[i].utot for snapshot in snapshots]) for i in range(2)]
    ugas = [np.array([snapshot.planets[i].ugas for snapshot in snapshots]) for i in range(2)]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # Create orbit plots
    for i in range(2):
        am = np.ma.masked_where(a[i] == 0., a[i])
        em = np.ma.masked_where(a[i] == 0., e[i])
        ax.plot(t, am, color = c[i], linestyle = '-')  #semi-major axes
        ax.plot(t, (1.-em)*am, color = c[i], linestyle = ':')  #periapsis
        ax.plot(t, (1.+em)*am, color = c[i], linestyle = ':')  #apoapsis

    ax.set_xlim([0,max(t)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize=15)
    plt.ylabel(r'$a\ [\mathrm{AU}]$', fontsize=15)

    fname = fdir + 'orbits.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    # Create orbit plots
    for i in range(2):
        if ref_system.exists:
            plt.plot([min(t), max(t)], [ref_system.a[i], ref_system.a[i]], linestyle = '--', color = c[i])

    fname = fdir + 'orbits_compare.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_compare.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()

#    fig = plt.figure(figsize=(6,10))
#    fig.subplots_adjust(hspace = 0, wspace = 0, bottom = 0.2, top = 0.95, right = 0.95, left = 0.15)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # Create orbit plots
    for i in range(2):
        am = np.ma.masked_where(a[i] == 0., a[i])
        ax.plot(t, am, color = c[i], linestyle = '-')  #semi-major axes

    #plt.xlim([0,max(t)])
    ax.set_xlim([0,max(t)])
    ax.set_ylim([0.1,0.2])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize=15)
    #plt.xlabel(r'$t\ [\mathrm{yr}]$', fontsize=15)
    plt.ylabel(r'$a\ [\mathrm{AU}]$', fontsize=15)

    fname = fdir + 'orbits_a.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_a.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    # Create orbit plots
    for i in range(2):
        if ref_system.exists:
            plt.plot([min(t), max(t)], [ref_system.a[i], ref_system.a[i]], linestyle = '--', color = c[i])

    fname = fdir + 'orbits_compare_a.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_compare_a.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()

    fig = plt.figure()
#    fig = plt.figure(figsize=(6,10))
#    fig.subplots_adjust(hspace = 0, wspace = 0, bottom = 0.2, top = 0.95, right = 0.95, left = 0.15)
    ax = fig.add_subplot(1,1,1)

    # Create orbit plots
    for i in range(2):
        em = np.ma.masked_where(a[i] == 0., e[i])
        ax.plot(t, em, color = c[i], linestyle = '-')  #semi-major axes

    #plt.xlim([0,max(t)])
    ax.set_xlim([0,max(t)])
    ax.set_ylim([0.,0.25])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize=15)
    #plt.xlabel(r'$t\ [\mathrm{yr}]$', fontsize=15)
    plt.ylabel(r'$e$', fontsize=15)

    fname = fdir + 'orbits_e.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_e.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    # Create orbit plots
    for i in range(2):
        if ref_system.exists:
            plt.plot([min(t), max(t)], [ref_system.e[i], ref_system.e[i]], linestyle = '--', color = 'k')

    fname = fdir + 'orbits_compare_e.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_compare_e.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()

    fig = plt.figure()
#    fig = plt.figure(figsize=(6,10))
#    fig.subplots_adjust(hspace = 0, wspace = 0, bottom = 0.2, top = 0.95, right = 0.95, left = 0.15)
    ax = fig.add_subplot(1,1,1)

    # Create orbit plots
    P_ratio = np.ma.masked_where(a[0] == 0., (a[0]/a[1])**1.5)
    P_ratio = [x if x > 1. else 1./x for x in P_ratio]
    ax.plot(t, P_ratio, color = 'k', linestyle = '-')  #semi-major axes

    #plt.xlim([0,max(t)])
    ax.set_ylim([1., 1.6])
    ax.set_xlim([0,max(t)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize=15)
    #plt.xlabel(r'$t\ [\mathrm{yr}]$', fontsize=15)
    plt.ylabel(r'$P_\mathrm{outer}/P_\mathrm{inner}$', fontsize=15)

    fname = fdir + 'orbits_ratio.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_ratio.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    # Create orbit plots
    if ref_system.exists:
        plt.plot([min(t), max(t)], [ref_system.pratio, ref_system.pratio], linestyle = '--', color = 'k')

    fname = fdir + 'orbits_compare_pratio.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'orbits_compare_pratio.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()



    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # Create mass plots
    for i in range(2):
        mm = np.ma.masked_where(a[i] == 0., m[i])
        plt.plot(t, (mm-mm[0])*ms/me, color = c[i], linestyle = '-')  #change in mass

    plt.plot(t, m_bound_gas, color = 'k', linestyle = ':')   #mass bound to system
    plt.plot(t, m_unbound_gas, color = 'k', linestyle = '-.')  #mass in ejecta

    plt.axis([0,max(t),min(min(m[0]-m[0][0]),min(m[1]-m[1][0]))*ms/me,max(max(m[0]-m[0][0]),max(m[1]-m[1][0]))*ms/me])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$\Delta M/\mathrm{M}_\oplus$', fontsize=15)

    fname = fdir + 'mass_transfer.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'mass_transfer.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for i in range(2):
        plt.plot(t, gas_mass_fraction[i], color = c[i], linestyle = '-')  #mass in planet 1
    plt.xlim([0,max(t)])
    plt.ylim([0,np.max(gas_mass_fraction)*1.1])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$m_\mathrm{gas}/M$', fontsize=15)

    fname = fdir + 'Gas_Mass_Fraction.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'Gas_Mass_Fraction.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for i in range(2):
        plt.plot(t, 1.-gas_mass_fraction[i]/gas_mass_fraction[i][0], color = c[i], linestyle = '-')  # % mass lost
    plt.xlim([0,max(t)])
    plt.ylim([0,1.])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$m_\mathrm{gas\ lost}/m_{gas}$', fontsize=15)

    fname = fdir + 'Gas_Mass_Lost.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'Gas_Mass_Lost.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for i in range(2):
        plt.plot(t, gas_mass[i], color = c[i], linestyle = '-')  #mass in planet 1
    plt.xlim([0,max(t)])
    plt.ylim([0,np.max(gas_mass)*1.1])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$m_\mathrm{gas}\ [M_\mathrm{E}]$', fontsize=15)

    fname = fdir + 'Gas_Mass.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'Gas_Mass.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for i in range(2):
        plt.plot(t, utot[i], color = c[i], linestyle = '-')  #mass in planet 1
    plt.xlim([0,max(t)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$u_\mathrm{tot}$', fontsize=15)

    fname = fdir + 'utot.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'utot.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for i in range(2):
        plt.plot(t, ugas[i], color = c[i], linestyle = '-')  #mass in planet 1
    plt.xlim([0,max(t)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$u_\mathrm{gas}$', fontsize=15)

    fname = fdir + 'ugas.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'ugas.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    plt.plot(t, planet_energy, color = 'k', linestyle = '-')
    plt.plot([0,max(t)],[0.,0.],linestyle='--')
    plt.xlim([0,max(t)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$E$', fontsize=15)

    fname = fdir + 'planet_energy.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'planet_energy.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # plot specific energy
    for i in range(2):
        plt.plot(t, specific_e[i], color = c[i], linestyle = '-')
        if ref_system.exists:
            plt.plot([min(t), max(t)], [ref_system.specific_e[i], ref_system.specific_e[i]], color = c[i], linestyle = '--')

    plt.plot(t, specific_e[0]+specific_e[1], color = 'k', linestyle = '-')  #mass in planet 1
    if ref_system.exists:
        plt.plot([min(t), max(t)], [sum(ref_system.specific_e), sum(ref_system.specific_e)], color = 'k', linestyle = '--')

    plt.xlim([0,max(t)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$E$', fontsize=15)

    fname = fdir + 'specific_orbital_energy.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'specific_orbital_energy.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # plot specific angular momentum
    for i in range(2):
        plt.plot(t, specific_h[i], color = c[i], linestyle = '-')
        if ref_system.exists:
            plt.plot([min(t), max(t)], [ref_system.specific_h[i], ref_system.specific_h[i]], color = c[i], linestyle = '--')

    plt.plot(t, specific_h[0]+specific_h[1], color = 'k', linestyle = '-')  #mass in planet 1
    if ref_system.exists:
        plt.plot([min(t), max(t)], [sum(ref_system.specific_h), sum(ref_system.specific_h)], color = 'k', linestyle = '--')

    plt.xlim([0,max(t)])
    plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    plt.ylabel(r'$h$', fontsize=15)

    fname = fdir + 'specific_orbital_angularmomentum.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'specific_orbital_angularmomentum.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
    plt.clf()



