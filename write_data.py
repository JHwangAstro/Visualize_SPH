
import math
import matplotlib.pyplot as plt
import os
import numpy as np

import orbit_calc


def write_bound_planet_aei(snapshots, aei_in, fdir, run_id, energy, P0, rad, c_rad, r_Hill, e_escape, e_binding):

    g   = 6.67390*math.pow(10,-8)
    yr  = 60.*60.*24.*365.24
    au = 1.496*10.**13.
    me = 5.972*math.pow(10.,27.)
    re = 6.371*10.**8.
    ms = 1.98855*math.pow(10.,33.)
    rs  = 6.9599*math.pow(10,10)
    ts  = math.sqrt(math.pow(rs,3)/(g*ms))

    hit_factor = 1.3

    core_rad = [_rad*_c_rad/re for _rad, _c_rad in zip(rad, c_rad)]

    aei = [x['aei'] for x in aei_in if x['aei']]
    sep = [x['sep']/re for x in aei_in if x['aei']]
#    aei2 = [x['aei'][1] for x in aei_in if x['aei'][1]]
    t = [snapshot.unitless_time for x, snapshot in zip(aei_in, snapshots) if x['aei']]

    orbits_dir = fdir + 'Orbits'
    if not os.path.exists(orbits_dir): os.makedirs(orbits_dir)

    # if there is a period where the planets are bound
    if len(aei) > 0:

        p1, p2 = snapshot.planets

        # write the header
        fname = orbits_dir + '/Run_' + run_id + '_bound_aei.txt'
        f = open(fname, 'w')
        columns  = ['m1', 'm1_gas', 'm2', 'm2_gas', 'a', 'e']#, 'i', 'n', 'g', 'M']
        for name in columns:
            f.write(name + '\t\t\t')
        f.write('\n')

        for aei_t in aei:
            f.write(str("%.6f"%round(p1.m*ms/me,6)) + '\t')
            f.write(str("%.6f"%round(p1.m*p1.gas_mass_fraction*ms/me,6)) + '\t')
            f.write(str("%.6f"%round(p2.m*ms/me,6)) + '\t')
            f.write(str("%.6f"%round(p2.m*p2.gas_mass_fraction*ms/me,6)) + '\t')
            for k in ['a', 'e']:#, 'i', 'g', 'n', 'M']:
                f.write(str("%.6f"%round(aei_t[k],6)) + '\t')
            f.write('\n')
        f.close()

        # plot planet orbits
        a = np.array([x['a'] for x in aei])
        e = np.array([x['e'] for x in aei])

        plt.plot([0, max(t)], [sum(core_rad), sum(core_rad)], color = 'r', linestyle = ':')  #apoapsis
        plt.plot(t, a*au/re, color = 'b', linestyle = '-')  #semi-major axes
        plt.plot(t, sep, color = 'k', linestyle = '--')
        plt.plot(t, (1.-e)*a*au/re, color = 'b', linestyle = ':')  #periapsis
        plt.plot(t, (1.+e)*a*au/re, color = 'b', linestyle = ':')  #apoapsis

        plt.xlim([0,max(t)])
        plt.yscale('log')
        plt.xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
        plt.ylabel(r'$a\ [R_\mathrm{E}]$', fontsize=15)

        fname = orbits_dir + '/Run_' + run_id + '_bound_planets_orbits.png'
        plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
        fname = orbits_dir + '/Run_' + run_id + '_bound_planets_orbits.ps'
        plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)
        plt.clf()


    # plot the orbital energy
    t = [snapshot.unitless_time for snapshot in snapshots]
    planet_energy = [x['e'] for x in aei_in]
    sep = [x['sep']/re for x in aei_in]

    f, ax1 = plt.subplots()

    ax1.plot(t, planet_energy, color = 'k', linestyle = '-')
    ax1.plot([0,max(t)], [0.,0.], linestyle='--')
    ax1.set_xlim([0,max(t)])
    ax1.set_xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    ax1.set_ylabel(r'$E$', fontsize=15)

    ax2 = ax1.twinx()
    ax2.plot(t, sep, color = 'k', linestyle = '--')
    ax2.set_ylabel(r'$d/R_\mathrm{E}$')
    f.tight_layout()

    fname = orbits_dir + '/Run_' + run_id + '_planet_energy.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = orbits_dir + '/Run_' + run_id + '_planet_energy.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()


    dcrit = max(r_Hill)*au/re*3.
    print 'dcrit: ', dcrit
    t = np.array([snapshot.unitless_time for snapshot, _d in zip(snapshots, sep) if _d < dcrit])
    e_orbit = np.array([snapshot.orbital_energy*ms*rs*rs/ts/ts for snapshot, _d in zip(snapshots, sep) if _d < dcrit])
    planet_energy = np.array([x['e'] for x, _d in zip(aei_in, sep) if _d < dcrit])
#    planet_energy = planet_energy - planet_energy[0]

    #_e_binding = np.array([sum(eb) for eb, _d in zip(e_binding, sep) if _d < dcrit])*ms*rs*rs/ts/ts
    #_e_binding = _e_binding - _e_binding[0]

    sep = [x['sep']/re for x, _d in zip(aei_in, sep) if _d < dcrit]
    ei = np.array([energy[k]['ei'] for k in sorted(energy) if min(t) < k*ts/yr/P0 < max(t)])
    ei = (ei - ei[0])*ms*rs*rs/ts/ts
    ei_t = [k*ts/yr/P0 for k in sorted(energy) if min(t) < k*ts/yr/P0 < max(t)]
    eg = np.array([energy[k]['eg'] for k in sorted(energy) if min(t) < k*ts/yr/P0 < max(t)])
    eg = (eg - eg[0])*ms*rs*rs/ts/ts
    ek = np.array([energy[k]['ek'] for k in sorted(energy) if min(t) < k*ts/yr/P0 < max(t)])
    ek = (ek - ek[0])*ms*rs*rs/ts/ts



    print min(ei), max(ei), min(planet_energy), max(planet_energy), e_escape
    t_in = np.array([snapshot.unitless_time for snapshot in snapshots])
    d = np.array([orbit['sep'] for orbit in aei_in])
    t_0, t_1 = orbit_calc.get_collision_times(t_in, d, max(r_Hill)*au, sum(rad)*hit_factor)[0]
    t_hit, t_leave = orbit_calc.get_collision_times(t_in, d, 2.*sum(rad), sum(rad)*hit_factor)[0]

    # get the time just before the first change in mass
    m = {s.unitless_time: [p.m for p in s.planets] for s in snapshots}
    t_first_dm = max([k for k in m if m[k][0] == m[min(m)][0] and m[k][1] == m[min(m)][1]])

    pe = {x['e']: _t for _t, x in zip(t_in, aei_in) if t_hit < _t < t_1}
    t_min = pe[min(pe)]

    f, ax1 = plt.subplots()

    ax1.plot(t, planet_energy, color = 'k', linestyle = '-')
    #ax1.plot(t, e_orbit - e_orbit[0], color = 'r', linestyle = '-')
    ax1.plot(ei_t, ei, color = 'r', linestyle = '-')
    #ax1.plot(ei_t, ek, color = 'b', linestyle = ':')
    #ax1.plot(ei_t, eg, color = 'k', linestyle = ':')
    #ax1.plot(t, _e_binding, color = 'b', linestyle = '-')
    for _e_escape in e_escape:
        ax1.plot([0,max(t)], [_e_escape, _e_escape], linestyle='--', color = 'r')
    ax1.set_xlim([0,max(t)])
    ax1.set_xlabel(r'$t/\mathrm{P}_\mathrm{0}$', fontsize=15)
    ax1.set_ylabel(r'$E$', fontsize=15)

    ax2 = ax1.twinx()

    ax2.plot([t_min, t_min], [min(sep), max(sep)], c = 'b', linestyle = '--')
    ax2.plot([t_first_dm, t_first_dm], [min(sep), max(sep)], c = 'b', linestyle = '--')

    ax2.plot(t, sep, color = 'k', linestyle = '--')
    ax2.plot([0, max(t)], [sum(rad)/re, sum(rad)/re], linestyle = '--', color = 'b')
    ax2.set_ylabel(r'$d/R_\mathrm{E}$')
    f.tight_layout()

    fname = orbits_dir + '/Run_' + run_id + '_planet_energy_ce.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = orbits_dir + '/Run_' + run_id + '_planet_energy_ce.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()







