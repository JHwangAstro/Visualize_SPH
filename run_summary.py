
import cPickle as pickle
import math
import matplotlib.pyplot as plt
import os
import numpy as np
import plots
from scipy.spatial.distance import cdist
import subprocess
import sys

import calcs
import density_calculation
import orbit_calc
import write_data
from nbody_orbits import NBodyOrbits
from reference_system import Ref_system
import run_plots

class RunSummary:

    def __init__(self, run_id, fname_restart, sum_dir, run_dir, postSPH_dir, collision_ic, star_mass, radii = None, core_rad = None, compare = None, old_restart = None):

        hit_factor = 1.3

        run = self.get_run(run_dir, run_dir + fname_restart, run_dir + old_restart)
        if not run: return None
        #if not hasattr(run, 'energy'): run.energy = run.read_energy(run_dir)
        run.energy = run.read_energy(run_dir)
        snapshots = [run.snapshots[k] for k in sorted(run.snapshots)]

        # directory for run specific summary
        fdir = sum_dir + 'Run' + run_id + '/'
        if not os.path.exists(fdir): os.makedirs(fdir)

        # set collision details
        self.d_min, self.e_c, r0, r1 = collision_ic
        self.r = [r0, r1]

        # load in snapshot values
        self.m_star = snapshots[0].star.m
        self.energy_initial = run.energy[[k for k in sorted(run.energy)][0]]
        self.energy_final = run.energy[[k for k in sorted(run.energy)][-1]]
        self.gas_binding_energy = [[p.gas_binding_energy for p in snapshot.planets] for snapshot in snapshots]
        self.binding_energy = [[p.binding_energy for p in snapshot.planets] for snapshot in snapshots]

        # calculate run specific limits
        self.e_escape = calcs.find_e_escape([p.m for p in snapshots[0].planets], [p.m for p in snapshots[-1].planets], self.r, star_mass)
        self.r_Hill = [(_m/star_mass/3.)**(1./3.)*sum(self.r)/2. for _m in [p.m for p in snapshots[0].planets]]

        # do physics
        self.orbit = orbit_calc.get_bound_planet_aei(snapshots, max(self.r_Hill), sum(radii)*hit_factor)
        self.status = self.get_status(snapshots, sum_dir, postSPH_dir, run_id, radii, core_rad)
        self.alpha_ce = self.calculate_alpha_ce(fdir, snapshots)
        self.t_collisions = self.get_collision_times(snapshots, radii)

        self.orbit_energy, self.energy_dissipated, self.d_min_actual, self.d_min_actual2, self.ec_actual, self.resolved = \
                self.get_energy_dissipated(snapshots, run.energy, run.innermost_period, radii)

        print 'this run is resolved? ', self.resolved

        self.t_post_sph_collision = None
        if 'stable' in self.status: self.t_post_sph_collision = self.get_t_post_sph_collision(postSPH_dir, [radii, [_c*_r for _c, _r in zip(core_rad, radii)]])

        print 'time to first collision post sph', self.t_post_sph_collision

        self.energy_conservation = math.fabs((self.energy_initial['et'] - self.energy_final['et'])/self.energy_initial['et'])

        self.initial, self.final = self.set_parameters(snapshots, collision_ic, star_mass)

        # Report
        print 'energy dissipated: ', self.energy_dissipated

        print 'initial densities:', self.initial['d']
        print 'final densities:', self.final['d']
        print 'final period ratio: ', self.final['period_ratio']
        print 'initial rho ratio: ', self.initial['d_ratio']
        print 'final rho ratio: ', self.final['d_ratio']

        write_data.write_bound_planet_aei(snapshots, self.orbit, sum_dir, run_id, run.energy, run.innermost_period, radii, core_rad, self.r_Hill, self.e_escape, self.binding_energy)
        run_plots.plot_energy(snapshots, sum_dir, run_id)
        print 'done plotting'

        # read in orbits from mercury output
        aei_files = ['Kepler0.aei', 'Kepler1.aei']
        ce_files = ['Kepler0.clo', 'Kepler1.clo']
        preSPH_orbits = NBodyOrbits(run_dir, aei_files, ce_files, star_mass, radii)
        postSPH_orbits = NBodyOrbits(postSPH_dir, aei_files, ce_files, star_mass)

        # remake plots for individual runs
        plots.energy_plots(fdir, run.energy, run.innermost_period)
        plots.orbit_plots(fdir, snapshots, Ref_system(compare))
        plots.combined_orbits(fdir, snapshots, Ref_system(compare), preSPH_orbits, postSPH_orbits)
        #self.a_range = orbit_calc.find_a_range(self.m1f, self.m2f, self.a1f, self.a2f, self.e1f, self.e2f, self.MS)




    def get_run(self, d, f, f_old):

        def update_run(old_run):

            del old_run.particle_params

            # delete unused indices that inflate size of pickled save
            for k, snapshot in old_run.snapshots.iteritems():
                for planet in snapshot.planets:
                    del planet.particles

            #pickle the new run with new name, f, in d
            with open(f, 'wb') as output:
                pickle.dump(old_run, output, pickle.HIGHEST_PROTOCOL)

            return old_run


        print '\n', d

        # Check if there is a restart file and return pickled object
        # Load in old snapshots and update empty snapshot dictionary
        # Update old save if modified time is later than new save
        if os.path.isfile(f) and (os.path.getmtime(f) > os.path.getmtime(f_old)):
            with open(f, 'rb') as input:
                return pickle.load(input)

        else:
            if os.path.isfile(f_old):
                with open(f_old, 'rb') as input:
                    old_run = pickle.load(input)
                    return update_run(old_run)

            # cannot find either restart file
            else:
                if os.path.isfile(d + 'out0000.sph'): print 'restart file missing from ', f, f_old
                return None




    def set_parameters(self, snapshots, ic, mstar):

        def create_snapshot_summary(snapshot, snapshot0):

            rs_to_re = 6.957*10.**10./(6.371*10.**8.)
            ms = 1.98855*math.pow(10.,33.)
            me = 5.9722*10.**27.

            mydict = {'t': snapshot.t, \
                      'm': [p.m for p in snapshot.planets], \
                      'dm': [(p.m - p0.m)*ms/me for p, p0 in zip(snapshot.planets, snapshot0.planets)], \
                      'a': [p.aei['a'] for p in snapshot.planets], \
                      'e': [p.aei['e'] for p in snapshot.planets], \
                      'i': [p.aei['i'] for p in snapshot.planets], \
                      'aei1': [snapshot.planets[0].aei['a'], snapshot.planets[0].aei['e']], \
                      'aei2': [snapshot.planets[1].aei['a'], snapshot.planets[1].aei['e']], \
                      'r': [p.r for p in snapshot.planets], \
                      'v': [p.v for p in snapshot.planets], \
                      'd': [density_calculation.get_density(p) for p in snapshot.planets], \
                      'mgas': [p.gas_mass_fraction for p in snapshot.planets], \
                      'status': self.status, \
                      'planets': snapshot.planets, \
                      'd_min_actual': self.d_min_actual, \
                      'd_min_actual2': self.d_min_actual2, \
                      'ec_actual': self.ec_actual, \
                      'e_dissipated': self.energy_dissipated, \
                      't_collisions': self.t_collisions, \
                      'e_conservation': self.energy_conservation}

            mydict['C'] = orbit_calc.calculate_AMD_run(mstar, mydict)
            mydict['C_min'] = orbit_calc.calculate_AMD_min_run(mstar, mydict)

            mydict['flipped'] = mydict['a'][0] > mydict['a'][1]

            i,j = 0,1
            if mydict['flipped']: i,j = 1,0
            mydict['period_ratio'] = (mydict['a'][j]/mydict['a'][i])**(3./2.)*((mstar+mydict['m'][j])/(mstar+mydict['m'][i]))**(-1./2.)
            if cdist([snapshot.planets[0].r], [snapshot.planets[1].r])*rs_to_re < 6: mydict['period_ratio'] = 1.

            mydict['d_ratio'] = mydict['d'][1]/mydict['d'][0]

#            mydict['status'] = self.status

            mydict['d_min'], mydict['e_c'], mydict['r1'], mydict['r2'] = ic

#            mydict['planets'] = [p for p in snapshot.planets]

            return mydict


        return create_snapshot_summary(snapshots[0], snapshots[0]), create_snapshot_summary(snapshots[-1], snapshots[0])




    def get_status(self, snapshots, fdir, ndir, run_id, radii, core_rad_frac):

        rs = 6.957*10.**10.
        au = 1.496*10.**13.

        #snapshots = [run.snapshots[k] for k in sorted(run.snapshots)]
        core_rad = [_rad*_core for _rad, _core in zip(radii, core_rad_frac)]
        m_i = [p.m for p in snapshots[0].planets]
        m_f = [p.m for p in snapshots[-1].planets]
        sep_f = cdist([snapshots[-1].planets[1].r], [snapshots[-1].planets[0].r])*rs

        status = None

        # determine if merged by looking at final separation, first if it exists
        # and then if the energy is less than the critical escape energy
        # check if the periapsis separation is less than the core sep
        max_aei_i = max([i for i, orbit in enumerate(self.orbit) if orbit['aei']] + [0])
        aei = self.orbit[max_aei_i]['aei']
        final_e = min([x['e'] for x in self.orbit[-10:]])

        # first use stricter energy requirement for bound
        #if final_aei[0] and final_e < self.e_escape[1]:
        print '__energy: ', final_e, self.e_escape
        if final_e < self.e_escape[0]:
            if aei:
                periapsis_sep = aei['a']*(1.-aei['e'])*au
                print 'periapsis_sep: ', periapsis_sep, sum(core_rad)
                if periapsis_sep > sum(core_rad) and m_f[0] > 0.7*m_i[0] and m_f[1] > 0.7*m_i[1]:
                    status = 'bound'
                else:
                    status = 'merged'
                print status
            if m_f[0] < 0.7*m_i[0] or m_f[1] < 0.7*m_i[1]: status = 'merged'

        # use stricter energy requirement for scattering, and make sure no core was shed
        elif final_e >= self.e_escape[1] and m_f[0] > 0.7*m_i[0] and m_f[1] > 0.7*m_i[1] and \
            sep_f > sum(core_rad):
            # if not bound planet pair check mercury to see if stable
            status = self.check_nbody(ndir, radii)
            if not status: print 'Need to run nbody calculation'

        # if in the middle of the energy bounds
        elif self.e_escape[0] <= final_e < self.e_escape[1]:
            print '__ in middle energy region'
            if m_f[0] < 0.7*m_i[0] or m_f[1] < 0.7*m_i[1] or sep_f < sum(core_rad):
                status = 'merged'
                print status
            else:
                print '__checking: ', sep_f, max(self.r_Hill)*au, final_e
                if sep_f > max(self.r_Hill)*au:
                    status = self.check_nbody(ndir, radii)
                    if not status: print 'Need to run nbody calculation'

                elif aei:
                    periapsis_sep = aei['a']*(1.-aei['e'])*au
                    print 'periapsis_sep: ', periapsis_sep, sum(core_rad)
                    if periapsis_sep > sum(core_rad):
                    #if sep_f > sum(core_rad):
                        status = 'bound'
                    else:
                        status = 'merged'

                else:
                    print '__final else'
                    status = self.check_nbody(ndir, radii)
                    if not status: print 'Need to run nbody calculation'

        else:
            status = 'merged'


        if not status:
            print 'error: no status'
            sys.exit()
            status = 'merged'

        print status

        return status




    def check_nbody(self, ndir, radii):

        au = 1.496*10.**13.#/(6.371*10.**8.)

        if os.path.isfile(ndir + 'Kepler0.clo'):
            if os.path.getmtime(ndir + 'Kepler0.clo') < os.path.getmtime(ndir + 'ce.out'):
                clo_files = [x for x in os.listdir(ndir) if '.clo' in x]
                for clo in clo_files:
                    os.rename(ndir + clo, ndir + clo + '_backup')

        if not os.path.isfile(ndir + 'Kepler0.clo'):
            if not os.path.isfile(ndir + 'ce.out'): return None
            #run analysis scripts
            os.chdir(ndir)
            subprocess.call([ndir+'/close6'])

        f = open(ndir + 'Kepler0.clo','rb')
        t_d = [[line.split()[0], line.split()[2]] for line in f if len(line.split()) > 12]
        t_d = [[float(y) for y in x] for x in t_d if x[1][0] != '*' and x[0][0] != '*']

        if len(t_d) > 0:
            if min([x[1] for x in t_d])*au < sum(radii): return 'unstable'

        return 'stable'




    def get_t_post_sph_collision(self, ndir, radii):

        au = 1.496*10.**13.#/(6.371*10.**8.)

        f = open(ndir + 'Kepler0.clo','rb')
        t_d = [[line.split()[0], line.split()[2]] for line in f if len(line.split()) > 12]
        t_d = [[float(y) for y in x] for x in t_d if x[1][0] != '*' and x[0][0] != '*']

        print 'minimum distance in nbody: ', min([x[1] for x in t_d])*au
        print 'core size: ', sum(radii[1])
        print 'total size: ', sum(radii[0])

        def find_col(mydict, d):
            if len(mydict) > 0: return min(mydict)
            return None

        keys = ['rad', 'core']

        return {k: find_col({x[0]: x[1] for x in t_d if x[1]*au < sum(rad)}, rad) for k, rad in zip(keys, radii)}




    def calculate_alpha_ce(self, fdir, snapshots):

        g   = 6.67390*math.pow(10,-8)                               # gravitational constant [cgs]
        rs  = 6.9599*math.pow(10,10)                                # unit of length [cm]
        ms  = 1.9891*math.pow(10,33)                                # unit of mass [g]
        ts  = math.sqrt(math.pow(rs,3)/(g*ms))                      # unit of time [s]

        if self.status != 'bound': return None

        def orbital_energy(p):
            p0, p1 = p
            mu = p0.m_core*p1.core.m_core/(p0.m_core+p1.m_core)
            ek = mu/2.*np.linalg.norm(p0.core.v-p1.core.v)**2.*ms*rs*rs/ts/ts
            eg = -g*mu*(p0.m_core+p1.m_core)/np.linalg.norm(p0.core.r-p1.core.r)/rs*ms*ms
            return ek+eg

        # calculate bound energy of planets
        binding_energy0 = orbital_energy(snapshots[0].planets)
        d_core_binding_energy = [orbital_energy(snapshot.planets) - binding_energy0 for snapshot in snapshots]

        # get list of gas-binding energy
        binding_energy0 = sum([p.binding_energy['gas']['e'] for p in snapshots[0].planets])
        d_gas_binding_energy = [(sum([p.binding_energy['gas']['e'] for p in snapshot.planets]) - binding_energy0)*ms*rs*rs/ts/ts for snapshot in snapshots]

        alpha = [-d2/d1 if math.fabs(d1) > 0. else 0. for d1, d2 in zip(d_core_binding_energy, d_gas_binding_energy)]

        t = np.array([snapshot.unitless_time for snapshot in snapshots])
        plt.plot(t, alpha, color = 'k', linestyle = '-')
        plt.xlabel(r'$t/\mathrm{P}_\mathrm{inner}$', fontsize = 15)
        plt.ylabel(r'$\alpha$')

        plt.xlim([0, max(t)])

        fname = fdir + 'alpha_ce.png'
        plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
        fname = fdir + 'alpha_ce.ps'
        plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

        plt.clf()

        # calculate efficiency
        return alpha




    def get_energy_dissipated(self, snapshots, energy, P0, radii):

        g   = 6.67390*math.pow(10,-8)
        ms = 1.98855*math.pow(10.,33.)
        rs  = 6.9599*math.pow(10,10)
        ts  = math.sqrt(math.pow(rs,3)/(g*ms))
        yr  = 60.*60.*24.*365.24
        au = 1.496*10.**13.

        hit_factor = 1.3

        t = np.array([s.unitless_time for s in snapshots])
        d = np.array([orbit['sep'] for orbit in self.orbit])

        # get the times bracketing the first collision
        t_0, t_1 = orbit_calc.get_collision_times(t, d, max(self.r_Hill)*au, sum(radii)*hit_factor)[0]
        t_hit, dummy = orbit_calc.get_collision_times(t, d, 2.*sum(radii), sum(radii)*hit_factor)[0]

        # get parameters for collision using the SPH calculations
        d_min = min([_d for _t, _d in zip(t, d) if t_0 < _t < t_1])
        d_min_2 = orbit_calc.find_d_min({s.unitless_time: orbit['sep'] for s, orbit in zip(snapshots, self.orbit) if t_0 <= s.unitless_time <= t_1})
            # [s.unitless_time for s in snapshots if t_0 <= s.unitless_time <= t_1], \
            # [orbit['sep'] for s, orbit in zip(snapshots, self.orbit) if t_0 <= s.unitless_time <= t_1])

        # get the time just before the first change in mass
        m = {s.unitless_time: [p.m for p in s.planets] for s in snapshots}
        t_first_dm = max([k for k in m if m[k][0] == m[min(m)][0] and m[k][1] == m[min(m)][1]])
        ec_t = {s.unitless_time: x['e'] for s, x in zip(snapshots, self.orbit)}
        ec_d_min = ec_t[t_first_dm]

        # binding energy
        e_binding = {_t: sum([p.binding_energy['tot']['k'] + p.binding_energy['tot']['w'] for p in s.planets]) for _t, s in zip(t, snapshots)}

        # internal energy
        ei = np.array([energy[k]['ei'] for k in sorted(energy) if t_0 < k*ts/yr/P0 < t_1])#*ms*rs*rs/ts/ts

        orbital_energies = {_t: orbit_calc.get_planet_energy(s.planets, s.star) for _t, s in zip(t, snapshots)}

        # calculate difference in internal energy in this time frame
        d_ei = 0.
        if len(ei) > 0: d_ei = ei[-1] - ei[0]

        # calculate dissipation and sources of dissipation
        planet_energy = {s.unitless_time: x['e'] for s, x in zip(snapshots, self.orbit)}
        e_dissipated_guess = (d_ei + e_binding[t_1] - e_binding[t_0])*ms*rs*rs/ts/ts
        e_dissipated = planet_energy[t_first_dm] - min([planet_energy[_t] for _t in planet_energy if t_hit<_t<t_1])
        de_orbit = (orbital_energies[t_1]['k'] + orbital_energies[t_1]['w']) - (orbital_energies[t_0]['k'] + orbital_energies[t_0]['w'])
        #alpha = (de_orbit + e_dissipated_guess)/de_orbit
        alpha = (e_dissipated - e_dissipated_guess)/e_dissipated

        print 'orbital/dissipated energies: ', e_dissipated, e_dissipated_guess, alpha, (e_binding[t_1] - e_binding[t_0])*ms*rs*rs/ts/ts, d_ei*ms*rs*rs/ts/ts

        # returns orbital_energies, dissipated energies, lowest snapshot dmin, interpolated dmin, if first pass is resolved
        return orbital_energies, e_dissipated, d_min, d_min_2, ec_d_min, (t_1 != max(t))
#        return (d_ei + e_binding[t_1] - e_binding[t_0])*ms*rs*rs/ts/ts, min(d_min)
        return 0.




    def get_collision_times(self, snapshots, radii):

        hit_factor = 1.3

        t = np.array([s.unitless_time for s in snapshots])
        d = np.array([orbit['sep'] for orbit in self.orbit])

        return orbit_calc.get_collision_times(t, d, sum(radii), sum(radii)*hit_factor)




    def get_mass_lost(self):

        t0, t1 = self.t_collisions[0]

        m_t = {orbit['t']: orbit['m'] for orbit in self.orbits}

        return m_t[t1]









