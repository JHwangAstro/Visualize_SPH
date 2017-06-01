
import cPickle as pickle
from datetime import datetime
import numpy as np
import os
import operator
import string
import math

# user defined libraries
from nbody import NBody
from particle import ParticleParameters
import plots
from snapshot import Snapshot
from reference_system import Ref_system


class Run:

    def __init__(self, args):

        # Create new snapshots
        # Create a dictionary of snapshots
        self.snapshot_keys = np.arange(args.nstart, args.numss, args.nskip)
        self.parent_numbers = {k: max([args.nstart, k-args.nskip]) for k in self.snapshot_keys}
        self.snapshots = {k: Snapshot(k) for k in self.snapshot_keys}

        # Get run parameters
        #self.read_sph_parameters(args.sdir, self.snapshots[args.nstart], args.nstart)
        self.read_sph_parameters(args.sdir, self.snapshots[args.nstart], 0)

        # Initialize empty_snapshots as all snapshots
        self.empty_snapshots = {k: snapshot for k, snapshot in self.snapshots.iteritems()}

        # Repopulate snapshots that exist in restart file
        if args.restart == 1: self.restart(args.sdir + args.rname)

        # Create snapshots in ascending order, using the previous snapshot
        for k, snapshot in sorted(self.empty_snapshots.items(), key=operator.itemgetter(0)):
            # Set snapshot filename
            fname = args.sdir + 'out' + format(k, "04") + '.sph'
            print 'snapshot: ', fname
            parent_snapshot = self.snapshots[self.parent_numbers[k]]
            snapshot.initialize(fname, args, parent_snapshot, self.particle_params)
            # Define initial internal energies, innermost period, hill radii of planets,
            # and outermost orbit of first snapshot.
            self.find_orbit_parameters(snapshot)

        self.energy = self.read_energy(args.sdir)
        self.write_restart(args.sdir + args.rname)
        self.summarize_run(args.fdir, args.compare, args.sdir)
        self.initialize_nbody(snapshot, args.mdir, args.ndir)
        self.initialize_thermal_evolution(snapshot)




    def read_sph_parameters(self, sdir, snapshot, first_snapshot):

        # Read in first snapshot
        fname = sdir + 'out' + format(first_snapshot, "04") + '.sph'
        header, particles = snapshot.read_snapshot(fname, {})
        u0 = np.array([particle.u for k, particle in particles.iteritems()])

        # Read in calculation parameters
        self.rho_e, self.rho_c, self.rho_m, self.rho_i, self.K_e, self.gamma_e = \
                np.loadtxt(sdir+'sph.rhoparams', skiprows=0, unpack=True, usecols=(0,1,2,3,4,5))
        self.rho_e2, self.rho_c2, self.rho_m2, self.rho_i2, self.K_e2, self.gamma_e2 = \
                np.loadtxt(sdir+'sph.rhoparams2', skiprows=0, unpack=True, usecols=(0,1,2,3,4,5))

        core = np.array([1 if ((particle.cc > 2) & (particle.rho > self.rho_e2)) | ((particle.cc < 3) & (particle.rho > self.rho_e)) \
                else 0 for k, particle in particles.iteritems()])
        core2 = np.array([1 if ((particle.cc > 2) & (particle.rho > self.rho_m2)) | ((particle.cc < 3) & (particle.rho > self.rho_m)) \
                else 0 for k, particle in particles.iteritems()])
        K = np.array([self.K_e2 if (particle.cc > 2) else self.K_e for k, particle in particles.iteritems()])
        gam = np.array([self.gamma_e2 if (particle.cc > 2) else self.gamma_e for k, particle in particles.iteritems()])

        # Read in particle data
        f = open(sdir+'sph.particle_comp0', 'rb')
        dtype = np.dtype("i4,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,i4")
        data = np.fromfile(f, dtype = dtype, count = header.ntot)
        f.close()

        self.particle_params = {i: ParticleParameters(i, dat, u0[i], core[i], core2[i], K[i], gam[i]) for i, dat in enumerate(data)}
        print 'particle_params: ', len(self.particle_params)




    def restart(self, f):

        # Check if there is a restart file and return pickled object
        # Load in old snapshots and update empty snapshot dictionary
        if os.path.isfile(f):
            with open(f, 'rb') as input:
                old_run = pickle.load(input)
                keys_a = set(self.snapshots.keys())
                keys_b = set(old_run.snapshots.keys())
                for k in keys_a & keys_b:
                    self.snapshots[k] = old_run.snapshots[k]
                    self.empty_snapshots.pop(k, None)
        else:
            print 'restart file missing, proceeding from beginning'
            return None




    def write_restart(self, f):

        with open(f, 'wb') as output:
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)




    def find_orbit_parameters(self, snapshot):

        au_to_rs = (1.496*10.**13.)/(6.957*10.**10.)

        # Sets various initial system values
        if not hasattr(self, 'innermost_period'):
            for planet in snapshot.planets:
                if not hasattr(self, 'innermost_period'):
                    self.innermost_period = planet.aei['a']**1.5
                else:
                    if planet.aei['a']**1.5 < self.innermost_period:
                        self.innermost_period = planet.aei['a']**1.5

        if not hasattr(self, 'hill_radius'):
            self.hill_radius = {string.ascii_lowercase[i]: \
                    planet.aei['a']*(1.-planet.aei['e'])*(planet.m/(3*snapshot.star.m))**(1./3.)*au_to_rs \
                    for i, planet in enumerate(snapshot.planets)}

        if not hasattr(self, 'outer_orbit'):
            for planet in snapshot.planets:
                if not hasattr(self, 'outer_orbit'):
                    self.outer_orbit = planet.aei['a']
                else:
                    if planet.aei['a'] > self.outer_orbit:
                        self.outer_orbit = planet.aei['a']
                self.outer_orbit = self.outer_orbit*au_to_rs

        snapshot.unitless_time = snapshot.set_unitless_time(self.innermost_period)
        snapshot.hill_radius = snapshot.set_hill_radius(self.hill_radius)
        snapshot.outer_orbit = snapshot.set_outer_orbit(self.outer_orbit)



    def summarize_run(self, fdir, ref_system_name, sdir):

        snapshots = [self.snapshots[k] for k in sorted(self.snapshots)]
        plots.orbit_plots(fdir, snapshots, Ref_system(ref_system_name))

        fname = sdir + 'out0000.sph'
        snapshots[0].create_planet_profiles(fname, self.particle_params)
        fname = sdir + 'out' + format(len(snapshots)-1, "04") + '.sph'
        snapshots[len(snapshots)-1].create_planet_profiles(fname, self.particle_params)

        self.write_summary_files(fdir, snapshots[0], snapshots[len(snapshots)-1])
        self.write_planet_profiles(fdir, snapshots[0], snapshots[len(snapshots)-1])




    def write_summary_files(self, fdir, initial_snapshot, final_snapshot):

        ms = 1.98855*math.pow(10.,33.)                          #Solar mass in g
        me = 5.972*math.pow(10.,27.)                            #Earth mass in g

        column_names = ['Planet', 'Mass', 'Semi-major axis', 'Eccentricity', 'Inclination', 'Argument of Periastron', 'Ascending Node', 'True Anomaly']

        initial_values = [[planet.name, planet.m, planet.aei['a'], planet.aei['e'], planet.aei['i'], planet.aei['g'], planet.aei['n'], planet.aei['M']] \
                for planet in initial_snapshot.planets]
        final_values = [[planet.name, planet.m, planet.aei['a'], planet.aei['e'], planet.aei['i'], planet.aei['g'], planet.aei['n'], planet.aei['M']] \
                for planet in final_snapshot.planets]

        f = open(fdir+'Summary.txt', 'w')
        f.write('Initial'+'\n')

        for name in column_names:
            f.write(name+'\t')
        f.write('\n')

        for value in initial_values:
            f.write(str(value)+'\t')
        f.write('\n')
        f.write('\n')


        f.write('Final'+'\n')

        for name in column_names:
            f.write(name+'\t')
        f.write('\n')

        for value in final_values:
            f.write(str(value)+'\t')
        f.write('\n')
        f.write('\n')


        column_names = ['dM_1', 'dM_2', 'dM_env', 'dM_eject', 'dE_1', 'dE_2', 'dE_tot', 'dh_1', 'dh_2', 'dh_tot']

        p = [planet for planet in initial_snapshot.planets]
        initial_values = [p[0].m*ms/me, p[1].m*ms/me, initial_snapshot.bound_mass, initial_snapshot.ejected_mass, \
                p[0].specific_e, p[1].specific_e, p[0].specific_e+p[1].specific_e, \
                p[0].specific_h, p[1].specific_h, p[0].specific_h+p[1].specific_h]

        p = [planet for planet in final_snapshot.planets]
        final_values = [p[0].m*ms/me, p[1].m*ms/me, final_snapshot.bound_mass, final_snapshot.ejected_mass, \
                p[0].specific_e, p[1].specific_e, p[0].specific_h+p[1].specific_h, \
                p[0].specific_h, p[1].specific_h, p[0].specific_h+p[1].specific_h]

        f.write('Change from Initial to Final\n')

        for name in column_names:
            f.write(name+'\t')
        f.write('\n')

        for val_i, val_f in zip(initial_values, final_values):
            f.write(str(val_f-val_i)+'\t')
        f.write('\n')
        f.write('\n')


        f.write('Change from Initial to Final in percentages\n')

        for name in column_names:
            f.write(name+'\t')
        f.write('\n')

        for val_i, val_f in zip(initial_values, final_values):
            f.write(str((val_f-val_i)/val_i)+'\t')
        f.write('\n')
        f.write('\n')



        for title, snapshot in zip(['Initial Spins', 'Final Spins'], [initial_snapshot, final_snapshot]):

            f.write(title+'\n')

            for p, prefix in zip([planet for planet in snapshot.planets], ['j1', 'j2']):

                column_names = [prefix, prefix+'_core', prefix+'_gas', prefix+'_specific', \
                        prefix+'_specific_core', prefix+'_specific_gas', prefix+'_vector_core', prefix + '_vector_gas']

                for name in column_names:
                    f.write(name+'\t')
                f.write('\n')

                values = [p.j_gas + p.j_core, p.j_core, p.j_gas, p.specific_j_gas + p.specific_j_core, \
                        p.specific_j_core, p.specific_j_gas, p.j_unit_vector_core, p.j_unit_vector_gas]

                for val in values:
                    f.write(str(val)+'\t')
                f.write('\n')
                f.write('\n')

        f.close()




    def write_planet_profiles(self, fdir, initial_snapshot, final_snapshot):

        fnames = ['initial_', 'final_']
        suffixes = ['_p.txt', '_r.txt']

        for snapshot, fname in zip([initial_snapshot, final_snapshot], fnames):
            for i, p in enumerate(snapshot.planets):
                for profile, suffix in zip([p.profile.r_profile, p.profile.p_profile], suffixes):

                    f = open(fdir + fname + str(i) + suffix, 'w')
                    keys = [k for k in profile]
                    header = str('')
                    for k in keys:
                        header = header + k + str('\t')

                    #_fname = fdir + fname + str(i) + suffix

                    my_list = [profile[k] for k in keys]
                    np.savetxt(fdir + fname + str(i) + suffix, np.transpose([profile[k] for k in keys]), fmt = '%1.4e', header = header)

                    # write header
                    #for k in keys:
                    #    f.write(k + '\t')

                    #for v in [profile[k] for k in keys]:
                    #    f.write('\n')
                    #    for _v in v:
                    #        f.write(str(_v) + '\t')




    def read_energy(self, sdir):

        energy_files = [x for x in os.listdir(sdir) if 'energy' in x]

        my_dict = {}

        for fname in energy_files:
            data = np.loadtxt(sdir + '/' + fname, skiprows=0, unpack = True, usecols=(0,1,2,3,4,5,6))
            for t, eg, ek, ei, et, s, a in zip(*data):
                my_dict[t] = {'eg': eg, 'ek': ek, 'ei': ei, 'et': et, 's': s, 'a': a}

        return my_dict




    def initialize_nbody(self, snapshot, mdir, ndir):

        nbody = NBody(snapshot)
        nbody.run_calculation(mdir, ndir)







