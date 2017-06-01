
import cPickle as pickle
import math
import numpy as np
from scipy.spatial.distance import cdist
import os

# user defined libraries
import calcs
import collision_plots
import mass_radius_tables
import orbit_calc
import plots
from run import Run
from run_summary import RunSummary
import status_specific_calcs
import summary_plots
import write_summary


class Summary:

    def __init__(self, args):

        # objects in args
        # ndir - universal folder for dynamical calculations
        # rname - name of restart file to read in
        # compare - observed system to compare to
        # ids - IDs of runs to aggregate
        # fname_ic - path to initial conditions file
        # sum_dir - directory to store aggregated output
        # sdir - list of run directories
        # ndir - list of directories for post-SPH dynamical calculations
        # star_mass - mass of host star

        # GENERALLY UNUSED
        # prefix - prefix of run directories
        # num_p - number of SPH planets
        # num_p2 - number of point-mass planets
        # bdir - base directory; normally unused
        # rdir - directory holding all runs to be analyzed

        # dictionary where k = run id and v = [d_min, E_collision]
        collision_parameters = self.create_collision_dictionary(args.fname_ic)

        self.cosmetics = self.init_cosmetics()

        # dictionary holding run summary object for each run
        # for each run directory, get summary object and if flagged; rerun analysis
        if os.path.isfile(args.sum_name) and not args.rerun:
            with open(args.sum_name, 'rb') as input:
                run_summary = pickle.load(input)
        #        run_summary = {k: v for k, v in run_summary.iteritems() if v.energy_dissipated > 0.}

        else:
            run_summary = {run_id: RunSummary(run_id, args.rname, args.sum_dir, run_directory, nbody_dir, \
                    collision_parameters[run_id], args.star_mass, args.structure['r'], args.structure['c'], args.compare, args.rname_old) \
                    for run_directory, nbody_dir, run_id in zip(args.sdir, args.ndir, args.ids)}

            run_summary = {k: v for k, v in run_summary.iteritems() if v.status}

            with open(args.sum_name, 'wb') as output:
                pickle.dump(run_summary, output, pickle.HIGHEST_PROTOCOL)

        self.mass_dmin = calcs.get_mass_dmin(args.relax_run)

        self.write_summary_files(args.sum_dir, run_summary, args.structure)

        self.e_conservation, self.e_efficiency, self.e_thermal_increase_unitless, self.e_thermal_increase_cgs = \
                self.energy_calculations(args.sum_dir, run_summary, args.sum_dir + 'energy_summary.txt', args.structure)

        # get models for runs if they have resolved
        print '______________Unresolved runs: '
        for k in run_summary:
            if not run_summary[k].resolved: print k
        self.models = self.get_models({k: run_summary[k] for k in run_summary if run_summary[k].resolved}, args.structure)

        self.status_specific_analyses(args.sum_dir, run_summary, args.star_mass, args.compare, args.structure, args.set_name)

        self.summarize_all_runs(args.sum_dir, run_summary, args.compare, args.structure, args.star_mass, args.set_name)




    def init_cosmetics(self):

        # edgecolors
        ec_dict  = {0: {'bound': 'k', 'merged': 'k', 'stable': 'k', 'unstable': 'k'},
                    1: {'bound': 'k', 'merged': 'k', 'stable': 'r', 'unstable': 'r'}}

        # facecolors
        fc_dict  = {'bound': 'k', 'merged': 'k', 'stable': 'y', 'unstable': 'none'}

        # markers
        m_dict = {'bound': '*', 'merged': 'x', 'stable': 'o', 'unstable': 'o'}

        return {'ec': ec_dict, 'fc': fc_dict, 'marker': m_dict}




    def set_cosmetics(self, runs):

        status = [runs[k].status for k in runs]
        ec = [self.cosmetics['ec'][runs[k].final['flipped']][runs[k].status] for k in runs]
        fc = [self.cosmetics['fc'][runs[k].status] for k in runs]
        markers = [self.cosmetics['marker'][runs[k].status] for k in runs]

        return ec, fc, markers




    def get_models(self, runs, structure):

        me = 5.972*10.**27
        ms = 1.989*10.**33


        print 'fitting energy dissipated as a function of d'

        # first prune data too close to core-core collisions
        # get the distance of closest approach from each run (cubic fit to SPH)
        r1, r2 = structure['r']
        c1, c2 = structure['c_sph']
        dcrit = 1.05*(r1*c1+r2*c2)/(r1+r2)

        d = [runs[k].d_min_actual2['cubic']/sum(structure['r']) for k in runs]
        e = [runs[k].ec_actual/structure['eb'] for k, _d in zip(runs, d) if _d > dcrit]
        de = [runs[k].energy_dissipated/structure['eb'] for k, _d in zip(runs, d) if _d > dcrit]
        d = [runs[k].d_min_actual2['cubic']/sum(structure['r']) for k, _d in zip(runs, d) if _d > dcrit]

        def f_de(x, a, gam, c):
            x1, x2 = x
            if gam < 1.: return 10.**30.
            model = a*np.exp(-gam*(x1 - dcrit)) + c
            return model

        model_de = calcs.fit_2d(de, d, e, f_de, [.5, 5., 0.], ['a', 'gamma', 'c'])



        print 'fitting nbody distance as a function of d'

        d = [runs[k].d_min_actual2['cubic']/sum(structure['r']) for k in runs]
        e = [runs[k].ec_actual/structure['eb'] for k in runs]
        d_nbody = [runs[k].d_min for k, _d in zip(runs, d)]
        e_nbody = [runs[k].e_c for k, _d in zip(runs, d)]

        def f_d_nbody(x, a, c):
            x1, x2 = x
            model = a*x1 + c
            return model

        model_d_nbody = calcs.fit_2d(d_nbody, d, e, f_d_nbody, [1., 0.], ['a', 'c'])



        print 'fitting planet-mass as a function of d'

        d = [runs[k].d_min_actual2['cubic']/sum(structure['r']) for k in runs if 'stable' in runs[k].status]
        e = [runs[k].ec_actual/structure['eb'] for k in runs if 'stable' in runs[k].status]
        dm = [[-runs[k].final['dm'][i] for k in runs if 'stable' in runs[k].status] for i in range(2)]
        #dm = [[orbit_calc.get_mass_lost(runs[k])[i]*ms/me for k in runs if 'stable' in runs[k].status] for i in range(2)]

        def f_mgas(x, a, gam, c):
            x1, x2 = x
            model = a*np.exp(-gam*(x1 - dcrit)) + c
            return model

        model_m = {i: calcs.fit_2d(p_dm, d, e, f_mgas, [0.01, -2., 0.], ['a', 'gamma', 'c']) for i, p_dm in enumerate(dm)}


        print 'fitting total mass as a function of d'

        d = [runs[k].d_min_actual2['cubic']/sum(structure['r']) for k in runs if 'stable' in runs[k].status]
        e = [runs[k].ec_actual/structure['eb'] for k in runs if 'stable' in runs[k].status]
        dm = [-sum(runs[k].final['dm']) for k in runs if 'stable' in runs[k].status]

        def f_mgas(x, a, gam, c):
            x1, x2 = x
            model = a*np.exp(-gam*(x1 - dcrit)) + c
            return model

        model_mtot = calcs.fit_2d(dm, d, e, f_mgas, [0.01, -2., 0.], ['a', 'gamma', 'c'])



        print 'fitting energy dissipated as a function of outside mass'
        # Fit using m_outside as a predictor
        d = [runs[k].d_min_actual2['cubic']/sum(structure['r']) for k in runs]
        m_outside = [self.mass_dmin(_d*sum(structure['r']))/me for _d, k in zip(d, runs) if _d > dcrit]
        e = [runs[k].ec_actual/structure['eb'] for k, _d in zip(runs, d) if _d > dcrit]
        de = [runs[k].energy_dissipated/structure['eb'] for k, _d in zip(runs, d) if _d > dcrit]

        def f_de_m(x, a, gam, c):
            x1, x2 = x
            model = a*np.exp(gam*x1) + c
            return model

        model_de_m = calcs.fit_2d(de, m_outside, e, f_de_m, [.5, 5., 0.], ['a', 'gamma', 'c'])


        print 'fitting planet-mass as a function of m'

        m_outside = [self.mass_dmin(_d*sum(structure['r'])) for _d, k in zip(d, runs) if 'stable' in runs[k].status]
        e = [runs[k].ec_actual/structure['eb'] for k in runs if 'stable' in runs[k].status]
        dm = [[-runs[k].final['dm'][i] for k in runs if 'stable' in runs[k].status] for i in range(2)]

        def f_mgas_m(x, a, gam, c):
            x1, x2 = x
            model = a*np.exp(gam*x1) + c
            return model

        model_m_m = {i: calcs.fit_2d(p_dm, m_outside, e, f_mgas, [0.01, 2., 0.], ['a', 'gamma', 'c']) for i, p_dm in enumerate(dm)}


        print 'fitting total mass as a function of m'

        m_outside = [self.mass_dmin(_d*sum(structure['r'])) for _d, k in zip(d, runs) if 'stable' in runs[k].status]
        e = [runs[k].ec_actual/structure['eb'] for k in runs if 'stable' in runs[k].status]
        dm = [-sum(runs[k].final['dm']) for k in runs if 'stable' in runs[k].status]

        def f_mgas_m(x, a, gam, c):
            x1, x2 = x
            model = a*np.exp(gam*x1) + c
            return model

        model_mtot_m = calcs.fit_2d(dm, m_outside, e, f_mgas, [0.01, 2., 0.], ['a', 'gamma', 'c'])


        return {'de_d': model_de, 'dnbody_d': model_d_nbody, 'pdm_d': model_m, 'mtot_d': model_mtot, \
                'de_m': model_de_m, 'pdm_m': model_m_m, 'mtot_m': model_mtot_m}




    def write_summary_files(self, fdir, runs, structure):

        me = 5.9722*10.**27.
        ms = 1.99*10.**33.

        # list of [title, key], where key is the key for the run_summary dictionary
        params_to_write = [['m', 'm'], ['mgas', 'mgas'], ['a', 'a'], ['e', 'e'], ['rho', 'd'], \
                ['P2/P1', 'period_ratio'] , ['d1/d2', 'd_ratio'], ['dmin', 'd_min'], ['ec', 'e_c'], ['status', 'status']]

        write_summary.write_system_parameters(fdir + 'system_parameters_final.txt', [(k, run.final) for k, run in runs.iteritems()], params_to_write)
        write_summary.write_system_parameters(fdir + 'system_parameters_initial.txt', [(k, run.initial) for k, run in runs.iteritems()], params_to_write)


        #__write properties
        params_to_write = [['dmin', 'd_min', 1., None], ['ec', 'e_c', 1., None], ['dm', 'dm', 1., None], \
                ['E_d', 'e_dissipated', structure['eb'], 'small'], ['dE', 'e_conservation', 1., 'sci'], ['status', 'status', 1., 'string']]

        initial = [(k, runs[k].initial) for k in runs]
        final = [(k, runs[k].final) for k in runs]

        write_summary.write_system_parameters_formatted(fdir + 'system_parameters_final_formatted.txt', final, params_to_write)
        write_summary.write_system_parameters_formatted(fdir + 'system_parameters_initial_formatted.txt', initial, params_to_write)
        write_summary.write_system_parameters_formatted(fdir + 'system_parameters_final_paper.txt', final, params_to_write, paper = 1)

        #__write results
        params_to_write = [['dmin', 'd_min', 1., None], ['ec', 'e_c', 1., None], \
                ['dminSPH', 'd_min_actual', sum(structure['r']), None], ['ecSPH', 'ec_actual', structure['eb'], None], \
                ['dm', 'dm', 2., None], \
                ['E_d', 'e_dissipated', structure['eb'], 'small'], ['dE', 'e_conservation', 1., 'sci'], ['status', 'status', 1., 'string']]

        write_summary.write_system_parameters_formatted(fdir + 'system_results_paper.txt', final, params_to_write, paper = 1)

        #__write observable results
        params_to_write = [['ae', 'aei1', 1., None], ['ae', 'aei2', 1., None], ['P2/P1', 'period_ratio', 1., None], \
                ['mgas', 'mgas', 1., None], ['rho', 'd', 1., None], ['d1/d2', 'd_ratio', 1., None]]

        write_summary.write_system_parameters_formatted(fdir + 'system_observables_final_formatted.txt', final, params_to_write)
        write_summary.write_system_parameters_formatted(fdir + 'system_observables_initial_formatted.txt', initial, params_to_write)
        write_summary.write_system_parameters_formatted(fdir + 'system_observables_final_paper.txt', final, params_to_write, paper = 1)

        #__write initial conditions
        params_to_write = [['dmin', 'd_min', 1., None], ['dminSPH', 'd_min_actual', sum(structure['r']), None], \
                ['ec', 'e_c', 1., None], ['ecSPH', 'ec_actual', structure['eb'], None], \
                ['a', 'a', 1., None], ['e', 'e', 1., None]]

        write_summary.write_system_parameters_formatted(fdir + 'system_ic_formatted.txt', initial, params_to_write)
        write_summary.write_system_parameters_formatted(fdir + 'system_ic_paper.txt', initial, params_to_write, paper = 1)

        #write_summary.write_system_parameters_formatted_extra(fdir + 'system_parameters_final_formatted_extra.txt', [(k, run.final) for k, run in runs.iteritems()])




    def create_collision_dictionary(self, fname):

        f = open(fname, 'rb')

        # returns a tuple with values [eta, Ec, r1, r2]
        def get_values(g):
            return [[format(int(x[0]), "03"), [float(x[3]), float(x[1]), float(x[7]), float(x[8])]] for x in [line.split() for line in g] if len(x) >= 7 and x[0] != 'Run']

        parameters_dict = {k: v for k, v in get_values(f)}
        f.close()

        return parameters_dict




    def energy_calculations(self, fdir, run_sums, fname, structure):

        g   = 6.67390*math.pow(10,-8)                               # gravitational constant [cgs]
        rs  = 6.9599*math.pow(10,10)                                # unit of length [cm]
        ms  = 1.9891*math.pow(10,33)                                # unit of mass [g]
        ts  = math.sqrt(math.pow(rs,3)/(g*ms))                      # unit of time [s]

        #e_initial = {k: {key: run.energy_initial[k] for key, run in run_sums.iteritems()} for k in run_sums.values()[0].energy_initial.keys()}
        #e_final = {k: {key: run.energy_final[k] for key, run in run_sums.iteritems()} for k in run_sums.values()[0].energy_initial.keys()}

        e_conservation = {k: math.fabs((run.energy_initial['et'] - run.energy_final['et'])/run.energy_initial['et']) for k, run in run_sums.iteritems()}
        e_efficiency = {k: math.fabs((run.energy_final['ei'] - run.energy_initial['ei'])/ \
                (run.energy_final['ek'] + run.energy_final['eg'] - run.energy_initial['ek'] - run.energy_initial['eg'])) \
                for k, run in run_sums.iteritems()}
        e_thermal_increase_unitless = {k: (run.energy_final['ei'] - run.energy_initial['ei'])/run.energy_initial['ei'] for k, run in run_sums.iteritems()}
        e_thermal_increase_cgs = {k: (run.energy_final['ei'] - run.energy_initial['ei'])*ms*rs*rs/ts/ts for k, run in run_sums.iteritems()}
        e_dissipation = {k: run.energy_dissipated for k, run in run_sums.iteritems()}
        resolved = {k: run.resolved for k, run in run_sums.iteritems()}
        #e_dissipation = {k: (sum([p.gas_binding_energy for p in run.final['planets']]) - sum([p.gas_binding_energy for p in run.initial['planets']]) + \
        #        e_final['ei'][k] - e_initial['ei'][k])*ms*rs*rs/ts/ts/structure['eb'] for k, run in run_sums.iteritems()}

        f = open(fname, 'w')
        f.write('run\tde_tot\talpha_e\tde_i_%\tde_i_cgs\te_d\n')
        for k in run_sums:
            f.write(str(k) + '\t' + \
                    str(e_conservation[k]) + '\t' + \
                    str(e_efficiency[k]) + '\t' + \
                    str(e_thermal_increase_unitless[k]) + '\t' + \
                    str(e_thermal_increase_cgs[k]) + '\t' + \
                    str(e_dissipation[k]) + '\t' + \
                    str(resolved[k]) + '\n')
        f.close

        return e_conservation, e_efficiency, e_thermal_increase_unitless, e_thermal_increase_cgs




    def status_specific_analyses(self, fdir, runs, mstar, ref_system_name, structure, title):

        scatterings = {k: runs[k] for k in runs if 'stable' in runs[k].status}
        mergers = {k: runs[k] for k in runs if runs[k].status == 'merger'}
        planet_pairs = {k: runs[k] for k in runs if runs[k].status == 'bound'}

        status_specific_calcs.scattering_calcs(fdir, scatterings, mstar, ref_system_name, structure, title)
        status_specific_calcs.merger_calcs(fdir, mergers)
        status_specific_calcs.planet_pair_calcs(fdir, planet_pairs)




    def summarize_all_runs(self, fdir, run_sums, ref_system_name, structure, m_star, set_name):

        g  = 6.67390*math.pow(10,-8)                                # gravitational constant [cgs]
        rs = 6.9599*math.pow(10,10)                                 # unit of length [cm]
        ms = 1.9891*math.pow(10,33)                                 # unit of mass [g]
        ts = math.sqrt(math.pow(rs,3)/(g*ms))                       # unit of time [s]
        me = 5.972*10.**27                                          # earth mass [g]

        # edgecolors
        ec_dict  = {0: {'bound': 'k', 'merged': 'k', 'stable': 'k', 'unstable': 'k'},
                    1: {'bound': 'k', 'merged': 'k', 'stable': 'r', 'unstable': 'r'}}

        # facecolors
        fc_dict  = {'bound': 'k', 'merged': 'k', 'stable': 'y', 'unstable': 'none'}

        # markers
        m_dict = {'bound': '*', 'merged': 'x', 'stable': 'o', 'unstable': 'o'}

        collision_plots.exploratory_analysis(fdir, run_sums, structure, ec_dict, fc_dict, m_dict, set_name)
        collision_plots.m_transfer_plots(fdir, run_sums, structure, self.models, ec_dict, fc_dict, m_dict, set_name)
        collision_plots.observed_plots(fdir, run_sums, structure, ref_system_name, ec_dict, fc_dict, m_dict, set_name)
        collision_plots.de_fits(fdir, run_sums, structure, self.models, self.mass_dmin, ec_dict, fc_dict, m_dict, set_name)
        collision_plots.predict_outcome(fdir, run_sums, structure, self.models, ec_dict, fc_dict, m_dict, set_name)













