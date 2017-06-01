
import numpy as np

from reference_system import Ref_system
import summary_plots




def exploratory_analysis(fdir, runs, structure, ec_dict, fc_dict, m_dict, set_name):

    #status = [run.status for k, run in run_sums.iteritems()]
    ec = [ec_dict[runs[k].final['flipped']][runs[k].status] for k in runs]
    fc = [fc_dict[runs[k].status] for k in runs]
    markers = [m_dict[runs[k].status] for k in runs]

    # generate lists
    d_nbody = [runs[k].d_min for k in runs]
    e_nbody = [runs[k].e_c for k in runs]
    d_SPH = [runs[k].d_min_actual2['cubic']/sum(structure['r']) for k in runs]
    e_SPH = [runs[k].ec_actual/structure['eb'] for k in runs]

    # compare collision characteristics from nbody with SPH
    fit2 = [sorted(d_nbody), sorted(d_nbody)]
    summary_plots.all_collisions(d_nbody, d_SPH, ec, fc, markers, fdir, 'all_collisions_nbody_SPH_d', \
            r'$d_\mathrm{nbody}/(R_1+R_2)$', r'$d_\mathrm{SPH}/(R_1+R_2)$', fit = fit2, title = set_name)

    # ec, ec_SPH
    fit2 = [sorted(e_nbody), sorted(e_nbody)]
    summary_plots.all_collisions(e_nbody, e_SPH, ec, fc, markers, fdir, 'all_collisions_nbody_SPH_e', \
            r'$\left(E_\mathrm{c}/E_\mathrm{b}\right)_{nbody}$', r'$\left(E_\mathrm{c}/E_\mathrm{b}\right)_\mathrm{SPH}$', fit = fit2, title = set_name)

    # d_min, energy plot
    summary_plots.all_collisions(d_SPH, e_SPH, ec, fc, markers, fdir, 'all_collisions_ic_SPH', \
            r'$d/(R_1+R_2)$', r'$E_\mathrm{c}/E_\mathrm{b}$', ymax_in = max(e_SPH) + 0.005, title = set_name)

    # d_min, energy plot (nbody)
    summary_plots.all_collisions(d_nbody, e_nbody, ec, fc, markers, fdir, 'all_collisions_ic_nbody', \
            r'$d/(R_1+R_2)$', r'$E_\mathrm{c}/E_\mathrm{b}$', title = set_name)




def m_transfer_plots(fdir, runs, structure, models, ec_dict, fc_dict, m_dict, set_name):

    r1, r2 = structure['r']
    c1, c2 = structure['c_sph']
    d_core_hit = [(r1*c1+r2*c2)/(r1+r2)]

    d_scatter = [runs[k].d_min_actual2['cubic']/sum(structure['r']) for k in runs if runs[k].resolved and 'stable' in runs[k].status]
    e_scatter = [runs[k].ec_actual/structure['eb'] for k in runs if runs[k].resolved and 'stable' in runs[k].status]
    dm_scatter = [[-runs[k].final['dm'][ind] for k in runs if runs[k].resolved and 'stable' in runs[k].status] for ind in range(2)]
    dm_sum_scatter = [-sum(runs[k].final['dm']) for k in runs if 'stable' in runs[k].status]

    ec_scatter = [ec_dict[runs[k].final['flipped']][runs[k].status] for k in runs if 'stable' in runs[k].status]
    fc_scatter = [fc_dict[runs[k].status] for k in runs if 'stable' in runs[k].status]
    markers_scatter = [m_dict[runs[k].status] for k in runs if 'stable' in runs[k].status]


    # mass lost as a function of d_min
    f_pdm_d = [models['pdm_d'][ind]['f'] for ind in range(2)]
    for p_ind, (_f, _dm) in enumerate(zip(f_pdm_d, dm_scatter)):
        d_grid = np.linspace(min(d_scatter), 1., 10000)
        fit2 = [d_grid, [_f(_d, 0.) for _d in d_grid]]
        summary_plots.all_collisions(d_scatter, _dm, ec_scatter, fc_scatter, markers_scatter, fdir, 'all_collisions_d_dpm' + str(p_ind + 1), \
                r'$d/(R_1+R_2)$', r'$\Delta m/M_\mathrm{E}$', fit = fit2, ymin_in = -0.00001, title = set_name)

    # d_min, dm_all plot
    func_mtot_d = models['mtot_d']['f']
    d_grid = np.linspace(min(d_scatter), 1., 10000)
    fit2 = [d_grid, [func_mtot_d(_d, 0.) for _d in d_grid]]
    summary_plots.all_collisions(d_scatter, dm_sum_scatter, ec_scatter, fc_scatter, markers_scatter, fdir, 'all_collisions_dm_sum_d', \
            r'$d/(R_1+R_2)$', r'$\Delta m_\mathrm{tot}/M_\mathrm{E}$', fit = fit2, vline = d_core_hit, ymin_in = -0.00001, title = set_name)




def observed_plots(fdir, runs, structure, ref_system_name, ec_dict, fc_dict, m_dict, set_name):

    if ref_system_name: ref_system = Ref_system(ref_system_name)

    ec = [ec_dict[runs[k].final['flipped']][runs[k].status] for k in runs]
    fc = [fc_dict[runs[k].status] for k in runs]
    markers = [m_dict[runs[k].status] for k in runs]

    d_SPH = [runs[k].d_min_actual2['cubic']/sum(structure['r']) for k in runs if runs[k].resolved]
    e_SPH = [runs[k].ec_actual/structure['eb'] for k in runs if runs[k].resolved]

    # generate dictionaries holding final and initial summary parameters
    keys = runs.values()[0].initial.keys()
    final = {}
    initial = {}
    for k in keys:
        final[k] = [runs[key].final[k] for key in runs if runs[key].resolved]
        initial[k] = [runs[key].initial[k] for key in runs if runs[key].resolved]
        if isinstance(final[k], list):
            if isinstance(final[k][0], list): final[k] = map(list, zip(*final[k]))
        if isinstance(initial[k], list):
            if isinstance(initial[k][0], list): initial[k] = map(list, zip(*initial[k]))


    # mgas, period ratio plot
    for mgas_i, mgas_f, fname in zip(initial['mgas'], final['mgas'], ['all_collisions_mgas_b_ratio', 'all_collisions_mgas_c_ratio']):

        arg0 = [mgas_i, initial['period_ratio']]
        arg1 = [['yline','r--'], ['xline','r--']]

        if ref_system_name:
            arg0 = [mgas_i, initial['period_ratio'], ref_system.pratio]
            arg1 = [['yline','r--'], ['xline','r--'], ['xline','k--']]

        summary_plots.all_collisions(mgas_f, final['period_ratio'], ec, fc, markers, fdir, fname, \
                r'$m_{\mathrm{gas},b}$', r'$P_c/P_b$', arg0, arg1, title = set_name)

    for i, xlabel, fname in zip(range(2), [r'$\rho_\mathrm{b}$', r'$\rho_\mathrm{c}$'], ['all_collisions_rhob_ratio', 'all_collisions_rhoc_ratio']):

        arg0 = [[initial['d'][i][0], initial['period_ratio'][0]], [initial['d'][i][0], initial['period_ratio'][0]]]
        arg1 = [['xyline','k--'], ['point', 'k']]
        if ref_system_name:
            arg0 = [[ref_system.d[i], ref_system.pratio], [initial['d'][i][0], initial['period_ratio'][0]], [initial['d'][i][0], initial['period_ratio'][0]]]
            arg1 = [['point','b'], ['xyline','k--'], ['point', 'k']]
        # rho_b, period ratio plot
        summary_plots.all_collisions(final['d'][i], final['period_ratio'], ec, fc, markers, fdir, fname, \
                xlabel, r'$P_c/P_b$', arg0, arg1, title = set_name)


    # d_min, period ratio
    summary_plots.all_collisions(d_SPH, final['period_ratio'], ec, fc, markers, fdir, 'all_collisions_d_p_ratio', \
            r'$d/(R_1+R_2)$', r'$P_2/P_1$', title = set_name)

    # ec, period ratio
    summary_plots.all_collisions(e_SPH, final['period_ratio'], ec, fc, markers, fdir, 'all_collisions_e_p_ratio', \
            r'$E_\mathrm{c}/E_\mathrm{b}$', r'$P_2/P_1$', title = set_name)


    # density ratio period ratio plot
    arg0 = [[initial['d_ratio'][0], initial['period_ratio'][0]], [initial['d_ratio'][0], initial['period_ratio'][0]]]
    arg1 = [['xyline','k--'], ['point', 'k']]
    if ref_system_name:
        arg0 = [[ref_system.rho_ratio, ref_system.pratio], [initial['d_ratio'][0], initial['period_ratio'][0]], [initial['d_ratio'][0], initial['period_ratio'][0]]]
        arg1 = [['point','b'], ['xyline','k--'], ['point', 'k']]
    summary_plots.all_collisions(final['d_ratio'], final['period_ratio'], ec, fc, markers, fdir, 'all_collisions_rho_ratio', \
            r'$\rho_\mathrm{b}/\rho_\mathrm{c}$', r'$P_c/P_b$', \
            arg0, arg1, title = set_name)


    # mass lost as a function of both energy and distance of closest approach
    for mgas_i, mgas_f, fname, fname2, ylabel in zip(initial['mgas'], final['mgas'], \
            ['all_collisions_mgas_b_d', 'all_collisions_mgas_c_d'], ['all_collisions_mgas_b_e', 'all_collisions_mgas_c_e'], \
            [r'$\Delta m_\mathrm{b}$', r'$\Delta m_\mathrm{c}$']):

        dmgas = [_mgas_i - _mgas_f for _mgas_i, _mgas_f in zip(mgas_i, mgas_f)]

        summary_plots.all_collisions(d_SPH, dmgas, ec, fc, markers, fdir, fname, \
                r'$d_\mathrm{min}$', ylabel, ymin_in = -0.00001, title = set_name)

        summary_plots.all_collisions(e_SPH, dmgas, ec, fc, markers, fdir, fname2, \
                r'$E_\mathrm{c}/E_\mathrm{b}$', ylabel, ymin_in = -0.00001, title = set_name)




def de_fits(fdir, runs, structure, models, mass_dmin, ec_dict, fc_dict, m_dict, set_name):

    me = 5.972*10.**27                                          # earth mass [g]

    # dictionaries for resolved runs
    ec = [ec_dict[runs[k].final['flipped']][runs[k].status] for k in runs if runs[k]]
    fc = [fc_dict[runs[k].status] for k in runs if runs[k].resolved]
    markers = [m_dict[runs[k].status] for k in runs if runs[k].resolved]

    m_outside = [mass_dmin(runs[k].d_min_actual2['cubic'])/me for k in runs if runs[k].resolved]
    de = [runs[k].energy_dissipated/structure['eb'] for k in runs if runs[k].resolved]
    d_SPH = [runs[k].d_min_actual2['cubic']/sum(structure['r']) for k in runs if runs[k].resolved]
    e_SPH = [runs[k].ec_actual/structure['eb'] for k in runs if runs[k].resolved]

    r1, r2 = structure['r']
    c1, c2 = structure['c_sph']
    d_core_hit = [(r1*c1+r2*c2)/(r1+r2)]
    d_grid = np.linspace(min(d_SPH), 1., 10000)

    #__Exploratory
    # d_min, energy dissipated
    summary_plots.all_collisions(d_SPH, de, ec, fc, markers, fdir, 'all_collisions_de_d', \
            r'$d/(R_1+R_2)$', r'$\Delta E_\mathrm{orb}/E_\mathrm{b}$', ymin_in = -0.00001, title = set_name)

    # ec, energy dissipated
    summary_plots.all_collisions(e_SPH, de, ec, fc, markers, fdir, 'all_collisions_de_e', \
            r'$E_\mathrm{c}/E_\mathrm{b}$', r'$\Delta E_\mathrm{orb}/E_\mathrm{b}$', ymin_in = -0.00001, title = set_name)

    #__Fit for energy dissipation
    # m_out, energy dissipated plot
    func = models['de_m']['f']
    m_grid = np.linspace(min(m_outside), max(m_outside), 10000)
    fit2 = [m_grid, [func(_m, 0.) for _m in m_grid]]
    summary_plots.all_collisions(m_outside, de, ec, fc, markers, fdir, 'all_collisions_de_m', \
            r'$m_\mathrm{out}/M_\mathrm{E}$', r'$\Delta E_\mathrm{orb}/E_\mathrm{b}$', xmin_in = -0.00001, ymin_in = -0.00001, title = set_name)

    # energy lost as a function of d_min and e
    func = models['de_d']['f']
    fit2 = [d_grid, [func(_d, 0.) for _d in d_grid]]
    summary_plots.all_collisions(d_SPH, de, ec, fc, markers, fdir, 'all_collisions_d_de', \
            r'$d/(R_1+R_2)$', r'$\Delta E_\mathrm{orb}/E_\mathrm{b}$', fit = fit2, ymin_in = -0.00001, vline = d_core_hit, title = set_name)




def predict_outcome(fdir, runs, structure, models, ec_dict, fc_dict, m_dict, set_name):

    r1, r2 = structure['r']
    c1, c2 = structure['c_sph']
    d_core_hit = [(r1*c1+r2*c2)/(r1+r2)]

    ec = [ec_dict[runs[k].final['flipped']][runs[k].status] for k in runs]
    fc = [fc_dict[runs[k].status] for k in runs]
    markers = [m_dict[runs[k].status] for k in runs]

    d_SPH = [runs[k].d_min_actual2['cubic']/sum(structure['r']) for k in runs]
    e_SPH = [runs[k].ec_actual/structure['eb'] for k in runs]
    de = [runs[k].energy_dissipated/structure['eb'] for k in runs]
    e_escape = [[_e/structure['eb'] for _e in runs[k].e_escape] for k in runs]

    func = models['de_d']['f']
    de_fit = [func(_d, _e) for _d, _e in zip(d_SPH, e_SPH)]

    fit2 = [np.linspace(min(de), max(de), 10000), np.linspace(min(de), max(de), 10000)]
    summary_plots.all_collisions(de, de_fit, ec, fc, markers, fdir, 'all_collisions_fit', \
            r'$\Delta E$', r'$\Delta E_\mathrm{fit}$', fit = fit2, ymin_in = -0.00001, title = set_name)

    # d_min, energy plot with fit (initial e_escape)
    d_grid = np.linspace(min(d_SPH), 1., 10000)
    fit2 = [d_grid, [func(_d, 0.) for _d in d_grid]]
    e0 = [_e - _e_escape[0] for _e, _e_escape in zip(e_SPH, e_escape)]
    summary_plots.all_collisions(d_SPH, e0, ec, fc, markers, fdir, 'all_collisions_ic_predict', \
            r'$d/(R_1+R_2)$', r'$(E_\mathrm{c}-E_\mathrm{esc})/E_\mathrm{b}$', fit = fit2, ymax_in = max(e0) + 0.005, vline = d_core_hit, title = set_name)

    # d_min, energy plot with fit (final e_escape)
    e0 = [_e - _e_escape[1] for _e, _e_escape in zip(e_SPH, e_escape)]
    summary_plots.all_collisions(d_SPH, e0, ec, fc, markers, fdir, 'all_collisions_ic_predict2', \
            r'$d/(R_1+R_2)$', r'$(E_\mathrm{c}-E_\mathrm{esc})/E_\mathrm{b}$', fit = fit2, ymax_in = max(e0) + 0.005, vline = d_core_hit, title = set_name)







