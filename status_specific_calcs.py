
import os

import orbit_calc
from reference_system import Ref_system
import summary_plots


def scattering_calcs(fdir, runs, mstar, ref_system_name, structure, title):

    if ref_system_name:
        ref_system = Ref_system(ref_system_name)
    else:
        ref_system = None

    g = 6.67390*10.**-8                               # gravitational constant [cgs]
    ms = 1.9891*10.**33                                # unit of mass [g]
    au = 1.496*10.**13.
    me = 5.9722*10.**27.

    scattering_dir = fdir + 'Scattering'
    if not os.path.exists(scattering_dir): os.makedirs(scattering_dir)

    # calculate initial and final AMD for all runs
    AMD_0 = ms**.5*me*mstar**.5*g**.5*au**.5*(sum(runs[min(runs)].initial['m']))
    AMD_i = [orbit_calc.calculate_AMD(mstar, runs[k].initial['m'], runs[k].initial['aei1'], runs[k].initial['aei2'], runs[k].initial['i'])/AMD_0 for k in runs]
    AMD_f = [orbit_calc.calculate_AMD(mstar, runs[k].final['m'], runs[k].final['aei1'], runs[k].final['aei2'], runs[k].final['i'])/AMD_0 for k in runs]
    AMD_norm = [(_f-_i) for _i, _f in zip(AMD_i, AMD_f)]

    # make delta AMD histograms
    _vline = None
    if ref_system: _vline = [[ref_system.AMD/(ms**.5*me*mstar**.5*g**.5*au**.5*(sum(runs[min(runs)].initial['m']))), 'g']]
    summary_plots.histogram_2d([AMD_i, AMD_f], scattering_dir, 'scatterings_AMD', r'$C/C_0$', r'$N$', _vline)

    # density ratio histogram
    density = [runs[k].final['d_ratio'] for k in runs]
    _vline = None
    if ref_system: _vline = [[ref_system.rho_ratio, 'g'], [runs[min(runs)].initial['d_ratio'], 'k']]
    summary_plots.histogram_1d(density, scattering_dir, 'scatterings_rho_ratio', r'$\rho_1/\rho_2$', r'$N$', _vline, n_bins = 10)

    # period ratio histogram
    pratio_i = [runs[k].initial['period_ratio'] for k in runs]
    pratio_f = [runs[k].final['period_ratio'] for k in runs]
    _vline = None
    if ref_system: _vline = [[ref_system.pratio, 'g']]
    summary_plots.histogram_1d(pratio_f, scattering_dir, 'scatterings_p_ratio', r'$P_2/P_1$', r'$N$', _vline, n_bins = 10)

    # AMD and period ratio scatter plot
    runs_i = [runs[k].initial for k in runs]
    runs_f = [runs[k].final for k in runs]
    summary_plots.stability_plot(scattering_dir, 'stability', runs, runs_i, runs_f, mstar, AMD_0, structure, None, title)


    # time to instability histogram




def merger_calcs(fdir, runs):

    pass




def planet_pair_calcs(fdir, runs):

    pass


