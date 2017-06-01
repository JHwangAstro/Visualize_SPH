
import math
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as op

import orbit_calc


def all_collisions(x, y, ec, fc, markers, fdir, f, \
        xlabel = r'$x$', ylabel = r'$y$', \
        compare = None, compare_style = None, \
        xmin_in = None, xmax_in = None, \
        ymin_in = None, ymax_in = None, \
        fit = None, vline = None, \
        title = None):

    def update_lim(vmin, vmax, new_v):

        new_vmax = np.max(new_v)*1.1**np.sign(np.max(new_v))
        if vmax: new_vmax = max([vmax, new_vmax])

        new_vmin = np.min(new_v)/1.1**np.sign(np.min(new_v))
        if vmin: new_vmin = min([vmin, new_vmin])

        return new_vmin, new_vmax


    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    # first plot the points with a facecolor
    for _x, _y, _ec, _fc, _markers in zip(x, y, ec, fc, markers):
        if _fc != 'none': ax1.scatter(_x, _y, marker = _markers, edgecolors = _ec, facecolors = _fc, s = 80)
    for _x, _y, _ec, _fc, _markers in zip(x, y, ec, fc, markers):
        if _fc == 'none': ax1.scatter(_x, _y, marker = _markers, edgecolors = _ec, facecolors = _fc, s = 80)

    xmin, xmax = update_lim(xmin_in, xmax_in, x)
    ymin, ymax = update_lim(ymin_in, ymax_in, y)

    if fit:
        if isinstance(fit[0][0], list):
            for _fit in fit: ax1.plot(_fit[0], _fit[1], linestyle = '--', color = 'k')
        else:
            ax1.plot(fit[0], fit[1], linestyle = '--', color = 'k')

    if vline:
        for _vline in vline:
            ax1.plot([_vline, _vline], [ymin, ymax], linestyle = ':', c = 'k')

    if title: ax1.set_title(title, loc = 'left')

    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])

    ax1.set_xlabel(xlabel, fontsize=15)
    ax1.set_ylabel(ylabel, fontsize=15)

    fname = fdir + f + '.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + f + '.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)


    if compare and compare_style:

        # update limits first
        for v, c_style in zip(compare, compare_style):

            code, style = c_style
            if code == 'yline': xmin, xmax = update_lim(xmin, xmax, v)
            if code == 'xline': ymin, ymax = update_lim(ymin, ymax, v)
            if code == 'xyline' or code == 'point':
                xmin, xmax = update_lim(xmin, xmax, v[0])
                ymin, ymax = update_lim(ymin, ymax, v[1])

        # make comparison plots
        for v, c_style in zip(compare, compare_style):

            code, style = c_style
            if code == 'yline': ax1.plot([xmin, xmax], [v, v], style)
            if code == 'xline': ax1.plot([v, v], [ymin, ymax], style)
            if code == 'point':
                if isinstance(v[0], list):
                    for _v in v:
                       ax1.scatter(_v[0], _v[1], color = style, s = 80)
                else:
                    ax1.scatter(v[0], v[1], color = style, s = 80)
            if code == 'xyline':
                for _x, _y, in zip(x, y):
                    if isinstance(v[0], list):
                        for _v in v:
                            ax1.plot([_x, _v[0]], [_y, _v[1]], style, alpha = 0.5)
                    else:
                        ax1.plot([_x, v[0]], [_y, v[1]], style, alpha = 0.5)

        ax1.set_xlim([xmin, xmax])
        ax1.set_ylim([ymin, ymax])

        fname = fdir + f + '_compare.png'
        plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
        fname = fdir + f + '_compare.ps'
        plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()




# draws AMD, period-ratio lines for all runs
def stability_plot(fdir, f, run_sums, runs_i, runs_f, mstar, C0, structure, ref_system = None, title = None):

    g = 6.67390*10.**-8                               # gravitational constant [cgs]
    ms = 1.9891*math.pow(10,33)                                # unit of mass [g]
    me = 5.972*10.**27
    au = 1.496*10.**13.

    C = ms**1.5*au**.5

    def h(r):
        return sum([mstar*r['m'][i]/(mstar + r['m'][i])*(g*(mstar + r['m'][i])*r['a'][i]*(1.-r['e'][i]**2.))**.5*math.cos(r['i'][i]*math.pi/180.) for i in range(2)])*C

    # angular momentum at 0 ecc, inner a, and pratio
    def h2(r, p):
        energy = -g/2.*sum([(mstar+_m)*_m/_a for _m, _a in zip(r['m'], r['a'])])*ms
        m1, m2 = r['m']
        a1 = -g/2./energy*((mstar+m1)*m1+(mstar+m2)*m2/p**(2./3.))*ms
        a2 = p**(2./3.)*a1
        return sum([mstar*r['m'][i]/(mstar + r['m'][i])*(g*(mstar + r['m'][i])*a)**.5 for i, a in zip(range(2), [a1, a2])])*C

    def h_elements(m, a, e, i):
        return sum([mstar*_m/(mstar + _m)*(g*(mstar + _m)*_a*(1.-_e**2.))**.5*math.cos(_i*math.pi/180.) for _m, _a, _e, _i in zip(m, a, e, i)])*C

    # angular momentum at 0 ecc, inner a, and pratio
    def h_elements_2(m, a, p):
        energy = -g/2.*sum([(mstar+_m)*_m/_a for _m, _a in zip(m, a)])*ms
        m1, m2 = m
        a1 = -g/2./energy*((mstar+m1)*m1+(mstar+m2)*m2/p**(2./3.))*ms
        a2 = p**(2./3.)*a1
        return sum([mstar*_m/(mstar + _m)*(g*(mstar + _m)*_a)**.5 for _m, _a in zip(m, [a1, a2])])*C


    # plot lines for each system
    max_AMD = 0.
    AMD_min = [[],[]]
    for runs, _c, _AMD_min in zip([runs_i, runs_f], ['g', 'r'], AMD_min):
        p_range = orbit_calc.find_p_ratios(runs, mstar)
        for run, _prange in zip(runs, p_range):
            # create list of possible pratios
            h_run = h(run)
            pratio_grid = np.linspace(_prange[0], _prange[1], 100)
            # for each pratio, calculate the AMD
            AMD_grid = [(h2(run, _p) - h_run)/C0 for _p in pratio_grid]
            max_AMD = max([max_AMD, max(AMD_grid)])
            plt.plot(AMD_grid, pratio_grid, linestyle = '-', color = _c, alpha = 0.5, zorder = 1)

            def AMD(x):
                if x < _prange[0]: return 1000.
                if x > _prange[1]: return 1000.
                return (h2(run, x) - h_run)/C0

            res = op.minimize(AMD, 1.)

            _AMD_min.append(AMD(res['x']))

    # plot coordinates of systems
    AMD_i = [orbit_calc.calculate_AMD_run(mstar, run)/C0 for run in runs_i]
    AMD_f = [orbit_calc.calculate_AMD_run(mstar, run)/C0 for run in runs_f]

    pratio_i = [(run['a'][1]/run['a'][0])**1.5 for run in runs_i]
    pratio_f = [(run['a'][1]/run['a'][0])**1.5 for run in runs_f]

    plt.scatter(AMD_i, pratio_i, color = 'g', edgecolor = 'k', s = 2., zorder = 5)
    plt.scatter(AMD_f, pratio_f, color = 'r', edgecolor = 'k', s = 2., zorder = 5)

    # plot comparison system
    if ref_system:
        p_range = orbit_calc.find_p_ratios2(ref_system.m, ref_system.a, ref_system.e, mstar)
        h_run = h_elements(ref_system.m, ref_system.a, ref_system.e, [0., 0.])

        pratio_grid = np.linspace(p_range[0], p_range[1], 100)

        # for each pratio, calculate the AMD
        AMD_grid = np.array([(h_elements_2(ref_system.m, ref_system.a, _p) - h_run)/C0 for _p in pratio_grid])

        plt.plot(AMD_grid[AMD_grid < max_AMD], pratio_grid[AMD_grid < max_AMD], linestyle = '-', color = 'b', zorder = 1)

        a1, a2 = ref_system.a
        pratio_ref = (a2/a1)**1.5
        AMD_ref = (h_elements_2(ref_system.m, ref_system.a, pratio_ref)-h_run)/C0

        plt.scatter(AMD_ref, pratio_ref, color = 'b', edgecolor = 'k', s = 2., zorder = 5)

        p_range = orbit_calc.find_p_ratios2(ref_system.m, ref_system.a, [0.01,0.01], mstar)
        h_run = h_elements(ref_system.m, ref_system.a, [0.01,0.], [0., 0.01])
        pratio_grid = np.linspace(p_range[0], p_range[1], 100)
        # for each pratio, calculate the AMD
        AMD_grid = np.array([(h_elements_2(ref_system.m, ref_system.a, _p) - h_run)/C0 for _p in pratio_grid])
        plt.plot(AMD_grid[AMD_grid < max_AMD], pratio_grid[AMD_grid < max_AMD], linestyle = '-', color = 'b', zorder = 1)

    plt.axvline(x=0., color = 'k', linestyle = '--')

    # set xlim
    plt.axis([min([0., min(AMD_i), min(AMD_f)])-100, max([max(AMD_i), max(AMD_f)])+100, \
            min([min(pratio_i), min(pratio_f)]) - 0.1, max([max(pratio_i), max(pratio_f)]) + 0.1])

    if title: plt.title(title, loc = 'left')
    plt.ylabel(r'$P_2/P_1$')
    plt.xlabel(r'$C/C_0$')

    #plt.xlim([-200,600])

    fname = fdir + '/' + f + '.png'
    plt.savefig(fname, facecolor = 'w', edgecolor = 'w', format = 'png', dpi = 200)
    fname = fdir + '/' + f + '.ps'
    plt.savefig(fname, facecolor = 'w', edgecolor = 'w', format = 'ps', dpi = 200)

    # plot stable
    delta = 2.*3.**.5
    for _a in [0.1]:

        R_Hill = 2.*_a/(2./(sum(runs_i[0]['m'])/3./mstar)**(1./3.)-delta)
        _a2 = _a + delta*R_Hill
        p_range = orbit_calc.find_p_ratios2(runs_i[0]['m'], [_a, _a2], [0., 0.], mstar)

        h_run = h_elements(runs_i[0]['m'], [_a, _a2], [0., 0.], [0., 0.])

        pratio_grid = np.linspace(p_range[0], p_range[1], 100)

        # for each pratio, calculate the AMD
        AMD_grid = np.array([(h_elements_2(runs_i[0]['m'], [_a, _a2], _p) - h_run)/C0 for _p in pratio_grid])

        plt.plot(AMD_grid[AMD_grid < max_AMD], pratio_grid[AMD_grid < max_AMD], linestyle = '-', color = 'k', zorder = 1)

    fname = fdir + '/' + f + '_Hill_Stable.png'
    plt.savefig(fname, facecolor = 'w', edgecolor = 'w', format = 'png', dpi = 200)
    fname = fdir + '/' + f + '_Hill_Stable.ps'
    plt.savefig(fname, facecolor = 'w', edgecolor = 'w', format = 'ps', dpi = 200)

    plt.clf()


    AMD_min_i, AMD_min_f = AMD_min
    d_SPH = [run_sums[k].d_min_actual2['cubic']/sum(structure['r']) for k in run_sums]
    d_AMD = [_f - _i for _i, _f in zip(AMD_min_i, AMD_min_f)]
    plt.scatter(d_SPH, d_AMD)

    fname = fdir + '/' + f + '_d_AMD.png'
    plt.savefig(fname, facecolor = 'w', edgecolor = 'w', format = 'png', dpi = 200)
    fname = fdir + '/' + f + '_d_AMD.ps'
    plt.savefig(fname, facecolor = 'w', edgecolor = 'w', format = 'ps', dpi = 200)

    plt.clf()






def all_collision_ellipse_plot():

    fig = plt.figure(0)
    ax = fig.add_subplot(111, aspect='equal')
    for _a1, _a2, _e1, _e2, _fc, _ec in zip(a1, a2, e1, e2, fc, ec):
        ell = Ellipse(xy = np.array([_a1, _a2]), width = 2.*_a1*_e1, height = 2.*_a2*_e2)
        ax.add_artist(ell)
        ell.set_clip_box(ax.bbox)
        ell.set_alpha(0.5)
        ell.set_facecolor(_fc)
        ell.set_edgecolor(_ec)

    if ref_system_name:
        ell = Ellipse(xy = np.array(rs.a), width = 2.*rs.a[0]*rs.e[0], height = 2.*rs.a[0]*rs.e[0])
        ax.add_artist(ell)
        ell.set_clip_box(ax.bbox)
        ell.set_alpha(0.5)
        ell.set_facecolor('b')
        ell.set_edgecolor('b')

    plt.xlim([0, max(np.array(a1)*(1.+np.array(e1)))])
    plt.ylim([0, max(np.array(a2)*(1.+np.array(e2)))])

    if ref_system_name:
        plt.xlim([0, max([rs.a[0]*(1.+rs.e[0]), max(np.array(a1)*(1.+np.array(e1)))])])
        plt.ylim([0, max([rs.a[1]*(1.+rs.e[1]), max(np.array(a2)*(1.+np.array(e2)))])])

    plt.xlabel(r'$a_1$', fontsize=15)
    plt.ylabel(r'$a_2$', fontsize=15)

    fname = fdir + 'all_collisions_a1_a2.png'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    fname = fdir + 'all_collisions_a1_a2.ps'
    plt.savefig(fname, facecolor='w', edgecolor='w', format='ps',dpi=200)

    plt.clf()




def histogram_1d(x, fdir, fname, xlabel, ylabel, vline = None, n_bins = 6, c = 'red'):

    plt.hist(x, histtype = 'stepfilled', bins = n_bins, edgecolor = c, linewidth = 2, facecolor='none')
    _x, _y, _ = plt.hist(x, histtype = 'stepfilled', bins = n_bins, edgecolor = c, linewidth = 2, facecolor='none')

    if vline:
        for _vline in vline:
            v, v_color = _vline
            plt.plot([v, v], [0., max(_x) + 1], linestyle = ':', c = v_color)

    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    plt.savefig(fdir + '/' + fname + '.png', facecolor='w', edgecolor='w', format='png', dpi=200)
    plt.savefig(fdir + '/' + fname + '.ps', facecolor='w', edgecolor='w', format='ps', dpi=200)
    plt.clf()






def histogram_2d(x, fdir, fname, xlabel, ylabel, vline = None, n_bins = 10, c = ['green', 'red'], alpha = 0.5):

    max_n = 0

    bins = np.linspace(np.min(x), np.max(x), n_bins)

    for _x, _c in zip(x, c):
        #plt.hist(_x, histtype = 'stepfilled', bins = n_bins, edgecolor = _c, linewidth = 2, facecolor='none')
        _n, _y, _ = plt.hist(_x, histtype = 'stepfilled', bins = bins, alpha = alpha, edgecolor = _c, linewidth = 2, facecolor = 'none')
        max_n = max([max_n, np.max(_n)])

    if vline:
        for _vline in vline:
            v, v_color = _vline
            plt.plot([v, v], [0., max_n + 1], linestyle = ':', c = v_color)

    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    plt.yticks(range(0, int(max_n) + 1))
    plt.savefig(fdir + '/' + fname + '.png', facecolor='w', edgecolor='w', format='png', dpi=200)
    plt.savefig(fdir + '/' + fname + '.ps', facecolor='w', edgecolor='w', format='ps', dpi=200)
    plt.clf()






