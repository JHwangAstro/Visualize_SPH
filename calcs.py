
import math
import matplotlib
matplotlib.use('Agg')
import corner
import emcee
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import scipy.optimize as op

def get_mass_dmin(profiles):

    # get mass(radius) relationship
    mass_radius = []
    for profile in profiles:

        m, r, p, d, t, mu, u = np.loadtxt(profile, unpack = True)

        # get gas envelope mass radius profile
        mass_radius.append({_r: _m for _r, _m, _d in zip(r, m, d) if _d < 1.})

    # convert to mass(d_min) for both planets
    max_m = [max(mr.values()) for mr in mass_radius]
    r_from_surface = [[max(mr) - _r for _r in sorted(mr)] for mr in mass_radius]
    m_outside = [[max(mr.values()) - mr[_r] for _r in sorted(mr)] for mr in mass_radius]
    mr_interp = [interp1d(_r, _m, kind = 'cubic') for _r, _m in zip(r_from_surface, m_outside)]

    # d must be in cm
    def mass_dmin(d):

        r_in = sum([max(mr) for mr in mass_radius]) - d

        if r_in < 0.: return 0.

        mass = 0.
        for mr, _m_outside, f in zip(mass_radius, m_outside, mr_interp):

            if r_in/2. > max(mr) - min(mr):
                mass = mass + max(_m_outside)
            else:
                mass = mass + f(r_in/2.)

        return mass

    return mass_dmin




def fit_2d(y, x1, x2, f, x_guess, keys):

    theta = curve_fit(f, [x1, x2], y, x_guess, maxfev = 1000000)[0]

    print 'best fit:\n',
    for k, _theta in zip(keys, theta):
        print k, ': ', _theta
    print '\n'

    def best_fit_f(x1, x2):

        return f([x1, x2], *theta)

    return {'f': best_fit_f, 'theta': {key: _theta for key, _theta in zip(keys, theta)}}




def get_posteriors(f, x, y, theta):

    big = 10.**30.
    ndim = len(theta)
    nwalkers = 4*ndim # must be even

    # log likelihood
    def ln_like(theta, x, y):

        #a, b, c, d0, z = theta
        #model = [f(_x, a, b, c, d0, z) for _x in x]
        #a, b, c, gam = theta
        #model = [f(_x, a, b, c) for _x in x]
        model = [f(_x, *theta) for _x in x]
        return -.5*np.sum([(_y - _model)**2. for _y, _model in zip(y, model)])
        #return -.5*np.sum([math.fabs(_y - _model) for _y, _model in zip(y, model)])


    # we are assuming uniform for now
    def ln_prior(theta):

        #a, b, c, d0, z = theta
        #if -big < a < big and 0. < d0 < 1. and 0. < b < big and -big < c < big and -big < z < big:
        #    return 0.0
        a, b, c, gam, d0 = theta
        if -big < a < 0. and -big < b < 0. and -big < c < big and d0 < min([_x[1] for _x in x]):
            return 0.0
        return -np.inf


    def ln_post(theta, x, y):

        if not np.isfinite(ln_prior(theta)): return -np.inf

        return ln_like(theta, x, y) + ln_prior(theta)


    nll = lambda *args: -ln_like(*args)
    result = op.minimize(nll, [theta], args = (x, y))


    #keys = ['A', 'B', 'C', 'd0', 'z']
    keys = ['A', 'B', 'C', r'$\gamma$', r'$d_0$']

    # add in some noise
    pos = [result['x'] + 1.e-7*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_post, args = (x, y))
    sampler.run_mcmc(pos, 5000)
    samples = sampler.chain[:, 500:, :].reshape((-1, ndim))

    # plot stuff
    fig = corner.corner(samples, labels = keys, truths = result['x'])
    fig.savefig("./triangle.png")

    return samples




# returns the energy required to sustain an orbit at the Hill Radius
def find_e_escape(mi, mf, a, mstar):

    g  = 6.67390*10**-8
    au = 1.496*10**13
    ms = 1.9891*10**33

    # find the largest Hill radius between the two planets
    r_Hill = max([(_m/mstar/3.)**(1./3.)*sum(a)/2. for _m in mi])
    r_Hill2 = max([(_m/mstar/3.)**(1./3.)*sum(a)/2. for _m in mf])

    # Calculate the energy required to escape the Hill radius,
    # and normalize to binding energy
    m0, m1 = mi
    m0f, m1f = mf
    #return -g*m0*m1/(1.5*r_Hill)*ms*ms/au
    return [-g*m0f*m1f/(r_Hill)*ms*ms/au, -g*m0f*m1f/r_Hill2/2.*ms*ms/au]




def fit_d_dissipation_old_2(d, ec, e, d_all):

    def func(x, a, b, c, gam, d0):

        _d, _ec = x

        b = 0.

        if (d0 > min(d_all)) or (d0 < 0.): return 1e38

#        model = a*(_d-d0)**gam + b*_ec + c
        model = a*_d**gam + b*_ec + c

        return model


    coef = curve_fit(func, [d, ec], e, [.01, 0., 0., -2., 0.0], maxfev = 1000000)[0]
    #coef = curve_fit(func, [d, ec], e, [-.1, -1., 1., .2, .45], maxfev = 1000000)[0]
    print 'best fit: ', coef

#    a, b, c, d0, z = coef

    return func, coef




def fit_d_dissipation_old(d, ec, e, d_all):

    def func(x, a, b, c, d0, z):

        _d, _ec = x

        # enforce boundaries
        #c = 0.
        if (a < 0.) or (b < 0.) or (d0 > min(d_all)) or (d0 < 0.): return 1e38
#        if (a < 0.) or (b < 0.) or (d0 < 0.) or (d0 > 1.): return 1e38
        #if (a > 0.) or (z > 0.) : return 1e38
#        if (d0 > min(d_all)): return 1e38

        model = a/(_d-d0)**b + c + z*_ec
        #model = a*np.exp(-b*(_d-d0)) + c + z*_ec
#        model = a*_d + z*_ec

        return model


    coef = curve_fit(func, [d, ec], e, [.001, 2., 0., .5, -4.], maxfev = 1000000)[0]
    print 'best fit: ', coef

    a, b, c, d0, z = coef

    return func, z, d0, coef







