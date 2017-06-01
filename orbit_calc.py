
import math
import numpy as np
from numpy import linalg as LA
import scipy.optimize as op
from scipy.interpolate import interp1d

import mercury_routines


def get_bound_planet_aei(snapshots, dcrit1, dcrit2):

    g   = 6.67390*math.pow(10,-8)                               # gravitational constant [cgs]
    rs  = 6.9599*math.pow(10,10)                                # unit of length [cm]
    ms  = 1.9891*math.pow(10,33)                                # unit of mass [g]
    ts  = math.sqrt(math.pow(rs,3)/(g*ms))                      # unit of time [s]
    au = 1.496*10.**13.

    def get_p_sep(s):

        p1, p2 = s.planets
        return LA.norm([p1.r-p2.r])

    t = np.array([s.unitless_time for s in snapshots])
    p_sep = np.array([get_p_sep(s)*rs for s in snapshots])
    r_t = {s.unitless_time: [LA.norm(p.r) for p in s.planets] for s in snapshots}

    # get distance of first collision
    t_0, t_1 = get_collision_times(t, p_sep, dcrit1*au, dcrit2)[0]

    d_t = {_d/rs: _t for _d, _t in zip(p_sep, t) if t_0<_t<t_1}
    if len(d_t) > 0:
        t_dmin = d_t[min(d_t)]
        r_dmin = r_t[t_dmin]


    bound_orbits = []



    for snapshot in snapshots:

        p1, p2 = snapshot.planets
        m_star = snapshot.star.m

        mu = p1.m*p2.m/(p1.m+p2.m)

        # attempt to normalize the effect of the host star; set these to 1. otherwise
        #r_cm = LA.norm((p1.r*p1.m+p2.r*p2.m)/(p1.m+p2.m))
        #v1_factor = (1.-m_star*(-1./LA.norm(p1.r)+1./LA.norm(r_dmin[0]))/LA.norm(p1.v)**2.)**.5
        #v2_factor = (1.-m_star*(-1./LA.norm(p2.r)+1./LA.norm(r_dmin[1]))/LA.norm(p2.v)**2.)**.5
        #v1_factor = 1.
        #v2_factor = 1.

        #v = LA.norm([p1.v*v1_factor - p2.v*v2_factor])
        v = LA.norm([p1.v - p2.v])
        r = LA.norm([p1.r - p2.r])

        e_kin = 0.5*mu*v**2.*ms*rs*rs/ts/ts
        e_pot = -g*p1.m*p2.m/r/rs*ms*ms
        # correcting for the host star's energy
        #e_correct = g*p1.m*m_star*(1./LA.norm(p1.r)+1./LA.norm(p2.r))*ms*ms/rs
        #print 'correction: ', LA.norm(p1.r), LA.norm(p2.r), r_dmin
        #e_correct = g*m_star*(-p1.m/LA.norm(p1.r)+p1.m/r_dmin[0]-p2.m/LA.norm(p2.r)+p2.m/r_dmin[1])*ms*ms/rs

        e_tot = e_kin + e_pot# + e_correct
    #    print e_kin, e_pot, e_correct

        if e_tot < 0:
            c1 = p1.core
            c2 = p2.core

            #Do calculations for bound planets
            com  = (p1.m*p1.r + p2.m*p2.r)/(p1.m+p2.m)
            vcom = (p1.m*p1.v + p2.m*p2.v)/(p1.m+p2.m)
            #aei1 = mercury_routines.calculate_orbital_elements(p2.m, p1.m, p1.r-com, c1.v - vcom)#v_cm1)
            #aei2 = mercury_routines.calculate_orbital_elements(p1.m, p2.m, p2.r-com, c2.v - vcom)#v_cm2)
            aei = calculate_ae_xyz([p1.m, p2.m], p1.r - p2.r, p1.v - p2.v)
            #print aei['a']*1.496*10.**13./(6.*10.**8)*(1.-aei['e'])
            bound_orbits.append({'t': snapshot.unitless_time, 'aei': aei, 'e': e_tot, 'sep': r*rs, 'm': [p1.m, p2.m]})

        else:
            bound_orbits.append({'t': snapshot.unitless_time, 'aei': None, 'e': e_tot, 'sep': r*rs, 'm': [p1.m, p2.m]})

    return bound_orbits




# m, r, v are the masses, positions, and velocities of the two bodies
def calculate_ae_xyz(m, r, v):

    m1, m2 = m

    mu = sum(m)
    red_mass = m1*m2/mu

    energy = np.linalg.norm(v)**2./2.-mu/np.linalg.norm(r)
    a = -mu/2./energy

    evec = ((np.linalg.norm(v)**2.-mu/np.linalg.norm(r))*r-np.dot(r, v)*v)/mu
    e = np.linalg.norm(evec)

    # Convert semi-major axis to AU
    a = a*0.00465

    return {'a': a, 'e': e}




def get_planet_energy(planets, star):

    g   = 6.67390*math.pow(10,-8)                               # gravitational constant [cgs]
    rs  = 6.9599*math.pow(10,10)                                # unit of length [cm]
    ms  = 1.9891*math.pow(10,33)                                # unit of mass [g]
    ts  = math.sqrt(math.pow(rs,3)/(g*ms))                      # unit of time [s]

    u = 0.
    k = np.sum([planet.m*LA.norm(planet.v)**2./2.*ms*rs*rs/ts/ts for planet in planets])
    w = -np.sum([g*star.m*planet.m/LA.norm(planet.r-star.r)*ms*ms/rs for planet in planets])

    p1, p2 = planets
    w = w - g*p1.m*p2.m/LA.norm(p1.r-p2.r)*ms*ms/rs


    return {'k': k, 'u': u, 'w': w}





# returns the initial and final times for each collision
def get_collision_times(t, d, dcrit, dhit):

    d_t = {_t: _d for _t, _d in zip(t, d)}

    # calculate time to first contact
    t_0 = [_t for _t, _d1, _d2 in zip(t[:-1], d[:-1], d[1:]) if _d2 < dcrit < _d1]
    if not t_0: t_0 = [min(t)]

    # calculate time at minimum distance
    def get_tmin(x, y):
        if len(x) == 0: return []
        return [_x for _x, _y1, _y2 in zip(x[:-1], y[:-1], y[1:]) if _y1 < _y2]

    t_min = [min(get_tmin(t[t>_t_0], d[t>_t_0])) if len(get_tmin(t[t>_t_0], d[t>_t_0])) > 0 else max(t) for _t_0 in t_0]

    # calculate time that ends first contact
    def get_t_1(x, y):
        if len(x) == 0: return []
        return [_x for _x, _y1, _y2 in zip(x[:-1], y[:-1], y[1:]) if _y1 > dcrit or _y2 < _y1]

    t_1 = [min(get_t_1(t[t>=_t_min], d[t>=_t_min])) if len(get_t_1(t[t>=_t_min], d[t>=_t_min])) > 0 else max(t) for _t_min in t_min]

#    if len(t_0) == 0: return [[None, None]]

    if len([1 for _t_min in t_min if d_t[_t_min] < dhit]) == 0: return [[None, None]]

    return [[_t_0, _t_1] for _t_0, _t_min, _t_1 in zip(t_0, t_min, t_1,) if d_t[_t_min] < dhit]




# returns minimum distance interpolating the incoming and outgoing trajectories
# input dictionary where time: separation of planets
def find_d_min(d):

    #def extrap1d(interpolator):
    #    xs = interpolator.x
    #    ys = interpolator.y

    #    def pointwise(x):
    #        if x < xs[0]:
    #            return ys[0] + (x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
    #        elif:
    #            return ys[-1] + (x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
    #        else:
    #            return interpolator(x)

    #    def ufunclike(xs):
    #        return array(map(pointwise, array(xs)))

    #    return ufunclike

    # get the times of the two minimum separations
    t_d = {d[t]: t for t in d}
    t_min = t_d[min(t_d)]
    t_d.pop(min(t_d), None)
    t_min2 = t_d[min(t_d)]
    t_lower = min([t_min, t_min2])

    d_in = {t: d[t] for t in d if t <= t_lower}
    d_out = {t: d[t] for t in d if t > t_lower}

    t1, t2 = sorted(d_in)[-2:]
    m1 = (d_in[t2] - d_in[t1])/(t2 - t1)
    b1 = d_in[t2] - m1*t2

    if len(d_out) > 1:
        print '____________Fine'
        t3, t4 = sorted(d_out)[:2]
        m2 = (d_out[t4] - d_out[t3])/(t4 - t3)
        b2 = d_out[t4] - m2*t4
    else:
        print '____________Not enough data'
        m2 = -1./m1
        b2 = -b1

    x_sol = (b2 - b1)/(m1 - m2)
    y_sol = m1*x_sol + b1

    f = interp1d(sorted(d), [d[t] for t in sorted(d)], kind = 'cubic')

    def f2(x):
        if x < min(d): return 10.**30.
        if x > max(d): return 10.**30.
        return f(x)

    bnds = ((t_min, t_min2),)
    if t_min > t_min2: bnds = ((t_min2, t_min),)

    x_sol2 = op.minimize(f2, (t_min+t_min2)/2., bounds = bnds)['x'][0]
    y_sol2 = f2(x_sol2)

    print y_sol, y_sol2, d[t_min], d[t_min2]

    return {'linear': y_sol, 'cubic': y_sol2}




def find_a_range(m1, m2, a1, a2, e1, e2, M):

    h_tot = M*m1/(M+m1)*((M+m1)*a1*(1.-e1**2.))**.5+M*m2/(M+m2)*((M+m2)*a2*(1.-e2**2.))**.5
    e_tot = -0.5*(M*m1/a1+M*m2/a2)
    max_e = .3

    #___As a function of a2, find a1

    # calculate angular momentum
    def h(a, e, m):
        return M*m/(M+m)*((M+m)*a*(1.-e**2.))**.5

    def find_min_a1(x):

        if x < 0: return 1000.
        a1_new = -M*m1/(e_tot*2.+M*m2/x)
        # ensure solution exists for min and max eccentricities
        h0 = h(a1_new, 0., m1) + h(x, 0., m2)
        h2 = h(a1_new, max_e, m1) + h(x, max_e, m2)
        if h0 < h_tot: return 1000.
        if h2 > h_tot: return 1000.
        if a1_new < 0: return 1000.
        return a1_new


    def find_max_a1(x):
        if x < 0: return 1000.
        a1_new = -M*m1/(e_tot*2.+M*m2/x)
        # ensure solution exists for min and max eccentricities
        h0 = h(a1_new, 0., m1) + h(x, 0., m2)
        h2 = h(a1_new, max_e, m1) + h(x, max_e, m2)
        if h0 < h_tot: return 1000.
        if h2 > h_tot: return 1000.
        if a1_new < 0: return 1000.
        return -a1_new


    # as a function of a2, find a1
    def find_min_a2(x):
        if x < 0: return 1000.
        a2_new = -M*m2/(e_tot*2.+M*m1/x)
        # ensure solution exists for min and max eccentricities
        h0 = h(a2_new, 0., m2) + h(x, 0., m1)
        h2 = h(a2_new, max_e, m2) + h(x, max_e, m1)
        if h0 < h_tot: return 1000.
        if h2 > h_tot: return 1000.
        if a2_new < 0: return 1000.
        return a2_new

    def find_max_a2(x):
        if x < 0: return 1000.
        a2_new = -M*m2/(e_tot*2.+M*m1/x)
        # ensure solution exists for min and max eccentricities
        h0 = h(a2_new, 0., m2) + h(x, 0., m1)
        h2 = h(a2_new, max_e, m2) + h(x, max_e, m1)
        if h0 < h_tot: return 1000.
        if h2 > h_tot: return 1000.
        if a2_new < 0: return 1000.
        return -a2_new

    min_a1 = find_min_a1(op.minimize(find_min_a1, a2).x[0])
    min_a2 = find_min_a2(op.minimize(find_min_a2, a1).x[0])
    max_a1 = -find_max_a1(min_a2)
    max_a2 = -find_max_a2(min_a1)

    return [[min_a1, max_a1], [min_a2, max_a2]]




def find_p_ratios(runs, M):

    m = [run['m'] for run in runs]
    a = [run['a'] for run in runs]
    e = [run['e'] for run in runs]

    a_range = [find_a_range(_m[0], _m[1], _a[0], _a[1], _e[0], _e[1], M) for _m, _a, _e in zip(m, a, e)]

    p_ratio = [[(_a[1][0]/_a[0][1])**1.5, (_a[1][1]/_a[0][0])**1.5] for _a in a_range]

    return p_ratio




def find_p_ratios2(m, a, e, M):

    a_range = find_a_range(m[0], m[1], a[0], a[1], e[0], e[1], M)

    p_ratio = [(a_range[1][0]/a_range[0][1])**1.5, (a_range[1][1]/a_range[0][0])**1.5]

    return p_ratio





def get_mass_lost(run):

    t0, t1 = run.t_collisions[0]

    m_t = {orbit['t']: [m0 - _m for m0, _m in zip(run.initial['m'], orbit['m'])] for orbit in run.orbit}

    return m_t[t1]




def calculate_AMD(mstar, m12, aei1, aei2, i = [0., 0.]):

    g = 6.67390*math.pow(10,-8)                               # gravitational constant [cgs]
    ms = 1.9891*math.pow(10,33)                                # unit of mass [g]
    au = 1.496*10.**13.

    a1, e1 = aei1
    a2, e2 = aei2
    m1, m2 = m12

    C = [m*a**.5*(1.-(1.-e**2)**.5*math.cos(_i*math.pi/180.)) for m, a, e, _i in zip([m1, m2], [a1, a2], [e1, e2], i)]

    return sum(C)*ms**1.5*mstar**.5*g**.5*au**.5




def calculate_AMD_run(mstar, run):

    g = 6.67390*math.pow(10,-8)                               # gravitational constant [cgs]
    ms = 1.9891*math.pow(10,33)                                # unit of mass [g]
    au = 1.496*10.**13.

    C = [m*a**.5*(1.-(1.-e**2)**.5*math.cos(i*math.pi/180.)) for m, a, e, i in zip(run['m'], run['a'], run['e'], run['i'])]

    return sum(C)*ms**1.5*mstar**.5*g**.5*au**.5




def calculate_AMD_min_run(mstar, run):

    g = 6.67390*10.**-8                               # gravitational constant [cgs]
    ms = 1.9891*math.pow(10,33)                                # unit of mass [g]
    au = 1.496*10.**13.
    me = 5.9722*10.**27.

    C = ms**1.5*au**.5
    C0 = ms**.5*me*mstar**.5*g**.5*au**.5*(sum(run['m']))

    h = sum([mstar*run['m'][i]/(mstar + run['m'][i])*(g*(mstar + run['m'][i])*run['a'][i]*(1.-run['e'][i]**2.))**.5*math.cos(run['i'][i]*math.pi/180.) for i in range(2)])*C

    # angular momentum at 0 ecc, inner a, and pratio
    def h2(r, p):
        energy = -g/2.*sum([(mstar+_m)*_m/_a for _m, _a in zip(r['m'], r['a'])])*ms
        m1, m2 = r['m']
        a1 = -g/2./energy*((mstar+m1)*m1+(mstar+m2)*m2/p**(2./3.))*ms
        a2 = p**(2./3.)*a1
        return sum([mstar*r['m'][i]/(mstar + r['m'][i])*(g*(mstar + r['m'][i])*a)**.5 for i, a in zip(range(2), [a1, a2])])*C

    p_range = find_p_ratios([run], mstar)[0]
    # create list of possible pratios

    def AMD(x):
        if x < p_range[0]: return 1000.
        if x > p_range[1]: return 1000.
        return (h2(run, x) - h)/C0

    res = op.minimize(AMD, 1.)

    return res['x']




