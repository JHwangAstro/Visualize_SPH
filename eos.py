
import math
from scipy.optimize import fsolve
from scipy.special import hyp2f1

def get_temperature(rho, u, mmw, K, gam, b, x, cc, ul0, uh0, u0):

    boltzcgs = 1.380658*10**-16.
    gravconst = 6.67390*10**-8.
    rs  = 6.9599*math.pow(10,10)                                # unit of length [cm]
    ms  = 1.9891*math.pow(10,33)                                # unit of mass [g]
    rho_m = 4.26
    gam_m = 1./0.549

    rhocgs = rho*ms/rs**3.
    ucgs = u*gravconst*ms/rs

    # core particles
    if not cc % 2: return 0.
    if cc % 2 and x == 1: return 0.

    # pure gas particle
    if not x:
        return (ucgs-K*rhocgs**(gam-1.)/(gam-1.))/b/boltzcgs*mmw
    # mixed gas-mantle particle
    else:

        ul = ucgs*ul0/u0
        uh = ucgs*uh0/u0

        # get component densities using pressure continuity and conservation of mass
        def p1(_u, _rho):
            if _rho < rho_m: return -1000000.
            return _u*(gam_m-1.)*_rho*(1.-rho_m/_rho)**gam_m/ \
                hyp2f1(1.-gam_m, -gam_m, 2.-gam_m, rho_m/_rho)

        def p2(_u, _rho):
            return _u*_rho/b + K*_rho**gam*(1.-_rho/(gam-1.)/b)

        dl, dh = get_component_densities(rhocgs, x, ul, uh, rho_m, p1, p2)
#        print dl, dh
#        print ucgs*ul0/u0, K*dl**(gam-1.)/(gam-1.)
        return max([0., (ucgs*ul0/u0-K*dl**(gam-1.)/(gam-1.))/b/boltzcgs*mmw])





def get_pressure(rho, u, mmw, K, gam, b, x, cc, ul0, uh0, u0):

    gravconst = 6.67390*10**-8.
    rs  = 6.9599*math.pow(10,10)                                # unit of length [cm]
    ms  = 1.9891*math.pow(10,33)                                # unit of mass [g]
    rhocgs = rho*ms/rs**3.
    ucgs = u*gravconst*ms/rs

    rho_c = 8.3
    gam_c = 1./0.528
    rho_m = 4.26
    gam_m = 1./0.549

    # core particle
    if not cc % 2:
        # pure iron-core
        if x == 1.: return ucgs*(gam_c-1.)*rhocgs*(1.-rho_c/rhocgs)**gam_c/ \
                hyp2f1(1.-gam_c, -gam_c, 2.-gam_c, rho_c/rhocgs)

        # pure silicate-mantle
        if x == 0.: return ucgs*(gam_m-1.)*rhocgs*(1.-rho_m/rhocgs)**gam_m/ \
                hyp2f1(1.-gam_m, -gam_m, 2.-gam_m, rho_m/rhocgs)

        # mixed core-mantle particle
        ul = ucgs*ul0/u0
        uh = ucgs*uh0/u0

        # get component densities using pressure continuity and conservation of mass
        def p1(_u, _rho):
            if _rho < rho_c: return -5000000.
            return _u*(gam_c-1.)*_rho*(1.-rho_c/_rho)**gam_c/ \
                hyp2f1(1.-gam_c, -gam_c, 2.-gam_c, rho_c/_rho)

        def p2(_u, _rho):
            if _rho < rho_m: return -1000000.
            return _u*(gam_m-1.)*_rho*(1.-rho_m/_rho)**gam_m/ \
                hyp2f1(1.-gam_m, -gam_m, 2.-gam_m, rho_m/_rho)

        dl, dh = get_component_densities(rhocgs, x, ul, uh, rho_c, p1, p2)

        # solve for density, choice of component should not matter
        return p1(uh, dh)

    # gas particle
    else:
        # pure silicate-mantle
        if x == 1.: return ucgs*(gam_m-1.)*rhocgs*(1.-rho_m/rhocgs)**gam_m/ \
                hyp2f1(1.-gam_m, -gam_m, 2.-gam_m, rho_m/rhocgs)

        # pure gas-particle
        if x == 0.: return rhocgs*ucgs/b + K*rhocgs**(gam-1.)/(gam-1.)*(1.-rhocgs/b)

        ul = ucgs*ul0/u0
        uh = ucgs*uh0/u0

        # get component densities using pressure continuity and conservation of mass
        def p1(_u, _rho):
            if _rho < rho_m: return -1000000.
            return _u*(gam_m-1.)*_rho*(1.-rho_m/_rho)**gam_m/ \
                hyp2f1(1.-gam_m, -gam_m, 2.-gam_m, rho_m/_rho)

        def p2(_u, _rho):
            return _u*_rho/b + K*_rho**gam*(1.-_rho/(gam-1.)/b)

        dl, dh = get_component_densities(rhocgs, x, ul, uh, rho_m, p1, p2)

        # solve for density, choice of component should not matter
        return p1(uh, dh)





def get_component_densities(rho, x, ul, uh, rho_crit, p1, p2):

    def calculate_heavy_density(rho_h):

        # conservation of mass
        rho_l = (1.-x)*rho*rho_h/(rho_h-x*rho)

        return p1(uh, rho_h) - p2(ul, rho_l)

    rho_h_solution = fsolve(calculate_heavy_density, max([rho, rho_crit*1.1]), xtol = 1e-10, factor = 1.)
    rho_l_solution = (1.-x)*rho*rho_h_solution/(rho_h_solution-x*rho)

    return rho_l_solution, rho_h_solution








def get_cooling_time(m, u, rho, T):

    gravconst = 6.67390*10**-8.
    rs  = 6.9599*math.pow(10,10)                                # unit of length [cm]
    ms  = 1.9891*math.pow(10,33)                                # unit of mass [g]
    rhocgs = rho*ms/rs**3.
    ucgs = u*gravconst*ms/rs
    sigma = 5.6704*10.**-5.

    A = 4.*math.pi*(3.*m*ms/4./math.pi/rhocgs)**(2./3.)
    uraddot = A*sigma*T**4

    return ucgs*m*ms/uraddot/60/60/24/365.2422    #in years
