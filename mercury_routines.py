import math
import numpy as np

# in:  r,v
# out: a,e,i,g,n,M
# adapted from mco_x2el subroutine in Mercury (Chambers 1999)
def calculate_orbital_elements(m_primary, m_secondary, r0, v0):

    G = 1.

    if m_secondary == 0: return -1

    gm = G*(m_primary+m_secondary)
    x = r0[0]
    y = r0[1]
    z = r0[2]
    u = v0[0]
    v = v0[1]
    w = v0[2]

    hx = y*w - z*v
    hy = z*u - x*w
    hz = x*v - y*u
    h2 = hx*hx + hy*hy + hz*hz
    v2 = u*u + v*v + w*w
    rv = x*u + y*v + z*w
    r = (x*x + y*y + z*z)**0.5
    h = h2**0.5
    s = h2/gm

    # inclination and node
    ci = hz/h
    if (abs(ci) < 1.):
        i = math.acos(ci)
        n = math.atan2(hx, -hy)
        if (n < 0.): n = n + 2.*math.pi
    else:
        if (ci > 0.): i = 0.
        if (ci < 0.): i = math.pi
        n = 0.

    # cccentricity and perihelion distance
    temp = 1.0 + s*(v2/gm - 2./r)
    if (temp < 0.):
        e = 0.
    else:
        e = temp**.5
    q = s/(1.+e)

    # true longitude
    if (hy != 0.):
        to = -hx/hy
        temp = (1.-ci)*to
        tmp2 = to*to
        true = math.atan2((y*(1.+tmp2*ci)-x*temp), (x*(tmp2+ci)-y*temp))
    else:
        true = atan2(y * ci, x)

    if (ci < 0): true = true + math.pi

    if (e < 0.):
        p = 0.
        l = 1
    else:
        ce = (v2*r-gm)/(e*gm)

        # mean anomaly for ellipse
        if (e < 1.):
            if (abs(ce) > 1.): ce = 1.*np.sign(ce)
            bige = math.acos(ce)
            if (rv < 0.): bige = 2.*math.pi - bige
            l = bige - e*math.sin(bige)
        # mean anomaly for hyperbola
        else:
            if (ce < 1.): ce = 1.
            bige = math.log(ce+(ce*ce-1.)**2.)
            if (rv < 0.): bige = -bige
            l = e*math.sinh(bige) - bige

        # Longitude of perihelion
        cf = (s-r)/(e*r)
        if (abs(cf) > 1): cf = 1.*np.sign(cf)
        f = math.acos(cf)
        if (rv < 0.): f = 2.*math.pi - f
        p = true - f
        p = (p+4.*math.pi) % (2.*math.pi)

    if (l < 0.): l = l + 2.*math.pi
    if (l > 2.*math.pi): l = l%(2.*math.pi)

    # Convert semi-major axis to AU
    a = q/(1.-e)*0.00465

    return {'a': a, 'e': e, 'i': i, 'g': p, 'n': n, 'M': l}



# adapted from mco_sine subroutine in Mercury (Chambers 1999)
def mco_sine(x):

    if (x > 0):
        x = x % (2.*math.pi)
    else:
        x = x % (2.*math.pi) + 2.*math.pi

    cx = math.cos(x)
    if (x > math.pi):
        sx = -(1.-cx*cx)**0.5
    else:
        sx = (1.-cx*cx)**0.5

    return cx, sx




def mco_kep(e, oldl):

    pi = 3.141592653589793
    twopi = 2.*pi
    piby2 = .5*pi

    # Reduce mean anomaly to lie in the range 0 < l < pi
    if (oldl >= 0):
        l = oldl%twopi
    else:
        l = (oldl%twopi)+twopi

    sign = 1.
    if (l > pi):
        l = twopi - l
        sign = -1.

    ome = 1. - e

    if (l >= .45 or e < .55):
    #
    # Regions A,B or C in Nijenhuis
    # -----------------------------
    #
    # Rough starting value for eccentric anomaly
        if (l < ome):
            u1 = ome
        else:
            if (l > (pi-1.-e)):
                u1 = (l+e*pi)/(1.+e)
            else:
                u1 = l + e

    # Improved value using Halley's method
        flag = (u1 > piby2)
        if flag:
            x = pi - u1
        else:
            x = u1

        x2 = x*x
        sn = x*(1. + x2*(-.16605 + x2*.00761))
        dsn = 1. + x2*(-.49815 + x2*.03805)
        if flag: dsn = -dsn
        f2 = e*sn
        f0 = u1 - f2 - l
        f1 = 1. - e*dsn
        u2 = u1 - f0/(f1 - .5*f0*f2/f1)
    else:
        #
        # Region D in Nijenhuis
        # ---------------------
        #
        # Rough starting value for eccentric anomaly
        z1 = 4.*e + .5
        p = ome / z1
        q = .5 * l / z1
        p2 = p*p
        z2 = math.exp(math.log((p2*p+q*q)**0.5+q)/1.5)
        u1 = 2.*q/(z2+p+p2/z2)
        #
        # Improved value using Newton's method
        z2 = u1*u1
        z3 = z2*z2
        u2 = u1-.075*u1*z3/(ome+z1*z2+.375*z3)
        u2 = l+e*u2*(3.-4.*u2*u2)
    #
    # Accurate value using 3rd-order version of Newton's method
    # N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!
    #
    # First get accurate values for u2 - sin(u2) and 1 - cos(u2)
    bigg = (u2 > piby2)
    if bigg:
        z3 = pi - u2
    else:
        z3 = u2

    big = (z3 > (.5*piby2))
    if big:
        x = piby2 - z3
    else:
        x = z3

    x2 = x*x
    ss = 1.
    cc = 1.

    ss = x*x2/6.*(1.-x2/20.*(1.-x2/42.*(1.-x2/72.*(1.-x2/110.*(1.-x2/156.*(1.-x2/210.*(1.-x2/272.)))))))
    cc = x2/2.*(1.-x2/12.*(1.-x2/30.*(1.-x2/56.*(1.-x2/ 90.*(1.-x2/132.*(1.-x2/182.*(1.-x2/240.*(1.-x2/306.))))))))

    if big:
        z1 = cc + z3 - 1.
        z2 = ss + z3 + 1. - piby2
    else:
        z1 = ss
        z2 = cc

    if bigg:
        z1 = 2.*u2 + z1 - pi
        z2 = 2. - z2

    f0 = l - u2*ome - e*z1
    f1 = ome + e*z2
    f2 = .5*e*(u2-z1)
    f3 = e/6.*(1.-z2)
    z1 = f0/f1
    z2 = f0/(f2*z1+f1)

    return sign*(u2+f0/((f3*z1+f2)*z2+f1))




# adapted from mco_el2x subroutine in Mercury (Chambers 1999)
def convert_el2x(aei, gm):

    a = aei['a']
    e = aei['e']
    i = aei['i']
    p = aei['g']
    n = aei['n']
    l = aei['M']

    # Change from longitude of perihelion to argument of perihelion
    g = p - n

    # Rotation factors
    ci, si = mco_sine(i)
    cg, sg = mco_sine(g)
    cn, sn = mco_sine(n)
    z1 = cg * cn
    z2 = cg * sn
    z3 = sg * cn
    z4 = sg * sn
    d11 = z1 - z4*ci
    d12 = z2 + z3*ci
    d13 = sg * si
    d21 = -z3 - z2*ci
    d22 = -z4 + z1*ci
    d23 = cg * si

    # Ellipse
    if (e < 1.):
        romes = (1.-e*e)**0.5
        temp = mco_kep(e, l)
        ce, se = mco_sine(temp)
        z1 = a*(ce - e)
        z2 = a*romes*se
        temp = (gm/a)**0.5/(1.-e*ce)
        z3 = -se*temp
        z4 = romes*ce*temp
    else:
        print 'not an ellipse!'
#        else:
# Parabola
#            if (e == 1.d0):
#                ce = orbel_zget(l)
#                z1 = q * (1.d0 - ce*ce)
#                z2 = 2.d0 * q * ce
#                z4 = sqrt(2.d0*gm/q) / (1.d0 + ce*ce)
#                z3 = -ce * z4
#            else:
# Hyperbola
#                romes = sqrt(e*e - 1.d0)
#                temp = orbel_fhybrid(e,l)
#                call mco_sinh (temp,se,ce)
#                z1 = a * (ce - e)
#                z2 = -a * romes * se
#                temp = sqrt(gm/abs(a)) / (e*ce - 1.d0)
#                z3 = -se * temp
#                z4 = romes * ce * temp

    x = d11 * z1  +  d21 * z2
    y = d12 * z1  +  d22 * z2
    z = d13 * z1  +  d23 * z2
    u = d11 * z3  +  d21 * z4
    v = d12 * z3  +  d22 * z4
    w = d13 * z3  +  d23 * z4

    return x, y, z, u, v, w





