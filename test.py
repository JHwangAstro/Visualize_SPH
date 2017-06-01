from datetime import datetime
import numpy as np

# Smoothing kernel
def w(r,h):
    u = r/h
    #if (u>=2.): return 1./r
    if (u<=1.): return -2./h*(u**2./3.-3.*u**4./20.+u**5./20.)+7./5./h
    if (u<=2.): return -1./15./r-(4.*u**2./3-u**3.+3*u**4./10.-u**4./30.-u**5./30.)/h+8./5./h
    return 0.

x = np.random.random(1000)*3.

starttime = datetime.now()

y = [w(i,1.) for i in x]

print 'first try: ', datetime.now()-starttime

startime = datetime.now()

h = 1.
u = x
np.piecewise(u, [u<=1., u<=2., u>2.], [lambda u: -2./h*(u**2./3.-3.*u**4./20.+u**5./20.)+7./5./h, lambda u: -1./15./u-(4.*u**2./3-u**3.+3*u**4./10.-u**4./30.-u**5./30.)/h+8./5./h, 0.])

print 'second try: ', datetime.now()-starttime
