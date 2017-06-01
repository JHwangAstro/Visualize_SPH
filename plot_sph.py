import numpy as np
import math
from scipy.spatial import cKDTree as KDTree
#from scipy.spatial.distance import cdist
from sklearn.metrics.pairwise import euclidean_distances
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess

from datetime import datetime

import os



def call_splash(sdir, idir, ss, new_min, new_max, prefix):

    limits = sdir+"splash.limits"

    # number of limits
    m = 12

    os.chdir(sdir)
    xmin = np.zeros(m)
    xmax = np.zeros(m)
    if os.path.exists(limits):
        data = np.loadtxt(limits)
        for i in range(len(data)):
            xmin[i] = data[i,0]
            xmax[i] = data[i,1]
    for i in range(m):
        if new_min[i] == 0.: new_min[i] = xmin[i]
        if new_max[i] == 0.: new_max[i] = xmax[i]
    np.savetxt(limits, np.c_[new_min, new_max])

    fname = 'out' + format(ss, "04") + '.sph'
    s = "printf '0\nq' | ~/splash/bin/jsplash -x 1 -y 2 -r 8 -vec 0 -c 0 -dev '"+ prefix + "_" + format(ss, "04") + ".png/png' " + fname

    subprocess.call(s, shell=True)



# Assume particles have an x, y, z, and h attribute
def plot_sph(particles, fname, size=[1000.,1000.], origin=[0.,0.], res=[500,500]):

    # Smoothing kernel
    def w(r,h):
        u = r/h
#        if (u<=1.): return -2./h*(u**2./3.-3.*u**4./20.+u**5./20.)+7./5./h
#        if (u<=2.): return -1./15./r-(4.*u**2./3-u**3.+3*u**4./10.-u**4./30.-u**5./30.)/h+8./5./h
#        return 0.

        u2 = u*u
        u3 = u2*u
        u4 = u3*u
        u5 = u4*u

        f1 =  -2./h*(u2/3.-3.*u4/20.+u5/20.)+7./5./h
        f2 =  -1./15.*r**(-1.)-(4.*u2/3.-u3+4*u4/15.-u5/30.)/h+8./5./h
#        f1 = (1.-1.5*u**2.+0.75*u**3.)/math.pi/h**3.
#        f2 = (0.25*(2.-u)**3.)/math.pi/h**3.
        answer = np.zeros(r.shape)
        answer[np.where(u<=2)] = f2[np.where(u<=2)]
        answer[np.where(u<=1)] = f1[np.where(u<=1)]
        return answer

#        return np.piecewise(u,r,h, [u<=1.
#        return  (u<=1)*(-2./h*(u**2./3.-3.*u**4./20.+u**5./20.)+7./5./h) + \
#                ((u>1)&(u<=2))*(-1./15./r-(4.*u**2./3-u**3.+3*u**4./10.-u**4./30.-u**5./30.)/h+8./5./h)

    nx = res[0]
    ny = res[1]

    x = np.array([particles[k].x for k in particles])
    y = np.array([particles[k].y for k in particles])
    h = np.array([particles[k].h for k in particles])
    m = np.array([particles[k].m for k in particles])

    width  = size[0]/nx
    height = size[1]/ny
    min_x  = origin[0]-size[0]/2.
    max_x  = origin[0]+size[0]/2.
    min_y  = origin[1]-size[1]/2.
    max_y  = origin[1]+size[1]/2.

    lower_x = (x-2.*h-min_x)//width
    upper_x = nx - (max_x-(x+2.*h))//width
    lower_y = (y-2.*h-min_y)//height
    upper_y = ny - (max_y-(y+2.*h))//height

    pixel_x = np.array([origin[0] - size[0]/2 + (i+0.5)*width  for i in xrange(nx)])
    pixel_y = np.array([origin[1] - size[1]/2 + (j+0.5)*height for j in xrange(ny)])
    pixel_xy = np.array([[pixel_x[i], pixel_y[j]] for j in xrange(ny) for i in xrange(nx)])

    ids = np.arange(len(x))[(lower_x < nx-1) & (upper_x > 0) & (lower_y < ny-1) & (upper_y > 0)]

    # Define boundaries after finding particle ids
    lower_x = np.maximum(0, lower_x.astype(int))
    upper_x = np.minimum(nx-1, upper_x.astype(int))
    lower_y = np.maximum(0, lower_y.astype(int))
    upper_y = np.minimum(ny-1, upper_y.astype(int))

    tstart = datetime.now()

    print 'checkpoint 2'

    data = np.full((nx, ny), 10.**-21., dtype = float)

    #r = [np.reshape([((pixel_x[i]-x[k])**2.+(pixel_y[j]-y[k])**2.)**0.5 \
    #    for i in xrange(lower_x[k], upper_x[k]) for j in xrange(lower_y[k], upper_y[k])], \
    #    (upper_x[k]-lower_x[k], upper_y[k]-lower_y[k])) for k in ids]
    for k in ids:
        pixel_xy = np.array([[pixel_x[i], pixel_y[j]] for i in xrange(lower_x[k], upper_x[k]) for j in xrange(lower_y[k],upper_y[k])])
        #r = np.reshape(cdist(np.array([[particles[k].x, particles[k].y]]), pixel_xy), [upper_y[k]-lower_y[k], upper_x[k]-lower_x[k]])
        #r = np.reshape(euclidean_distances(np.array([[particles[k].x, particles[k].y]]), pixel_xy), [upper_y[k]-lower_y[k], upper_x[k]-lower_x[k]])
        X = np.array([particles[k].x, particles[k].y])
        Y = pixel_xy
        r = np.sqrt(np.dot(X, X))-(np.dot(X, Y)+np.dot(X,Y))+np.dot(Y, Y)
        #data += w(r,h[k])*m[k]
        #r = np.reshape([((pixel_x[i]-x[k])**2.+(pixel_y[j]-y[k])**2.)**0.5 \
        #    for i in xrange(lower_x[k], upper_x[k]) for j in xrange(lower_y[k], upper_y[k])], \
        #    (upper_x[k]-lower_x[k], upper_y[k]-lower_y[k]))
        #data[lower_x[k]:upper_x[k], lower_y[k]:upper_y[k]] += w(r,h[k])*m[k]
        data[lower_y[k]:upper_y[k], lower_x[k]:upper_x[k]] += w(r,h[k])*m[k]
    print 'checkpoint 3'

    #data = np.full((nx, ny), 10.**-21., dtype = float)

    #for k in ids:
    #    data[lower_x[k]:upper_x[k], lower_y[k]:upper_y[k]] += z[k]

    data2 = data.flatten()
    data2 = data2[data2>10.**-20.]
    max_z = 10.**-2.5#np.amax(data2)
    min_z = 10.**-9.#np.amin(data2)

    plt.imshow(data, cmap = 'afmhot', norm=mpl.colors.LogNorm(vmin=min_z, vmax=max_z), interpolation='bicubic', origin = 'lower', \
            extent = [min_x, max_x, min_y, max_y])
    plt.colorbar()
    plt.show()
    plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
    plt.clf()






