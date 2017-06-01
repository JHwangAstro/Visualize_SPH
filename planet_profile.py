

import math
import numpy as np

class PlanetProfile:

    def __init__(self, particles, x):

        nbins = 100

        r = np.array([np.linalg.norm(particle.r - x) for particle in particles])
        m = np.array([particle.m for particle in particles])
        u = np.array([particle.u for particle in particles])
        T = np.array([particle.T for particle in particles])
        p = np.array([particle.p for particle in particles])
        j = np.array([particle.j for particle in particles])
        ind = np.array([i for i in range(len(particles))])

        # for negative p-values use p of the closest neighbor
        rp = {_r: _p for _r, _p in zip(r[p>0.], p[p>0.])}
        for i, _r, _p in zip(ind[p<0.], r[p<0.], p[p<0.]):
            p[i] = rp[rp.keys()[np.abs(rp.keys() - _r).argmin()]]


        np.savetxt('/jhwang/SPH_Visual/debug_p.txt', np.c_[r,m,u,T,p], fmt = '%1.4e')

        my_dict = {'r': r, 'm': m, 'u': u, 'T': T, 'p': p, 'j': j}

        rbins = np.linspace(0., np.max(r), nbins+1)

        print 'checking pressure: ', min(p), max(p)
        pbins = np.logspace(math.log10(min(p[p>0.])), math.log10(max(p[p>0.])), nbins+1)

        self.r_profile = {k: [sum(v[(r0<r) & (r<r1)]*m[(r0<r) & (r<r1)])/sum(m[(r0<r) & (r<r1)]) \
                for r0, r1 in zip(rbins[0:-1], rbins[1:]) if len(m[(r0<r) & (r<r1)]) > 0] for k, v in my_dict.iteritems()}
        self.p_profile = {k: [sum(v[(p0<p) & (p<p1)]*m[(p0<p) & (p<p1)])/sum(m[(p0<p) & (p<p1)]) \
                for p0, p1 in zip(pbins[0:-1], pbins[1:]) if len(m[(p0<p) & (p<p1)]) > 0] for k, v in my_dict.iteritems()}




#    def get_temperature(self, p):

#        kb = 1.38064852*10**-16.    #erg/K

        # for core particles
#        if not p.cc % 2: return 0.

#        if return (p.u - p.K*p.rho**(p.gam-1.)/(p.gam-1.))*(p.mmw/(p.beta*kb))











