
import numpy as np
import os
import subprocess

class NBodyOrbits:

    def __init__(self, ndir, fname_aei, fname_ce, MS, radii = None):

        au_to_re = 1.496*10.**13./(6.371*10.**8.)

        self.exists = True
        self.nbody_Period_Ratio = -1.

        # if there are mercury files
        if not os.path.isfile(ndir + fname_aei[0]) and not os.path.isfile(ndir + 'xv.out'):
            self.exists = False
            return None

        if not os.path.isfile(ndir + fname_ce[0]) and not os.path.isfile(ndir + 'ce.out'):
            self.exists = False
            return None

        # check to see if xv.out is newer if it exists
        if os.path.isfile(ndir + fname_aei[0]):
            if os.path.isfile(ndir + 'xv.out'):
                if os.path.getmtime(ndir + fname_aei[0]) + 3600. < os.path.getmtime(ndir + 'xv.out'):
                    aei_files = [x for x in os.listdir(ndir) if '.aei' in x]
                    for aei in aei_files:
                        os.rename(ndir + aei, ndir + aei + '_backup')
                    os.chdir(ndir)
                    subprocess.call([ndir + '/element6'])
        else:
            os.chdir(ndir)
            subprocess.call([ndir + '/element6'])

        if os.path.isfile(ndir + fname_ce[0]):
            if os.path.isfile(ndir + 'ce.out'):
                # currently not rerunning
                if os.path.getmtime(ndir + fname_ce[0]) + 3600000. < os.path.getmtime(ndir + 'ce.out'):
                    clo_files = [x for x in os.listdir(ndir) if '.clo' in x]
                    for clo in clo_files:
                        os.rename(ndir + clo, ndir + clo + '_backup')
                    os.chdir(ndir)
                    subprocess.call([ndir + '/close6'])
        else:
            os.chdir(ndir)
            subprocess.call([ndir + '/close6'])

        self.t, self.a, self.e, self.i, self.m = self.get_orbital_params([ndir + f for f in fname_aei])

        t_collision = max(self.t[0])
        if os.path.isfile(ndir + fname_ce[0]):
            t_collision = self.get_collision_time([ndir + f for f in fname_ce], radii)

        if radii:
            self.a = [np.ma.masked_where(t > t_collision, a) for t, a in zip(self.t, self.a)]
            self.e = [np.ma.masked_where(t > t_collision, e) for t, e in zip(self.t, self.e)]
            self.i = [np.ma.masked_where(t > t_collision, i) for t, i in zip(self.t, self.i)]
            self.m = [np.ma.masked_where(t > t_collision, m) for t, m in zip(self.t, self.m)]
            self.t = [np.ma.masked_where(t > t_collision, t) for t in self.t]

        self.e = [np.ma.masked_where(a < 0., e) for e, a in zip(self.e, self.a)]
        self.i = [np.ma.masked_where(a < 0., i) for i, a in zip(self.i, self.a)]
        self.m = [np.ma.masked_where(a < 0., m) for m, a in zip(self.m, self.a)]
        self.t = [np.ma.masked_where(a < 0., t) for t, a in zip(self.t, self.a)]
        self.a = [np.ma.masked_where(a < 0., a) for a in self.a]

        self.p_ratio = [(a1/a0)**1.5*((MS+m1)/(MS+m0))**-.5 for a0, a1, m0, m1 in zip(self.a[0], self.a[1], self.m[0], self.m[1])]
        self.p_ratio = np.array([x if x > 1. else 1./x for x in self.p_ratio])

#        if self.nbody_a0 < 0. or self.nbody_a1 < 0.:
#            self.nbody_Period_Ratio = -1.




    def get_orbital_params(self, fname):

        planets = [np.loadtxt(f, skiprows = 4, usecols = (0,1,2,3,7)) for f in fname]
        planets = [map(list, zip(*p)) for p in planets]
        orbits = map(list, zip(*[[np.array(elements) for elements in p] for p in planets]))
        return orbits




    def get_collision_time(self, fname, r):

        au_to_cm = 1.496*10.**13.

        #Read in data files
        close_encounters = [np.loadtxt(f, skiprows = 4, usecols = (0,2)) for f in fname]
        t = [np.array([float(x[0]) for x in ce]) for ce in close_encounters]
        d = [np.array([float(x[1]) for x in ce]) for ce in close_encounters]

        if len(t[0]) > 0:

            for _t, _d in zip(t[0], d[0]):

                if isinstance(r, list):
                    if _d*au_to_cm < r[0]+r[1]: return _t
                else:
                    return _t

            return None








