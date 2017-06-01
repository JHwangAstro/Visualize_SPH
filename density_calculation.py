
import math
import numpy as np
from scipy import interpolate

def read_table(fname):

    f = open(fname, 'rb')

    me = []
    r = []

    for line in f:
        x = line.split()
        me.append(float(x[0]))
        r.append(float(x[1]))

    return me, r




def get_density(planet):

    ms_to_me = 1.99*10.**33./(5.9722*10.**27.)
    me_to_g = 5.9722*10.**27.
    ms_to_g = 1.99*10.**33.
    re_to_cm = 6.371*10.**8.

    base = '/jhwang/SPH_Visual/'
    fnames = ['Tables/2Mc_R_Me.txt', 'Tables/5Mc_R_Me.txt', 'Tables/10Mc_R_Me.txt']
    fnames = [base + f for f in fnames]
    core_masses = [2., 5., 10.]

    env_masses = {}
    rad = {}
    for fname, mc in zip(fnames, core_masses):
        me, r = read_table(fname)
        env_masses[mc] = me
        for _me, _r in zip(me, r):
            rad[(mc, _me)] = float(_r)

    p_mc = planet.core_mass_fraction*planet.m*ms_to_me
    p_me = planet.gas_mass_fraction*planet.m*ms_to_me

    # if planet core mass is outside of table range, return -1
    if not min(core_masses) < p_mc < max(core_masses): return -1.

    # find bracketing core masses:
    mc0 = max(np.array(core_masses)[core_masses < p_mc])
    mc1 = min(np.array(core_masses)[core_masses > p_mc])
    me00 = max(np.array(env_masses[mc0])[env_masses[mc0] < p_me])
    me01 = min(np.array(env_masses[mc0])[env_masses[mc0] > p_me])
    me10 = max(np.array(env_masses[mc1])[env_masses[mc1] < p_me])
    me11 = min(np.array(env_masses[mc1])[env_masses[mc1] > p_me])

    r00 = rad[(mc0, me00)]
    r01 = rad[(mc0, me01)]
    r10 = rad[(mc1, me10)]
    r11 = rad[(mc1, me11)]

    r0 = (me01-p_me)/(me01-me00)*r00 + (p_me-me00)/(me01-me00)*r01
    r1 = (me11-p_me)/(me11-me10)*r10 + (p_me-me10)/(me11-me10)*r11

    p_r = (mc1-p_mc)/(mc1-mc0)*r0 + (p_mc-mc0)/(mc1-mc0)*r1

#    if not mc0 or not mc1 or not me0 or not me1: return None

    return planet.m*ms_to_g/(4.*math.pi*(p_r*re_to_cm)**3./3.)






