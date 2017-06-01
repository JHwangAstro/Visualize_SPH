import numpy as np
import math
import random
import matplotlib.pyplot as plt
import mercury_routines

from planet_profile import PlanetProfile


class Planet:

    def __init__(self, particles, star):

        G = 1.0

        print len(particles)

        self.particles = [p.i for p in particles]

        self.m = np.sum([p.m for p in particles])
        self.m_gas = np.sum([p.m*(1.-p.xmix) for p in particles if p.cc % 2 and p.xmix < 1.])
        self.m_core = np.sum([p.m*p.xmix for p in particles if p.cc % 2 and p.xmix > 0.] + [p.m for p in particles if not p.cc % 2])

        self.utot = np.sum([p.m*p.u for p in particles])
        self.ugas = np.sum([p.m*(1.-p.xmix)*p.u for p in particles if p.cc % 2])
        self.mu = G*(star.m+self.m)

        self.r = self.calculate_center_of_mass(particles)
        self.v = self.calculate_velocity(particles)
        # temporary code to find largest distance of core particle:
        print '!!!!!!!max core particle: ',  max([np.linalg.norm(p.r-self.r) for p in particles if p.xmix > 0.])


        self.j_gas, self.specific_j_gas, self.e_rot_gas, self.specific_e_rot_gas, self.j_unit_vector_gas = \
                self.calculate_angular_momentum([(p.m*(1.-p.xmix), p.v - self.v, p.r) for p in particles if p.cc % 2 and p.xmix < 1.])
        self.j_core, self.specific_j_core, self.e_rot_core, self.specific_e_rot_core, self.j_unit_vector_core = \
                self.calculate_angular_momentum([(p.m*p.xmix, p.v - self.v, p.r) for p in particles if p.cc % 2 and p.xmix > 0.] + \
                [(p.m, p.v - self.v, p.r - self.r) for p in particles if not p.cc % 2])

        self.aei = self.calculate_orbital_elements(self.r-star.r, self.v)

        self.E = -self.mu/(2.*self.aei['a'])
        self.h = math.sqrt(self.mu*self.aei['a']*(1.-self.aei['e']**2.))
        self.specific_e = -0.5*(star.m*self.m)/self.aei['a']
        self.specific_h = self.m*star.m/(self.m+star.m)*math.sqrt(self.aei['a']*(1.-self.aei['e']**2.)*(star.m+self.m))

        self.gas_mass_fraction = self.m_gas/self.m
        self.core_mass_fraction = self.m_core/self.m#1.-self.gas_mass_fraction

        self.gas_binding_energy = self.calculate_binding_energy_old([p for p in particles if p.cc % 2], star, gas = 1)
        self.binding_energy = self.calculate_binding_energy_old(particles, star, gas = 0)
        self.energy = self.calculate_energy(particles)

        #self.quickplot(temp, self.aei, self.r)




    # calculate binding energy taking into account the other planet (p)
    def recalculate_binding_energy(self, planet, particles, star):

        self.binding_energy = self.calculate_binding_energy(particles, star, planet)

        print 'new binding energy: ', self.binding_energy['gas']['e'], self.binding_energy['tot']['e']
        print 'new gas binding energy all components: ', self.binding_energy['gas']
        print 'new tot binding energy all components: ', self.binding_energy['tot']





    def create_profile(self, particles):

        for particle in particles: particle.get_T()
        for particle in particles: particle.get_p()
        for p in particles: p.calculate_angular_momentum(self.v, self.r, self.j_unit_vector_gas)


        self.profile = PlanetProfile(particles, self.r)




    def quickplot(self, aei1, aei2, r):

        N = 1000
        theta = [i*2.*math.pi/N for i in range(N)]
        c = ['r', 'b']
        for j, aei in enumerate([aei1, aei2]):
            x,y,z,u,v,w = mercury_routines.convert_el2x(aei, self.mu)
            plt.plot(x, y, color = c[j], marker = 'o', markersize = 4., linestyle = '')
        plt.plot(r[0]*0.00465, r[1]*0.00465, color='k', marker = 'o', markersize = 1., linestyle = '')
        plt.show()




    def calculate_velocity(self, particles):

        if self.m == 0: return -1
        vx = sum([particle.m*particle.vx for particle in particles])/self.m
        vy = sum([particle.m*particle.vy for particle in particles])/self.m
        vz = sum([particle.m*particle.vz for particle in particles])/self.m

        return np.array([vx,vy,vz])




    def calculate_angular_momentum(self, mvr):

        g   = 6.67390*math.pow(10,-8)                               # gravitational constant [cgs]
        rs  = 6.9599*math.pow(10,10)                                # unit of length [cm]
        ms  = 1.9891*math.pow(10,33)                                # unit of mass [g]
        ts  = math.sqrt(math.pow(rs,3)/(g*ms))                      # unit of time [s]

        if len(mvr) > 0:
            random.shuffle(mvr)
            m, v, r = zip(*mvr)
            n = len(v)/2

            # calculate spin axis by taking cross products of pairs of particle velocities
            spin_vectors = np.array([np.cross(np.array(v1), np.array(v2), axis = 0) for v1, v2 in zip(v[0:n], v[n:2*n])])
            spin_vector = np.sum(vector/np.linalg.norm(vector) for vector in spin_vectors)
            spin_vector = spin_vector/np.linalg.norm(spin_vector)

            # calculate the angular momentum by projecting onto the spin axis
            d_proj = np.array([np.cross(_r - (self.r + spin_vector), _r - self.r)/np.linalg.norm(spin_vector) for _r in r])
            #d_proj = np.array([math.fabs(np.cross(np.array(_r), np.array(_r) - spin_vector))/np.linalg.norm(spin_vector) for _r in r])

            e_rot = 0.5*np.sum([_m*np.linalg.norm(np.cross(d, _v))**2./np.linalg.norm(d)**2. for _m, d, _v in zip(m, d_proj, v)])*ms*rs*rs/ts/ts

            l = np.sum([_m*np.linalg.norm(np.cross(d,_v)) for _m, _v, d in zip(m, v, d_proj)])*ms*rs*rs/ts

            return l, l/np.sum(m)/ms, e_rot, e_rot/np.sum(m)/ms, spin_vector

        else:
            return 0., 0., 0., 0., 0.




    def calculate_center_of_mass(self, particles):

        if self.m == 0: return -1
        x = sum([particle.m*particle.x for particle in particles])/self.m
        y = sum([particle.m*particle.y for particle in particles])/self.m
        z = sum([particle.m*particle.z for particle in particles])/self.m

        return np.array([x,y,z])




    def add_particle(self, particle):

        for i in range(3):
            self.r[i] = (self.r[i]*self.m+particle.r[i]*particle.m)/(self.m+particle.m)
            self.v[i] = (self.v[i]*self.m+particle.v[i]*particle.m)/(self.m+particle.m)
        self.m = self.m+particle.m
        self.particles.append(particle.i)




    # in:  r,v
    # out: a,e,i,g,n,M
    # adapted from mco_x2el subroutine in Mercury (Chambers 1999)
    def calculate_orbital_elements(self, r0, v0):

        if self.m == 0: return -1

        gm = self.mu
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
                bige = math.log(ce+sqrt(ce*ce-1.))
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




    def add_core(self, core):

        self.core = core




    def assign_name(self, name):

        self.name = name




    def calculate_energy(self, particles):

        k = np.sum([p.m*np.linalg.norm(p.v)**2./2. for p in particles])
        w = np.sum([p.m*p.grpot for p in particles])
        u = np.sum([p.m*p.u for p in particles])

        return {'k': k, 'w': w, 'u': u}




    def calculate_binding_energy(self, particles, star, planet):

        G = 1.

        dv_planet = [np.linalg.norm(p.v - self.v) for p in particles]
        dr_star = [np.linalg.norm(p.r - star.r) for p in particles]
        dr_planet2 = [np.linalg.norm(p.r - planet.r) for p in particles]

        k = [p.m*dv**2./2. for dv, p in zip(dv_planet, particles)]
        w = [p.m*(p.grpot + G*star.m/r + G*planet.m/r2) for r, r2, p in zip(dr_star, dr_planet2, particles)]
        u = [p.m*p.u for p in particles]

        e = [_k + _w + _u for _k, _w, _u in zip(k, w, u)]

        def f(x):
            return np.sum([_x*(1.-p.xmix) for _x, p in zip(x, particles) if p.cc % 2])

        return {'gas': {'k': f(k), 'w': f(w), 'u': f(u), 'e': f(e)}, \
                'tot': {'k': np.sum(k), 'w': np.sum(w), 'u': np.sum(u), 'e': np.sum(e)}}




    def calculate_binding_energy_old(self, particles, star, gas):

        G = 1.

        dv_planet = [np.linalg.norm(particle.v - self.v) for particle in particles]
        dr_star = [np.linalg.norm(particle.r - star.r) for particle in particles]

        e_binding = [particle.m*(dv**2./2.+particle.grpot+G*star.m/r) for dv, r, particle in zip(dv_planet, dr_star, particles)]

        if gas: return np.sum([e*(1.-particle.xmix) for e, particle in zip(e_binding, particles) if particle.cc % 2])
        return np.sum(e_binding)




    #in:  r,v
    #out: a,e,i,g,n,M
    def calculate_orbital_elements_old(self, r, v):
        if self.m == 0: return -1

        eps = math.pow(10.,-6)

        h=np.cross(r,v)
        norm=np.cross([0,0,1],h)

        evec = ((np.linalg.norm(v)**2.-self.mu/np.linalg.norm(r))*r-np.dot(r,v)*v)/self.mu
        e = np.linalg.norm(evec)

        energy = np.linalg.norm(v)**2./2.-self.mu/np.linalg.norm(r)

        if math.fabs(e-1.0)>eps:
            a = -self.mu/(2.*energy)
            p = a*(1.-e**2.)
        else:
            a = math.inf
            p = np.linalg.norm(h)**2./self.mu

        i = math.acos(h[2]/np.linalg.norm(h))

        n = math.acos(norm[0]/np.linalg.norm(norm))
        if norm[1]<0:
            n = 2.*math.pi-n

        g = math.acos(np.dot(norm,evec)/(np.linalg.norm(norm)*e))
        if evec[2]<0.0:
            g = 2.*math.pi-g

        M = math.acos(np.dot(evec,r)/(e*np.linalg.norm(r)))
        if np.dot(r,v)<0.0:
            M = 2.*math.pi-M

        # Convert semi-major axis to AU
        a = a*0.00465

        return {'a': a, 'e': e, 'i': i, 'g': g, 'n': n, 'M': M}



class PlanetProf:

    def __init__(self,Planet):

        nbins = 100
        self.r,self.rho,self.m,self.T,self.P,self.u,self.gamm,self.cc = Calcs.Generate_Profile(nbins,Planet.com,Planet.Particles)




























