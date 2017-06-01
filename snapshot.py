
from datetime import datetime
import math
import numpy.linalg as LA
import numpy as np
from scipy.spatial.distance import cdist
import string
import struct

from header import Header
from planet import Planet
from particle import Particle
import plot_sph

import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline


# TODO
# Change updatecore to assign core particles as having a certain mass or x
# add in restartrad as final snapshot

# Notes: Assumes that non-hydro particles are stored at the end of the output files


class Snapshot:

    def __init__(self, n):

        self.n = n # snapshot number



    def initialize(self, fname, args, parent, particle_params):

        header, particles = self.read_snapshot(fname, particle_params)
        self.t = header.time

        self.star, self.star_index, self.other_planets, self.planets, self.bound_mass, self.accreted_mass, self.ejected_mass \
                = self.create_system(particles, parent)

        # float containing orbital energy
        self.orbital_energy = self.calculate_orbital_energy(particles)

        print self.orbital_energy

        # recalculate binding energies now that we have both planets
        p, p2 = self.planets
        p.recalculate_binding_energy(p2, [particles[k] for k in p.particles], self.star)
        p2.recalculate_binding_energy(p, [particles[k] for k in p2.particles], self.star)
        #for p, p2 in zip(self.planets, self.planets[::-1]):
        #    p.recalculate_binding_energy(p2, [particles[k] for k in p.particles], self.star)

        if args.rewrite: self.write_snapshot(header, [planet.core for planet in self.planets], particles, args.sdir)
        if np.min([planet.m for planet in self.planets]) > 0.: self.planet_energy = self.calculate_planet_energy(self.planets)
#        self.plot_snapshot(args.sdir, args.idir, particles)




    def create_planet_profiles(self, fname, particle_params):

        header, particles = self.read_snapshot(fname, particle_params)

        for planet in self.planets:
            planet.create_profile([particles[k] for k in planet.particles])




    #read_snapshot imports the data from each snapshot and
    #outputs a list of particle objects and a header object
    #in:  particle_params
    #out: header, particles
    def read_snapshot(self, fname, particle_params):

        # Read snapshot
        f = open(fname, 'rb')
        header = Header(f)
        dtype = np.dtype("i4,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8, \
                f8,f8,f8,f8,f8,f8,i4,f8,i4")
        data = np.fromfile(f, dtype = dtype, count = header.ntot)
        f.close()

        # Create hydro-particle objects
        particles = {k: Particle(k, data[k], particle_params) for k, particle_params in particle_params.iteritems()}

        # Fill in non-hydro particles with dummy parameters
        # This assumes all non-hydro particles are at the end of the file
        for i in range(len(particle_params), header.ntot):
            particles[i] = Particle(i, data[i], None)

        return header, particles



    def create_system(self, particles, parent):

        # If no parent, iterate through each particle assuming the most massive
        # point-mass particle is the star, the others are other planets,
        # half of the non-point particles belong to planet 1,
        # and the remaining belong to planet 2.

        # Checks if parent is itself
        if parent.t == self.t:

            star = 0
            other_planets = {}
            for k, particle in particles.iteritems():
                if particle.u <= 0.:
                    if not star:
                        star = particle
                        star_index = k
                    else:
                        if particle.m > star.m:
                            star = particle
                            star_index = k
                    other_planets[k] = particle

            # Remove point-particles from particles
            for k in other_planets:
                del particles[k]
            del other_planets[star_index]

            n = len(particles)# - n_other_planets
            planets = [Planet([particles[k] for k in range(0,n/2)], star), \
                    Planet([particles[k] for k in range(n/2,n)], star)]

            bound_mass = 0.
            accreted_mass = 0.
            ejected_mass = 0.

        # If there is a parent, use the parent's particle bounds to update
        # particle bound properties.
        else:
            star = parent.star
            star_index = parent.star_index
            other_planets = parent.other_planets

            for k in other_planets:
                particles.pop(k, None)
            particles.pop(star_index, None)

            planets, bound_mass, accreted_mass, ejected_mass = self.assign_particles(particles, parent.planets, star)

        # Create core objects and assign names
        for j, planet in enumerate(planets):
            planet.assign_name(string.ascii_lowercase[j+1])
            planet.add_core(Planet([particles[k] for k in planet.particles if particles[k].core], star))

        return star, star_index, other_planets, planets, bound_mass, accreted_mass, ejected_mass




    def assign_particles(self, particles, planets, star):

        # Create planet cores with new particle properties
        planets = [Planet([particles[i] for i in planet.particles if particles[i].core2], star) for planet in planets]

        # Physical constants
        G  = 1.0                                                # Gravitational constant in code units
        rs = 1.0                                                # Radius of star in code units
        ms = 1.98855*math.pow(10.,33.)                          # Solar mass in g
        me = 5.972*math.pow(10.,27.)                            # Earth mass in g
        au_to_rs = 1./.00465
        time_to_years  = 5.05000318388e-05                      # Convert from code units to years

        # calculate the distances, velocities and energies from each particle to the respective deep cores
        ind = np.array([particle.i  for k, particle in particles.iteritems() if not particle.core2])
        m   = np.array([particle.m  for k, particle in particles.iteritems() if not particle.core2])
        x   = np.array([particle.x-star.x  for k, particle in particles.iteritems() if not particle.core2])
        y   = np.array([particle.y-star.y  for k, particle in particles.iteritems() if not particle.core2])
        z   = np.array([particle.z-star.z  for k, particle in particles.iteritems() if not particle.core2])
        vx  = np.array([particle.vx for k, particle in particles.iteritems() if not particle.core2])
        vy  = np.array([particle.vy for k, particle in particles.iteritems() if not particle.core2])
        vz  = np.array([particle.vz for k, particle in particles.iteritems() if not particle.core2])
        grpot = np.array([particle.grpot for k, particle in particles.iteritems() if not particle.core2])

        # For each planet, calculate if the particle is bound in ascending distance
        xyz  = np.array([particle.r-star.r for k, particle in particles.iteritems() if not particle.core2])
        vxyz = np.array([particle.v for k, particle in particles.iteritems() if not particle.core2])

        # calculate the distances, velocities and energies of each particle
        r = (x**2 + y**2 + z**2)**0.5
        v = (vx**2 + vy**2 + vz**2)**0.5
        Ecom = m*(v**2/2. + grpot)

        # determine the state of the particle relative to the entire system
        ejecta = {particle.i for k, particle in particles.iteritems()}
        bound  = set([i for i in ind[Ecom < 0.]])
        for planet in planets: bound = bound | set([i for i in planet.particles])
        accreted = set([i for i in ind[r < rs]])
        ejecta = ejecta - bound - accreted
        bound = bound - accreted

        E = {}
        for i, planet in enumerate(planets):
            x_unit = np.linalg.norm(planet.r-star.r)
            theta = -math.atan2(planet.r[1] - star.r[1], planet.r[0] - star.r[0])
            Rz = np.array([[math.cos(theta), -math.sin(theta), 0.], [math.sin(theta), math.cos(theta), 0.], [0., 0., 1.]])
            planet_rot = np.dot(Rz, (planet.r-star.r).transpose())
            phi = -math.atan2(planet_rot[2], planet_rot[0])
            Ry = np.array([[math.cos(phi), 0., -math.sin(phi)], [0., 1., 0.], [math.sin(phi), 0., math.cos(phi)]])
            planet_rot = np.dot(Ry, planet_rot)/x_unit
            xyz_rot = np.array(np.dot(np.dot(Ry, Rz), xyz.transpose())/x_unit).transpose()

            dx = np.array([x[0] for x in cdist(xyz_rot, [planet_rot])])
            dv = np.array([x[0] for x in cdist(vxyz, [planet.v])])/x_unit/time_to_years*planet.aei['a']**1.5
            v2 = np.array([np.linalg.norm(x)**2. for x in dv])

            # in order of distance from the core, determine if a particle is bound using the Jacobi Constant
            for j in np.argsort(dx):
                mu1 = star.m/(star.m+planet.m)
                mu2 = planet.m/(star.m+planet.m)
                #CJ_L1 = 3.+3.**(4./3.)*mu2**(2./3.)-10.*mu2/3.
                CJ_L1 = 3.+3.**(4./3.)*mu2**(2./3.)-10.*mu2/3.
                CJ = (xyz_rot[j,0]**2.+xyz_rot[j,1]**2.)+2.*(mu1/r[j]*x_unit+mu2/dx[j])-v2[j]
                if CJ > CJ_L1: planet.add_particle(particles[ind[j]])
                E[(i,ind[j])] = m[j]*(dv[j]**2./2.0-G*planet.m/dx[j])

        if len(planets) == 1:
            bound_to_planet = set([i for i in planets[0].particles])
            bound = bound - bound_to_planet[0] - accreted
            bound_to_planet = [bound_to_planet - accreted]


        if len(planets) == 2:
            bound_to_planet = [set([i for i in planet.particles]) for planet in planets]
            intersection = bound_to_planet[0] & bound_to_planet[1]
            e21 = set([i for i in intersection if E[(0,i)] > E[(1,i)]])
            bound_to_planet[0] = bound_to_planet[0] - e21
            bound_to_planet[1] = bound_to_planet[1] - (intersection - e21)
            bound = bound - bound_to_planet[0] - bound_to_planet[1] - accreted
            bound_to_planet = [bound_to_planet[i] - accreted for i in range(2)]

        print 'n particles bound: ', len(bound_to_planet[0]), len(bound_to_planet[1]), len(bound), len(accreted), len(ejecta)
        print 'n particles: ', len(bound_to_planet[0])+len(bound_to_planet[1])+len(bound)+len(accreted)+len(ejecta)

        # Return updated objects
        planets = [Planet([particles[j] for j in bound_to_planet[i]], star) for i in range(len(planets))]

        # report masses in units of me
        bound_mass = np.sum([particles[k].m for k in list(bound)])*ms/me
        accreted_mass = (np.sum([particles[k].m for k in list(accreted)])-star.m)*ms/me
        ejected_mass = np.sum([particles[k].m for k in list(ejecta)])*ms/me

        print 'masses: ', self.n, [planet.m*ms/me for planet in planets], bound_mass, accreted_mass, ejected_mass

        return planets, bound_mass, accreted_mass, ejected_mass




    def calculate_planet_energy(self, planets):

        G = 1.0
        n_planets = len(planets)

        v_cm = [(planets[0].m*planets[0].v[i]+planets[1].m*planets[1].v[i])/(planets[0].m+planets[1].m) for i in range(3)]
        v_cm = [[planet.v[i]-v_cm[i] for i in range(3)] for planet in planets]
        v_cm2 = [sum([v_cm[i][j]**2. for j in range(3)]) for i in range(n_planets)]
        r = sum([(planets[1].r[i]-planets[0].r[i])**2. for i in range(3)])**0.5
        K_cm = 0.5*sum([planets[i].m*v_cm2[i] for i in range(n_planets)])
        P_cm = -G*planets[0].m*planets[1].m/r

        return K_cm+P_cm




    def write_snapshot(self, header, planets, particles, sdir):

        for i, planet in enumerate(planets):

            # Set snapshot filename
            fname = sdir + string.ascii_lowercase[i] + '_new_out' + format(self.n, "04") + '.sph'

            #_______________________________________________________________WRITES DATA FILE(S)
            f = open(fname, 'wb')
            dum = 124
            f.write(struct.pack('i',dum))                       #dummy variable
            f.write(struct.pack('i',header.ntot))               #total number of particles
            f.write(struct.pack('i',header.nnopt))              #optimal number of neighbors
            f.write(struct.pack('d',header.hmin))               #smoothing length for black holes
            f.write(struct.pack('d',header.hmax))               #smoothing length for black holes
            f.write(struct.pack('d',header.sep0))               #initial separation of two stars
            f.write(struct.pack('d',header.tf))                 #time the code stops in code units
            f.write(struct.pack('d',header.dtout))              #time interval between writing outxxx.sph files
            f.write(struct.pack('i',header.nout))               #xxx value of outputxxx.sph
            f.write(struct.pack('i',header.nit))                #number of iterations completed
            f.write(struct.pack('d',header.time))               #current time in code units
            f.write(struct.pack('i',header.nav))                #chooses artificial viscosity scheme
            f.write(struct.pack('d',header.alpha))              #artificial viscosity parameter
            f.write(struct.pack('d',header.beta))               #artificial viscosity parameter
            f.write(struct.pack('d',header.tskip))              #??
            f.write(struct.pack('i',header.ngr))                #gravity boundary conditions
            f.write(struct.pack('i',header.nrelax))             #relaxation flag
            f.write(struct.pack('d',header.trelax))             #timescale for artificial drag force
            f.write(struct.pack('d',header.dt))                 #current timestep
            f.write(struct.pack('d',header.omega2))             #square of the orbital velocity
            f.write(struct.pack('i',dum))                       #dummy variable

            dum = 140
            for k, particle in particles.iteritems():
                x1 = (particle.x-planet.r[0])
                y1 = (particle.y-planet.r[1])
                z1 = (particle.z-planet.r[2])
                f.write(struct.pack('i',dum))                     #dummy variable
                f.write(struct.pack('d',x1))                    #x position
                f.write(struct.pack('d',y1))                    #y position
                f.write(struct.pack('d',z1))                    #z position
                f.write(struct.pack('d',particle.m))          #mass of the particle
                f.write(struct.pack('d',particle.h))          #smoothing length
                f.write(struct.pack('d',particle.rho))        #density
                f.write(struct.pack('d',particle.vx))         #x-velocity
                f.write(struct.pack('d',particle.vy))         #y-velocity
                f.write(struct.pack('d',particle.vz))         #z-velocity
                f.write(struct.pack('d',particle.vxdot))      #x-acceleration
                f.write(struct.pack('d',particle.vydot))      #y-acceleration
                f.write(struct.pack('d',particle.vzdot))      #z-acceleration
                f.write(struct.pack('d',particle.u))          #potential
                f.write(struct.pack('d',particle.udot))       #change in potential
                f.write(struct.pack('d',particle.grpot))      #gravitational potential
                f.write(struct.pack('d',particle.mmw))        #mean molecular weight
                f.write(struct.pack('i',particle.cc))         #flag for particle type
                f.write(struct.pack('d',particle.divv))       #divergence of velocities
                f.write(struct.pack('i',dum))                     #dummy variable

            f.close()




    def set_unitless_time(self, innermost_period):

        print 'time: ', self.t/innermost_period
        return self.t/innermost_period



    def set_hill_radius(self, hill_radius):

        return hill_radius



    def set_outer_orbit(self, outer_orbit):

        return outer_orbit



    def plot_snapshot(self, sdir, idir, particles):

        n = 1000

        u  = np.array([particle.u for k, particle in particles.iteritems()])
        u0 = np.array([particle.u0 for k, particle in particles.iteritems()])
        du = u[0:len(u0)]/u0

        p, x = np.histogram(du, bins=n)
        x = np.array([(x[i]+x[i+1])/2. for i in range(len(x)-1)])

        plt.plot(x, p)
        fname = idir + 'u_pdf' + format(self.n, "04") + '.png'
        plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
        plt.clf()

        # for each planet, create a pdf of the internal energy of
        # gas particles scaled to the initial
        for planet in self.planets:

            u = np.array([particles[i].u for i in planet.particles if not particles[i].core])
            u0 = np.array([particles[i].u0 for i in planet.particles if not particles[i].core])

            rho = np.array([particles[i].rho for i in planet.particles if not particles[i].core])

            du = u/u0
            p, x = np.histogram(du, bins=n)
            x = np.array([(x[i]+x[i+1])/2. for i in range(len(x)-1)])

            plt.plot(x, p)
            fname = idir + planet.name + '_ugas_pdf' + format(self.n, "04") + '.png'
            plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
            plt.clf()

            T = np.array([particles[i].T for i in planet.particles if not particles[i].core])

            p, x = np.histogram(np.log10(T), bins=n)
            x = np.array([(x[i]+x[i+1])/2. for i in range(len(x)-1)])

            plt.plot(x, p)
            fname = idir + planet.name + '_T_pdf' + format(self.n, "04") + '.png'
            plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
            plt.clf()

            t_cool = np.array([particles[i].t_cool for i in planet.particles if not particles[i].core])

            p, x = np.histogram(np.log10(t_cool), bins=n)
            x = np.array([(x[i]+x[i+1])/2. for i in range(len(x)-1)])

            plt.plot(x, p)
            fname = idir + planet.name + '_t_cool_pdf' + format(self.n, "04") + '.png'
            plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
            plt.clf()

            heatmap, xedges, yedges = np.histogram2d(np.log10(T), np.log10(t_cool), bins = 100)
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

            plt.imshow(heatmap.T, extent = extent, origin = 'lower', aspect = 'auto', cmap = 'afmhot')
            fname = idir + planet.name + '_T_tcool_heatmap' + format(self.n, "04") + '.png'
            plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
            plt.clf()

            heatmap, xedges, yedges = np.histogram2d(np.log10(u), np.log10(t_cool), bins = 100)
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

            plt.imshow(heatmap.T, extent = extent, origin = 'lower', aspect = 'auto', cmap = 'afmhot')
            fname = idir + planet.name + '_u_tcool_heatmap' + format(self.n, "04") + '.png'
            plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
            plt.clf()


            heatmap, xedges, yedges = np.histogram2d(np.log10(rho), np.log10(t_cool), bins = 100)
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

            plt.imshow(heatmap.T, extent = extent, origin = 'lower', aspect = 'auto', cmap = 'afmhot')
            fname = idir + planet.name + '_rho_tcool_heatmap' + format(self.n, "04") + '.png'
            plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
            plt.clf()


            heatmap, xedges, yedges = np.histogram2d(np.log10(rho), np.log10(T), bins = 100)
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

            plt.imshow(heatmap.T, extent = extent, origin = 'lower', aspect = 'auto', cmap = 'afmhot')
            fname = idir + planet.name + '_rho_T_heatmap' + format(self.n, "04") + '.png'
            plt.savefig(fname, facecolor='w', edgecolor='w', format='png',dpi=200)
            plt.clf()




    def calculate_orbital_energy(self, particles):

        return np.sum([p.m*(LA.norm(p.v)**2./2.-self.star.m/LA.norm(p.r-self.star.r)) for k, p in particles.iteritems()])# if p.i != self.star_index])







