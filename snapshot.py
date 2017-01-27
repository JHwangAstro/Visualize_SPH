from header import Header
from particle import Particle
from planet import Planet
import numpy as np

class Snapshot:

    def __init__(self, n, p):

        self.n = n # snapshot number
        self.p = n # parent snapshot number
        self.t = -1



    def initialize(self, args, parent):

        header, particles = self.read_snapshot(args.sdir, self.n)
        self.t = header.time
        self.create_system(args.num_p2, particles, parent)
        # If there are two planets, track the mutual energy
        if np.min([planet.m for planet in self.planets]) > 0.:
            self.planet_energy = self.calculate_planet_energy(self.planets)
        #finish this
        #plot_shapshot(idir,i,Planets,Envelope,Ejecta,Star)



    #read_snapshot imports the data from each snapshot and
    #outputs a list of particle objects and a header object
    #in:  sdir,ss
    #out: header,particles
    def read_snapshot(self, sdir, ss):

        # Set snapshot filename
        fname = sdir + 'out000' + repr(ss) + '.sph'
        if ss >= 10:   fname = sdir + 'out00' + repr(ss) + '.sph'
        if ss >= 100:  fname = sdir + 'out0'  + repr(ss) + '.sph'
        if ss >= 1000: fname = sdir + 'out'   + repr(ss) + '.sph'

        # Read snapshot
        f = open(fname, 'rb')
        header = Header(f)
        particles = [Particle(i,f) for i in range(header.ntot)]
        f.close()

        return header, particles



    def create_system(self, n_other_planets, particles, parent):

        # If no parent, iterate through each particle assuming the most massive
        # point-mass particle is the star, the others are other planets,
        # half of the non-point particles belong to planet 1,
        # and the remaining belong to planet 2.
        if parent.t == self.t:
            star = 0
            self.other_planets = []
            for particle in particles:
                if particle.u <= 0.:
                    if not star:
                        self.star = particle
                    else:
                        if particle.m > self.star.m: self.star = particle
                    particles.remove(particle)
                    self.other_planets.append(particle)
            self.other_planets.remove(self.star)

            n = len(particles) - n_other_planets
            self.planets = [Planet(particles[0:n/2], self.star), Planet(particles[n/2:n], self.star)]
            self.gas = []
        # If there is a parent, use the parent's particle bounds to update
        # particle bound properties.
        else:
            self.planets, self.gas = self.assign_particles(particles, parent.planets, parent.star)
            self.star = parent.star
            self.other_planets = parent.other_planets

        # Create core objects
        for planet in self.planets:
            planet.add_core(Planet([particles[i] for i in planet.particles if particles[i].cc % 1 == 0], self.star))
            print planet.gas_mass_fraction, planet.core_mass_fraction



    def assign_particles(self, particles, planets, star):

        # Update planets with new particle properties
        planets = [Planet([particles[i] for i in planet.particles], star) for planet in planets]

        for particle in particles:
            particle.determine_bound(planets, star)

        # Return updated objects
        bound_particles = [[particle for particle in particles if particle.bound == i] for i in range(len(planets))]
        planets  = [Planet(bound_particles[i], star) for i in range(len(planets))]
        free_particles = [particle for particle in particles if particle.bound > 1]

        return [planets, free_particles]



    def calculate_planet_energy(self, planets):

        G = 1.0
        n_planets = len(planets)

        v_cm = np.array([(planets[0].m*planets[0].v[i]+planets[1].m*planets[1].v[i])/(planets[0].m+planets[1].m) for i in range(3)])
        v_cm = [np.array([planet.v[i]-v_cm[i] for i in range(3)]) for planet in planets]
        v_cm2 = [np.sum(np.array([v_cm[i][j]**2. for j in range(3)])) for i in range(n_planets)]
        r = np.sqrt(np.sum(np.array([(planets[1].r[i]-planets[0].r[i])**2. for i in range(3)])))
        K_cm = 0.5*np.sum([planets[i].m*v_cm2[i] for i in range(n_planets)])
        P_cm = -G*planets[0].m*planets[1].m/r

        return K_cm+P_cm



    def set_unitless_time(self, innermost_period):

        self.unitless_time = self.t/innermost_period



