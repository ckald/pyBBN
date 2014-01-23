import numpy
from common import momentum_to_index, GRID, STATISTICS, PARAMS, theta, benchmark
from ds import D
from scipy import integrate
import numericalunits as nu


class INTERACTIONS:
    DECAY = 'decay'


class Interaction:

    particles = []
    in_particles = []
    out_particles = []
    decoupling_temperature = 0.
    symmetry_factor = 1.
    constant = 1/64 / numpy.pi**3
    K1 = 0
    K2 = 0

    def __init__(self, *args, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.particles = self.in_particles + self.out_particles
        self.collision = numpy.vectorize(self.collision)

    def calculate(self):
        if PARAMS.T < self.decoupling_temperature and not self.in_particles[0].in_equilibrium:
            print ".",
            return 0

        self.in_particles[0].distribution += PARAMS.dt * self.collision(GRID.TEMPLATE)

        print
        print self.in_particles[0].density

    def collision(self, p0):
        tmp = self.constant / p0 / self.particles[0].energy(p0) * self.symmetry_factor

        tmp *= integrate.dblquad(
            lambda p1, p2: self.integrand(p=[p0, p1, p2, p0 - p1 - p2]),
            GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
            lambda x: GRID.MIN_MOMENTUM, lambda x: GRID.MAX_MOMENTUM,
        )[0]
        return tmp

    def integrand(self, p=[]):
        E = []
        m = []
        for i, particle in enumerate(self.particles):
            E.append(particle.energy(p[i]))
            m.append(particle.mass)

        integrand = theta(E[0] - E[1] - E[2] > m[3])
        if integrand != 0.:
            d = D(p=p, E=E, m=m, K1=self.K1, K2=self.K2)
            integrand *= d
            if integrand != 0.:
                f = self.F(in_p=p[:1], out_p=p[1:])
                integrand *= (
                    f * d * p[1] / E[1] * p[2] / E[2]
                )

        return integrand

    def F(self, in_p=[], out_p=[]):

        def mult_them(out_particles, in_particles, out_p, in_p):
            temp = 1.
            for i, particle in enumerate(out_particles):
                temp *= particle.dist_value(out_p[i])
            for i, particle in enumerate(in_particles):
                temp *= 1. - particle.eta * particle.dist_value(in_p[i])

            return temp

        f = mult_them(self.out_particles, self.in_particles, out_p, in_p)\
            - mult_them(self.in_particles, self.out_particles, in_p, out_p)

        return f
