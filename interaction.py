# -*- coding: utf-8 -*-
import numpy
from common import GRID, PARAMS, theta, UNITS
from ds import D
from scipy import integrate


class INTERACTIONS:
    DECAY = 'decay'


class Interaction:

    """
        \begin{multline}
            I_{coll}(t,p_i) = \frac{1}{2 E_i} \sum_{reactions} \int \cdots \int
                \prod_{j} \frac{d^3 p_{in,j}}{(2 \pi)^3 2 E_{in,j}}
                \prod_{k} \frac{d^3 p_{out,k}}{(2 \pi)^3 2 E_{out,k}} \\
                S |\mathcal{M}|^2 \mathcal{F}(\{f_\alpha\}) (2 \pi)^4 \\
                \delta^4(\sum_{m} p_{in,m}^\mu - \sum_{n} p_{out,n}^\mu)
        \end{multline}
    """

    particles = []
    in_particles = []
    out_particles = []
    decoupling_temperature = 0.
    symmetry_factor = 1.
    constant = 1./64. / numpy.pi**3
    K1 = 0
    K2 = 0
    s = 1

    def __init__(self, *args, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.particles = self.in_particles + self.out_particles
        self.collision = numpy.vectorize(self.collision, otypes=[numpy.float_])
        self.s = 1 if len(self.in_particles) == 2 else -1

    def calculate(self):
        if PARAMS.T > self.decoupling_temperature and not self.in_particles[0].in_equilibrium:

            self.particles[0].F_f.append(
                lambda p0, p1, p2:
                self.quad_integrand(p=[p0, p1, p2, 0], method=self.variable_integrand)
            )
            self.particles[0].F_1.append(
                lambda p0, p1, p2:
                self.quad_integrand(p=[p0, p1, p2, 0], method=self.const_integrand)
            )

    def collision(self, p0):
        tmp = self.constant / p0 / self.particles[0].energy_normalized(p0)

        integral = integrate.dblquad(
            lambda p1, p2: (
                self.particles[0].distribution(p0)
                * self.quad_integrand(p=[p0, p1, p2, 0], method=self.variable_integrand)
                + self.quad_integrand(p=[p0, p1, p2, 0], method=self.const_integrand)
            ),
            GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
            lambda x: GRID.MIN_MOMENTUM, lambda x: GRID.MAX_MOMENTUM
        )
        tmp *= integral[0] * PARAMS.m**5 / PARAMS.x**6 / PARAMS.H
        # print p0 / UNITS.MeV, '\t', \
        #     tmp * PARAMS.dx, '\t',\
        #     integral[1] / integral[0]

        return tmp

    def calculate_impulses(self, p=[]):
        E = []
        m = []
        for i, particle in enumerate(self.particles):
            E.append(particle.energy_normalized(p[i]))
            m.append(particle.mass_normalized)
        E[3] = E[0] + self.s * E[1] - E[2]
        p[3] = numpy.sqrt(numpy.abs(E[3]**2 - m[3]**2))
        return p, E, m

    def quad_integrand(self, p=[], method=None):
        p, E, m = self.calculate_impulses(p)
        integrand = theta(E[3] - m[3]) * self.symmetry_factor
        if integrand == 0:
            return integrand

        integrand *= D(p=p, E=E, m=m, K1=self.K1, K2=self.K2)
        if integrand == 0:
            return integrand

        return integrand * method(p=p, E=E, m=m)

    def variable_integrand(self, p=[], E=[], m=[]):
        f = \
            - self.in_particles[0].eta * \
            self.F_A(in_p=p[:len(self.in_particles)], out_p=p[len(self.in_particles):])\
            - self.F_B(in_p=p[:len(self.in_particles)], out_p=p[len(self.in_particles):])
        integrand = f * p[1] / E[1] * p[2] / E[2]

        return integrand

    def const_integrand(self, p=[], E=[], m=[]):
        f = self.F_A(in_p=p[:len(self.in_particles)], out_p=p[len(self.in_particles):])
        integrand = f * p[1] / E[1] * p[2] / E[2]

        return integrand

    # def F(self, in_p=[], out_p=[]):

    #     def mult_them(out_particles, in_particles, out_p, in_p):
    #         temp = 1.
    #         for i, particle in enumerate(out_particles):
    #             temp *= particle.distribution(out_p[i])
    #         for i, particle in enumerate(in_particles):
    #             temp *= 1. - particle.eta * particle.distribution(in_p[i])

    #         return temp

    #     f = mult_them(self.out_particles, self.in_particles, out_p, in_p)\
    #         - mult_them(self.in_particles, self.out_particles, in_p, out_p)

    #     return f

    # F == f (±A - B) + A
    def F_A(self, in_p=[], out_p=[], index=0):
        temp = 1.
        for i, particle in enumerate(self.out_particles):
            temp *= particle.distribution(out_p[i])
        for i, particle in enumerate(self.in_particles):
            if i != index:
                temp *= 1. - particle.eta * particle.distribution(in_p[i])
        return temp

    def F_B(self, in_p=[], out_p=[], index=0):
        temp = 1.
        for i, particle in enumerate(self.in_particles):
            if i != index:
                temp *= particle.distribution(in_p[i])
        for i, particle in enumerate(self.out_particles):
                temp *= 1. - particle.eta * particle.distribution(out_p[i])
        return temp
