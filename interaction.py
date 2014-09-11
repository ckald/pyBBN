# -*- coding: utf-8 -*-
import numpy
from common import GRID, PARAMS, theta
from ds import D
from scipy import integrate


class INTERACTIONS:
    DECAY = 'decay'


"""
== Boltzmann collision integral ==

\begin{multline}
    I_{coll}(t,p_i) = \frac{1}{2 E_i} \sum_{reactions} \int \cdots \int
        \prod_{j} \frac{d^3 p_{in,j}}{(2 \pi)^3 2 E_{in,j}}
        \prod_{k} \frac{d^3 p_{out,k}}{(2 \pi)^3 2 E_{out,k}} \\
        S |\mathcal{M}|^2 \mathcal{F}(\{f_\alpha\}) (2 \pi)^4 \\
        \delta^4(\sum_{m} p_{in,m}^\mu - \sum_{n} p_{out,n}^\mu)
\end{multline}
"""


class Interaction:

    """ Main class used for calculation of the non-equilibrium dynamics of the particles """

    # Incoming particles
    in_particles = []
    # Outgoing particles
    out_particles = []
    # All particles involved
    particles = []
    # Temperature when the typical interaction time exceeds the Hubble expansion time
    decoupling_temperature = 0.
    # Interaction symmetry factor
    symmetry_factor = 1.
    # Constant shared by all 4-particle interactions
    constant = 1./64. / numpy.pi**3

    """ Four-particle interactions of the interest can all be rewritten in a form

        \begin{equation}
            |\mathcal{M}|^2 = \sum_{\{i \neq j \neq k \neq l\}} K_1 (p_i \cdot p_j)\
                 + K_2 m_i m_j (p_k \cdot p_l)
        \end{equation} """
    K1 = 0.
    K2 = 0.
    """ 2-to-2 interactions and 1-to-3 decays can be generalized by introducing a multiplier `s` \
        for the momentum of the second particle: $p_1 + s p_2 = p_3 + p_4 $ """
    s = 1.

    def __init__(self, *args, **kwargs):
        """ Init """
        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.particles = self.in_particles + self.out_particles
        # self.collision = numpy.vectorize(self.collision, otypes=[numpy.float_])
        self.s = 1. if len(self.in_particles) == 2 else -1.

    def calculate(self):
        """
        Calculate collision integral constants and save them to the first involved particle
        """
        if PARAMS.T > self.decoupling_temperature and not self.in_particles[0].in_equilibrium:
            self.particles[0].F_f.append(
                lambda p0, p1, p2:
                self.quad_integrand(p=[p0, p1, p2, 0], method=self.variable_integrand)
            )
            self.particles[0].F_1.append(
                lambda p0, p1, p2:
                self.quad_integrand(p=[p0, p1, p2, 0], method=self.const_integrand)
            )

    # def collision(self, p0):
    #     """ TODO: what is this? Probably not used at the moment """
    #     tmp = self.constant / p0 / self.particles[0].energy_normalized(p0)

    #     integral = integrate.dblquad(
    #         lambda p1, p2: (
    #             self.particles[0].distribution(p0)
    #             * self.quad_integrand(p=[p0, p1, p2, 0], method=self.variable_integrand)
    #             + self.quad_integrand(p=[p0, p1, p2, 0], method=self.const_integrand)
    #         ),
    #         GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
    #         lambda x: GRID.MIN_MOMENTUM, lambda x: GRID.MAX_MOMENTUM
    #     )
    #     tmp *= integral[0] * PARAMS.m**5 / PARAMS.x**6 / PARAMS.H
    #     # print p0 / UNITS.MeV, '\t', \
    #     #     tmp * PARAMS.dx, '\t',\
    #     #     integral[1] / integral[0]

    #     return tmp

    def calculate_impulses(self, p=[]):
        """ Helper procedure that caches energies and normalized masses of particles """
        E = []
        m = []
        for i, particle in enumerate(self.particles):
            E.append(particle.energy_normalized(p[i]))
            m.append(particle.mass_normalized)
        E[3] = E[0] + self.s * E[1] - E[2]
        p[3] = numpy.sqrt(numpy.abs(E[3]**2 - m[3]**2))
        return p, E, m

    def quad_integrand(self, p=[], method=None):
        """ Total collision integral interior with performance optimizations """
        p, E, m = self.calculate_impulses(p)
        integrand = theta(E[3] - m[3]) * self.symmetry_factor
        if integrand == 0:
            return 0

        integrand *= D(p=p, E=E, m=m, K1=self.K1, K2=self.K2)
        if integrand == 0:
            return 0

        return -integrand * method(p=p, E=E, m=m) * p[1] / E[1] * p[2] / E[2]

    def variable_integrand(self, p=[], E=[], m=[]):
        """ Variable part of the collision integral ready for integration """
        return \
            self.F_A(in_p=p[:len(self.in_particles)],
                     out_p=p[len(self.in_particles):])\
            - self.in_particles[0].eta * self.F_B(in_p=p[:len(self.in_particles)],
                                                  out_p=p[len(self.in_particles):])

    def const_integrand(self, p=[], E=[], m=[]):
        """ Constant part of the collision integral ready for integration """
        return - self.F_B(in_p=p[:len(self.in_particles)], out_p=p[len(self.in_particles):])

    """
    \begin{equation}
        -\mathcal{F}(f) = f_1 f_2 (1 \pm f_3) (1 \pm f_4) - f_3 f_4 (1 \pm f_1) (1 \pm f_2)
    \end{equation}

    \begin{equation}
        -\mathcal{F}(f) = f_1 (f_2 (1 \pm f_3) (1 \pm f_4) \pm f_3 f_4 (1 \pm f_2)) \
        - f_3 f_4 (1 \pm f_2)
    \end{equation}

    \begin{equation}
        -\mathcal{F}(f) = f_1 (\mathcal{F}_A \pm \mathcal{F}_B) - \mathcal{F}_B
    \end{equation}

    \begin{equation}
        \mathcal{F}_A = f_2 (1 \pm f_3) (1 \pm f_4)
    \end{equation}
    \begin{equation}
        \mathcal{F}_B = f_3 f_4 (1 \pm f_2)
    \end{equation}

    """
    def F_A(self, in_p=[], out_p=[], index=0):
        """ Forward reaction distribution functional term without the `index`-th particle """
        temp = 1.
        for i, particle in enumerate(self.out_particles):
            if i != index:
                temp *= particle.distribution(out_p[i])
        for i, particle in enumerate(self.in_particles):
            temp *= 1. - particle.eta * particle.distribution(in_p[i])
        return temp

    def F_B(self, in_p=[], out_p=[], index=0):
        """ Backward reaction distribution functional term without the `index`-th particle  """
        temp = 1.
        for i, particle in enumerate(self.in_particles):
            temp *= particle.distribution(in_p[i])
        for i, particle in enumerate(self.out_particles):
            if i != index:
                temp *= 1. - particle.eta * particle.distribution(out_p[i])
        return temp
