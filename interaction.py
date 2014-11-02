# -*- coding: utf-8 -*-
import numpy
from common import PARAMS, theta
from ds import D, Db1, Db2


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

    integrals = []

    in_particles = []  # Incoming particles
    out_particles = []  # Outgoing particles
    particles = []  # All particles involved
    # Temperature when the typical interaction time exceeds the Hubble expansion time
    decoupling_temperature = 0.
    # Interaction symmetry factor
    symmetry_factor = 1.

    """ Four-particle interactions of the interest can all be rewritten in a form

        \begin{equation}
            |\mathcal{M}|^2 = \sum_{\{i \neq j \neq k \neq l\}} K_1 (p_i \cdot p_j) (p_k \cdot p_l)\
                 + K_2 m_i m_j (p_k \cdot p_l)
        \end{equation} """
    K1 = 0.
    K2 = 0.
    order = (0, 1, 2, 3)
    """ 2-to-2 interactions and 1-to-3 decays can be generalized by introducing a multiplier `s` \
        for the momentum of the second particle: $p_1 + s p_2 = p_3 + p_4 $ """
    s = 1.

    def __init__(self, *args, **kwargs):
        """ Init """

        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.integrals = []
        self.particles = self.in_particles + self.out_particles

        # Remember all particle species we've already considered to avoid double-counting
        accounted_particles = set()

        # Step 1. Permute all in_particles
        in_particles = self.in_particles
        out_particles = self.out_particles
        order = self.order
        accounted_particles = self.init_integrals(in_particles, out_particles,
                                                  order, accounted_particles)

        # Step 2. Turn to the backward process
        order = order[len(self.in_particles):] + order[:len(self.in_particles)]

        # Step 3. Permute all new in_particles (former out_particle)
        accounted_particles = self.init_integrals(out_particles, in_particles,
                                                  order, accounted_particles)

    def init_integrals(self, in_particles, out_particles, order, accounted_particles):
        for i, particle in enumerate(in_particles):

            if particle in accounted_particles:
                # Skip already accounted species
                continue

            # Add interaction integrals by putting each incoming particle as the first one
            self.integrals.append(Integral(
                in_particles=in_particles[i:i+1] + in_particles[:i] + in_particles[i+1:],
                out_particles=out_particles,
                decoupling_temperature=self.decoupling_temperature,
                K1=self.K1,
                K2=self.K2,
                order=order[i:i+1] + order[:i] + order[i+1:]
            ))
            accounted_particles.add(particle)
        return accounted_particles

    def calculate(self):
        """
        Calculate collision integral constants and save them to the first involved particle
        """

        for integral in self.integrals:
            integral.calculate()


class Integral:

    """ Main class used for calculation of the non-equilibrium dynamics of the particles """

    in_particles = []  # Incoming particles
    out_particles = []  # Outgoing particles
    particles = []  # All particles involved
    # Temperature when the typical interaction time exceeds the Hubble expansion time
    decoupling_temperature = 0.
    # Interaction symmetry factor
    symmetry_factor = 1.

    """ Four-particle interactions of the interest can all be rewritten in a form

        \begin{equation}
            |\mathcal{M}|^2 = \sum_{\{i \neq j \neq k \neq l\}} K_1 (p_i \cdot p_j) (p_k \cdot p_l)\
                 + K_2 m_i m_j (p_k \cdot p_l)
        \end{equation} """
    K1 = 0.
    K2 = 0.
    order = (0, 1, 2, 3)
    """ 2-to-2 interactions and 1-to-3 decays can be generalized by introducing a multiplier `s` \
        for the momentum of the second particle: $p_1 + s p_2 = p_3 + p_4 $ """
    s = 1.

    def __init__(self, *args, **kwargs):
        """ Init """
        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.particles = self.in_particles + self.out_particles
        self.s = 1. if len(self.in_particles) == 2 else -1.

    def __str__(self):
        return " + ".join([p.name for p in self.in_particles]) \
            + " ⟶  " + " + ".join([p.name for p in self.out_particles])

    def __repr__(self):
        return self.__str__()

    def calculate(self):
        """
        Calculate collision integral constants and save them to the first involved particle
        """
        if PARAMS.T > self.decoupling_temperature and not self.in_particles[0].in_equilibrium:

            self.particles[0].collision_integrands.append(
                lambda p0, p1, p2: self.integrand(p=[p0, p1, p2, 0], order=self.order)
            )

    def calculate_impulses(self, p=[]):
        """ Helper procedure that caches energies and normalized masses of particles """
        E = []
        m = []
        for i, particle in enumerate(self.particles):
            E.append(particle.energy_normalized(p[i]))
            m.append(particle.mass_normalized)
        E[3] = sum(self.in_values(E)) - sum(self.out_values(E))
        p[3] = numpy.sqrt(numpy.abs(E[3]**2 - m[3]**2))
        return p, E, m

    def integrand(self, p=[], order=(0, 1, 2, 3)):
        """ Total collision integral interior with performance optimizations """
        p, E, m = self.calculate_impulses(p)
        integrand = theta(E[3] - m[3])

        if integrand == 0:
            return 0

        if p[0] != 0:
            integrand *= D(p=p, E=E, m=m, K1=self.K1, K2=self.K2, order=order) / p[0] / E[0]
        else:
            integrand *= Db1(*p[1:]) * m[1] * (E[2] * E[3] + Db2(*p[1:]))

        if integrand == 0:
            return 0

        for i in [1, 2, 3]:
            if m[i] != 0:
                integrand *= p[i] / E[i]

        integrand *= (self.F_B(p) - self.F_A(p))

        return integrand

    def F_f(self, p=[]):
        """ Variable part of the collision integral ready for integration """
        return -1. * (
            self.F_A(p=p, skip_index=0) + self.in_particles[0].eta * self.F_B(p=p, skip_index=0)
        )

    def F_1(self, p=[]):
        """ Constant part of the collision integral ready for integration """
        return self.F_B(p=p, skip_index=0)

    """
    \begin{equation}
        \mathcal{F}(f) = f_3 f_4 (1 \pm f_1) (1 \pm f_2) - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
    \end{equation}

    \begin{equation}
        \mathcal{F}(f) = f_1 (\mp f_3 f_4 (1 \pm f_2) - f_2 (1 \pm f_3) (1 \pm f_4)) \
        + f_3 f_4 (1 \pm f_2)
    \end{equation}

    \begin{equation}
        \mathcal{F}(f) = \mathcal{F}_B^{(1)} - f_1 (\mathcal{F}_A^{(1)} \mp \mathcal{F}_B^{(1)})
    \end{equation}
    """
    def F_A(self, p=[], skip_index=None):
        """
        Forward reaction distribution functional term

        \begin{equation}
            \mathcal{F}_A = f_1 f_2 (1 \pm f_3) (1 \pm f_4)
        \end{equation}

        :param skip_index: Particle to skip in the expression
        """
        in_p = self.in_values(p)
        out_p = self.out_values(p)
        temp = 1.
        for i, particle in enumerate(self.in_particles):
            if i != skip_index:
                temp *= particle.distribution(in_p[i])
        for i, particle in enumerate(self.out_particles):
            temp *= 1. - particle.eta * particle.distribution(out_p[i])
        return temp

    def F_B(self, p=[], skip_index=None):
        """
        Backward reaction distribution functional term

        \begin{equation}
            \mathcal{F}_B = f_3 f_4 (1 \pm f_1) (1 \pm f_2)
        \end{equation}

        :param skip_index: Particle to skip in the expression
        """
        in_p = self.in_values(p)
        out_p = self.out_values(p)
        temp = 1.
        for i, particle in enumerate(self.out_particles):
            temp *= particle.distribution(out_p[i])
        for i, particle in enumerate(self.in_particles):
            if i != skip_index:
                temp *= 1. - particle.eta * particle.distribution(in_p[i])
        return temp

    def in_values(self, p=[]):
        return p[:len(self.in_particles)]

    def out_values(self, p=[]):
        return p[len(self.in_particles):]
