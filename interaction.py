# -*- coding: utf-8 -*-
import copy
import numpy
from common import PARAMS, theta, GRID, CONST
from common.utils import benchmark
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


class M(object):

    """ == Matrix element ==
        All four-particle interactions of the interest can be rewritten in a form

        \begin{equation}
            |\mathcal{M}|^2 = \sum_{\{i \neq j \neq k \neq l\}} K_1 (p_i \cdot p_j) (p_k \cdot p_l)\
                 + K_2 m_i m_j (p_k \cdot p_l)
        \end{equation} """
    K1 = 0.
    K2 = 0.
    # Order defines the values of the $(i, j, k, l)$ indices
    order = (0, 1, 2, 3)

    def __init__(self, *args, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def __str__(self):
        return "K1={: .2e}, K2={: .2e}, {}".format(self.K1, self.K2, self.order)


class WeakM(M):
    def __init__(self, *args, **kwargs):
        super(WeakM, self).__init__(*args, **kwargs)

        self.K1 *= 32 * CONST.G_F**2
        self.K2 *= 32 * CONST.G_F**2


class Interaction:

    """ == Interaction ==
        Helper class that takes care of creating all necessary `Integral`s for the actual\
        interaction.
    """

    name = "Particle interaction"

    integrals = []

    in_particles = []  # Incoming particles
    out_particles = []  # Outgoing particles
    particles = []  # All particles involved
    # Temperature when the typical interaction time exceeds the Hubble expansion time
    decoupling_temperature = 0.
    # Interaction symmetry factor
    symmetry_factor = 1.

    # Matrix elements of the interaction
    Ms = []

    """ 2-to-2 interactions and 1-to-3 decays can be generalized by introducing a multiplier `s` \
        for the momentum of the second particle: $p_1 + s p_2 = p_3 + p_4 $ """
    s = 1.

    def __init__(self, *args, **kwargs):
        """ Init """

        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.integrals = []
        self.particles = self.in_particles + self.out_particles

        Ms = copy.copy(self.Ms)

        # Remember all particle species we've already considered to avoid double-counting
        accounted_particles = set()

        # Step 1. Permute all in_particles
        in_particles = self.in_particles
        out_particles = self.out_particles
        accounted_particles = self.init_integrals(in_particles, out_particles,
                                                  Ms, accounted_particles)

        # Step 2. Turn to the backward process
        for M in Ms:
            M.order = M.order[len(self.in_particles):] + M.order[:len(self.in_particles)]

        # Step 3. Permute all new in_particles (former out_particle)
        accounted_particles = self.init_integrals(out_particles, in_particles,
                                                  Ms, accounted_particles)

        print self

    def __str__(self):
        return self.name + "\n\t" + "\n\t".join([str(integral) for integral in self.integrals])

    def init_integrals(self, in_particles, out_particles, Ms, accounted_particles):
        for i, particle in enumerate(in_particles):

            if particle in accounted_particles:
                # Skip already accounted species
                continue

            particle_Ms = copy.copy(Ms)
            for M in particle_Ms:
                M.order = M.order[i:i+1] + M.order[:i] + M.order[i+1:]

            # Add interaction integrals by putting each incoming particle as the first one
            self.integrals.append(Integral(
                in_particles=in_particles[i:i+1] + in_particles[:i] + in_particles[i+1:],
                out_particles=out_particles,
                decoupling_temperature=self.decoupling_temperature,
                Ms=particle_Ms
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

    Ms = []
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
        return " + ".join([p.symbol for p in self.in_particles]) \
            + " ⟶  " + " + ".join([p.symbol for p in self.out_particles]) \
            + "\t({})".format(len(self.Ms))

    def __repr__(self):
        return self.__str__()

    def calculate(self):
        """
        Calculate collision integral constants and save them to the first involved particle
        """
        if PARAMS.T > self.decoupling_temperature and not self.in_particles[0].in_equilibrium:

            self.particles[0].collision_integrands.append(self)

    def calculate_impulses(self, p=[]):
        """ Helper procedure that caches energies and normalized masses of particles """
        p = (p + [0., 0., 0., 0.])[:4]
        E = []
        m = []
        for i, particle in enumerate(self.particles):
            E.append(particle.energy_normalized(p[i]))
            m.append(particle.mass_normalized)
        E[3] = sum(self.in_values(E)) - sum(self.out_values(E))
        p[3] = numpy.sqrt(numpy.abs(E[3]**2 - m[3]**2))
        return p, E, m

    def integrand(self, p=[], F_A=True, F_B=True, F_1=False, F_f=False):
        """ Total collision integral interior with performance optimizations """

        p, E, m = self.calculate_impulses(p)
        if not self.in_bounds(p, E, m):
            return 0.

        integrand = 1./64. / numpy.pi**3 * PARAMS.m**5 / PARAMS.x**6 / PARAMS.H

        ds = 0.
        for M in self.Ms:
            if p[0] != 0:
                ds += D(p=p, E=E, m=m, K1=M.K1, K2=M.K2, order=M.order) / p[0] / E[0]
            else:
                ds += Db1(*p[1:]) + m[1] * (E[2] * E[3] + Db2(*p[1:]))

        integrand *= ds

        for i in [1, 2, 3]:
            if m[i] != 0:
                integrand *= p[i] / E[i]

        if integrand == 0:
            return 0

        fau = 0
        if F_B:
            fau += self.F_B(p)
        if F_A:
            fau -= self.F_A(p)
        if F_1:
            fau += self.F_1(p)
        if F_f:
            fau += self.F_f(p)

        integrand *= fau

        return integrand

    def in_bounds(self, p=[], E=None, m=None):
        if not E or not m:
            p, E, m = self.calculate_impulses(p)

        q1, q2 = (p[0], p[1]) if p[0] > p[1] else (p[1], p[0])
        q3, q4 = (p[2], p[3]) if p[2] > p[3] else (p[3], p[2])

        is_in = E[3] >= m[3] \
            and q1 <= q2 + q3 + q4 \
            and q3 <= q1 + q2 + q4

        return is_in

    def bounds(self, p0):
        points = []
        for p1 in GRID.TEMPLATE:
            points.append((p1, self.lower_bound(p0, p1),))
            points.append((p1, self.upper_bound(p0, p1),))

        return points

    def lower_bound(self, p0, p1):

        index = 0
        while index < GRID.MOMENTUM_SAMPLES and not self.in_bounds([p0, p1, GRID.TEMPLATE[index], 0]):
            index += 1

        if index == GRID.MOMENTUM_SAMPLES:
            return GRID.MIN_MOMENTUM

        return GRID.TEMPLATE[index]

    def upper_bound(self, p0, p1):

        index = int((min(p0 + p1, GRID.MAX_MOMENTUM) - GRID.MIN_MOMENTUM) / GRID.MOMENTUM_STEP)

        while index >= 0 and not self.in_bounds([p0, p1, GRID.TEMPLATE[index], 0]):
            index -= 1

        if index == -1:
            return GRID.MIN_MOMENTUM

        return GRID.TEMPLATE[index]

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
            if skip_index is None or i != skip_index:
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
            if skip_index is None or i != skip_index:
                temp *= 1. - particle.eta * particle.distribution(in_p[i])
        return temp

    def in_values(self, p=[]):
        return p[:len(self.in_particles)]

    def out_values(self, p=[]):
        return p[len(self.in_particles):]
