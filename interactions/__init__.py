# -*- coding: utf-8 -*-
import copy
from common import CONST
from common.utils import PicklableObject
from interactions.scattering import Integral


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

    def __init__(self, **kwargs):
        """ Configure matrix element, check that it makes sense """
        for key in kwargs:
            setattr(self, key, kwargs[key])

        if len(set(self.order)) != len(self.order):
            raise Exception("Meaningless order of momenta of the matrix element: {}"
                            .format(self.order))

    def __str__(self):
        """ String-like representation of the matrix element """
        return "K1={: .2e}, K2={: .2e}, {}".format(self.K1, self.K2, self.order)


class WeakM(M):

    """ == Weak interactions matrix element ==
        Weak processes usually include a common factor of $32 G_F^2$ """

    def __init__(self, *args, **kwargs):
        super(WeakM, self).__init__(*args, **kwargs)

        self.K1 *= 32 * CONST.G_F**2
        self.K2 *= 32 * CONST.G_F**2


class Interaction(PicklableObject):

    """
    == Interaction ==
    Helper class that takes care of creating all necessary `Integral`s for the actual interaction.
    """

    # Human-friendly interaction identifier
    name = "Particle interaction"

    integrals = []

    in_particles = []  # Incoming particles
    out_particles = []  # Outgoing particles
    particles = []  # All particles involved

    # Temperature when the typical interaction time exceeds the Hubble expansion time
    decoupling_temperature = 0.

    # Matrix elements of the interaction
    Ms = []

    def __init__(self, **kwargs):
        """ Create an `Integral` object for all particle species involved in the interaction.

            Precise expressions for all integrals can be derived by permuting all particle-related\
            functions in the distribution functions, matrix elements, momenta, energy and mass\
            arrays.

            To avoid double-counting, one should create an integral for each particle specie only \
            once.
        """

        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.integrals = []
        self.particles = self.in_particles + self.out_particles

        Ms = copy.deepcopy(self.Ms)

        # Remember all particle species we've already considered to avoid double-counting
        accounted_particles = set()

        # === Interaction integrals initialization strategy ===

        # ==== 1. Permute all in_particles ====
        in_particles = self.in_particles
        out_particles = self.out_particles
        accounted_particles = self.init_integrals(in_particles, out_particles,
                                                  Ms, accounted_particles)

        # ==== 2. Turn to the backward process ====
        # `(0, 1, 2, 3) -> (2, 3, 0, 1)`
        for M in Ms:
            M.order = M.order[len(self.in_particles):] + M.order[:len(self.in_particles)]

        # ==== 3. Permute all new `in_particles` (former `out_particles`) ====
        accounted_particles = self.init_integrals(out_particles, in_particles,
                                                  Ms, accounted_particles)

    def __str__(self):
        """ String-like representation of the interaction and all its integral """
        return self.name + "\n\t" + "\n\t".join([str(integral) for integral in self.integrals])

    def init_integrals(self, in_particles, out_particles, Ms, accounted_particles):
        """ Starting from the single representation of the interaction, create integrals objects\
            to be computed for every particle specie involved """
        for i, particle in enumerate(in_particles):

            # Skip already accounted species
            if particle in accounted_particles:
                continue

            particle_Ms = copy.deepcopy(Ms)

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

    def initialize(self):
        """ Proxy method """

        for integral in self.integrals:
            integral.initialize()
