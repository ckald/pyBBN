# -*- coding: utf-8 -*-
import copy
import numpy
from common import GRID, integrators
from common.utils import PicklableObject
from interactions.four_particle import FourParticleIntegral


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
            self.integrals.append(FourParticleIntegral(
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
