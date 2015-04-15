# -*- coding: utf-8 -*-
import copy
import itertools
from collections import namedtuple, Counter
from common.utils import PicklableObject
from interactions.four_particle import FourParticleM


"""
## Boltzmann collision integral

\begin{multline}
    I_{coll}(t,p_i) = \frac{1}{2 E_i} \sum_{reactions} \int \cdots \int
        \prod_{j} \frac{d^3 p_{in,j}}{(2 \pi)^3 2 E_{in,j}}
        \prod_{k} \frac{d^3 p_{out,k}}{(2 \pi)^3 2 E_{out,k}} \\
        S |\mathcal{M}|^2 \mathcal{F}(\{f_\alpha\}) (2 \pi)^4 \\
        \delta^4(\sum_{m} p_{in,m}^\mu - \sum_{n} p_{out,n}^\mu)
\end{multline}
"""

IntegralItem = namedtuple('IntegralItem', ['side', 'specie', 'antiparticle', 'index', 'crossed'])


class Interaction(PicklableObject):

    """
    ## Interaction
    Helper class that takes care of creating all necessary `Integral`s for the actual interaction.
    """

    # Human-friendly interaction identifier
    name = "Particle interaction"

    integrals = None
    particles = None  # All particles involved

    # Temperature when the typical interaction time exceeds the Hubble expansion time
    decoupling_temperature = 0.

    # Matrix elements of the interaction
    Ms = None

    def __init__(self, **kwargs):
        """ Create an `Integral` object for all particle species involved in the interaction.

            Precise expressions for all integrals can be derived by permuting all particle-related\
            functions in the distribution functions, matrix elements, momenta, energy and mass\
            arrays.

            To avoid double-counting, one should create an integral for each particle specie only \
            once.

            Additionally, one has to generate the integrals for all the crossing processes by\

        """
        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.integrals = []

        for reaction in self.reactions_map():
            lhs = [item for item in reaction if item.side == -1]
            item = lhs[0]

            # Skip collision integrals for anti-particles
            if item.antiparticle:
                continue

            self.init_integrals(reaction, item)

        self.integrals = IntegralSet(self.integrals).integrals

    def __str__(self):
        """ String-like representation of the interaction and all its integral """
        return self.name + "\n\t" + "\n\t".join([str(integral) for integral in self.integrals])

    @classmethod
    def cross_particle(cls, item):
        return item._replace(
            crossed=not item.crossed,
            antiparticle=not item.antiparticle if not item.specie.majorana else False,
            side=-item.side
        )

    def permutations(self, particles, antiparticles):
        return itertools.permutations(
            [
                IntegralItem(side=-1,
                             specie=particle,
                             antiparticle=antiparticle,
                             index=i,
                             crossed=False)
                for i, (particle, antiparticle)
                in enumerate(zip(particles[0], antiparticles[0]))
            ] + [
                IntegralItem(side=1,
                             specie=particle,
                             antiparticle=antiparticle,
                             index=i+len(particles[0]),
                             crossed=False)
                for i, (particle, antiparticle)
                in enumerate(zip(particles[1], antiparticles[1]))
            ]
        )

    def reactions_map(self):
        """ Return all relevant collision integrals involving `count` particles """
        reactions = {}

        perms = self.permutations(self.particles, self.antiparticles)

        # if any(specie.majorana for specie in itertools.chain.from_iterable(self.particles)):
        #     print "There is Majorana in da house!"
        #     conjugated_antiparticles = (
        #         tuple(not antip if not specie.majorana else False
        #               for specie, antip in zip(self.particles[0], self.antiparticles[0])),
        #         tuple(not antip if not specie.majorana else False
        #               for specie, antip in zip(self.particles[1], self.antiparticles[1]))
        #     )
        #     perms = itertools.chain(perms,
        #                             self.permutations(self.particles, conjugated_antiparticles))

        # Generate all possible permutations of particles
        for perm in perms:
            # Generate possible reactions by inserting an arrow at every position
            for i in range(sum(map(len, self.particles)) - 1):

                # Cross particles according to their position relative to the reaction arrow
                reaction = tuple(
                    Interaction.cross_particle(item)
                    if (item.side > 0 and j <= i) or (item.side < 0 and j > i)
                    else item
                    for j, item in enumerate(perm)
                )

                lhs = [item for item in reaction if item.side == -1]
                rhs = [item for item in reaction if item.side == 1]

                lhs = lhs[0:1] + sorted(lhs[1:])
                rhs = sorted(rhs)

                if lhs[0].antiparticle:
                    continue

                # Skip kinematically impossible reactions
                if (
                    (len(lhs) == 1 and lhs[0].specie.mass <= sum(item.specie.mass for item in rhs))
                    or
                    (len(rhs) == 1 and rhs[0].specie.mass <= sum(item.specie.mass for item in lhs))
                ):
                    continue

                k = tuple(map(tuple, lhs + rhs))
                if k not in reactions:
                    reactions[k] = tuple(lhs + rhs)

        return reactions.values()

    def init_integrals(self, reaction, item):
        """ Starting from the single representation of the interaction, create integrals objects\
            to be computed for every particle specie involved """

        particle = item.specie
        particle_Ms = copy.deepcopy(self.Ms)

        for M in particle_Ms:
            if M.__class__ is FourParticleM:
                map = {item.index: i for i, item in enumerate(reaction)}
                M.apply_order(tuple(map[val] for val in M.order), reaction)

        # Add interaction integrals by putting each incoming particle as the first one
        self.integrals.append(self.integral(
            particle=particle,
            reaction=reaction,
            decoupling_temperature=self.decoupling_temperature,
            Ms=particle_Ms,
        ))

    def initialize(self):
        """ Proxy method """

        for integral in self.integrals:
            integral.initialize()


class IntegralSet(PicklableObject):

    integrals = None

    def __init__(self, new_integrals):
        self.integrals = []
        for new_integral in new_integrals:
            self.append(new_integral)

    def species(self, integral):
        return (
            Counter(item.specie for item in integral.reaction if item.side == -1),
            Counter(item.specie for item in integral.reaction if item.side == 1)
        )

    def append(self, new_integral):
        old_integral = None

        new_species = self.species(new_integral)

        # Check if there is already a similar integral
        for integral in self.integrals:
            if new_integral.particle == integral.particle and new_species == self.species(integral):
                old_integral = integral
                break

        # If there is no, just append a new one
        if not old_integral:
            self.integrals.append(new_integral)
            return

        # Otherwise, check all matrix elements
        Ms = []
        for new_M in new_integral.Ms:
            reduced = False
            for old_M in old_integral.Ms:
                if old_M.stackable(new_M):
                    reduced = True
                    old_M += new_M
                    break

            if not reduced:
                Ms.append(new_M)

        if Ms:
            old_integral.Ms = tuple(list(old_integral.Ms) + Ms)
