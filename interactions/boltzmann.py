# -*- coding: utf-8 -*-
import numpy
from common import integrators


class BoltzmannIntegral(object):

    """ ## Integral
        Representation of the concrete collision integral for a specific particle \
        `Integral.reaction[0]` """

    reaction = None  # All particles involved
    sides = None

    """ ### Crossed particles in the integral

        If the Boltzmann integral is obtained as a crossing of some other process, the mass of the\
        crossed fermion in the matrix element has to change sign.

        For example, compare the matrix elements of reactions

        \begin{align}
            \nu + e &\to \nu + e \\\\
            \nu + \overline{\nu} &\to e + e^+
        \end{align}

        In particular, this has something to do with the `K2` term in the `FourParticleIntegral`.
    """

    # Temperature when the typical interaction time exceeds the Hubble expansion time
    washout_temperature = 0.
    constant = 0.

    """ Four-particle interactions of the interest can all be rewritten in a form

        \begin{equation}
            |\mathcal{M}|^2 = \sum_{\{i \neq j \neq k \neq l\}} K_1 (p_i \cdot p_j) (p_k \cdot p_l)\
                 + K_2 m_i m_j (p_k \cdot p_l)
        \end{equation} """

    Ms = None

    """ Grids corresponding to particles integrated over """
    grids = None

    def __init__(self, **kwargs):
        """ Update self with configuration `kwargs`, construct particles list and \
            energy conservation law of the integral. """

        for key in kwargs:
            setattr(self, key, kwargs[key])

        if not self.Ms:
            self.Ms = []

        self.sides = tuple([item.side for item in self.reaction])

    def __str__(self):
        """ String-like representation of the integral. Corresponds to the first particle """
        return (
            " + ".join([p.specie.symbol + ('\'' if p.antiparticle else '')
                        for p in self.reaction if p.side == -1])
            + " ⟶  "
            + " + ".join([p.specie.symbol + ('\'' if p.antiparticle else '')
                          for p in self.reaction if p.side == 1])
            + "\t({})".format(', '.join([str(M) for M in self.Ms]))
        )

    def __repr__(self):
        return self.__str__()

    def initialize(self):
        """
        Initialize collision integral constants and save them to the first involved particle
        """
        raise NotImplementedError()

    def rates(self):
        forward_integral = numpy.vectorize(lambda p0: p0**2 / (2 * numpy.pi)**3
                                           * self.integrate(p0, self.F_A)[0])

        backward_integral = numpy.vectorize(lambda p0: p0**2 / (2 * numpy.pi)**3
                                            * self.integrate(p0, self.F_B)[0])

        grid = self.particle.grid

        forward_rate, _ = integrators.integrate_1D(forward_integral,
                                                   (grid.MIN_MOMENTUM, grid.MAX_MOMENTUM))

        backward_rate, _ = integrators.integrate_1D(backward_integral,
                                                    (grid.MIN_MOMENTUM, grid.MAX_MOMENTUM))

        return -forward_rate, backward_rate

    def rate(self):
        forward_rate, backward_rate = self.rates()
        return backward_rate - forward_rate

    @staticmethod
    def integrate(p0, integrand, bounds=None, kwargs=None):
        raise NotImplementedError()

    def integrand(self, *args, **kwargs):
        """ Collision integral interior. """
        raise NotImplementedError()
