# -*- coding: utf-8 -*-

"""
# Common

This file contains constants and utilities shared by all other modules in the project.
"""
import numpy
import numericalunits as nu
from utils import PicklableObject


class UNITS(object):

    __slots__ = ['use_numericalunits', 'eV', 'keV', 'MeV', 'GeV', 'TeV', 's', 'kg', 'm', 'N']

    """ ## Units
        As we use natural units in the project, all units from `numericalunits` except energy units\
        are useless. Here some useful units are defined in terms of `GeV`s. """

    use_numericalunits = False
    eV = nu.eV if use_numericalunits else 1e-9

    @classmethod
    def reset_units(cls):
        UNITS.keV = UNITS.eV * 1e3
        UNITS.MeV = UNITS.keV * 1e3
        UNITS.GeV = UNITS.MeV * 1e3
        UNITS.TeV = UNITS.GeV * 1e3
        UNITS.s = 1. / 6.58 * 1e25 / UNITS.GeV
        UNITS.kg = 1e27 / 1.8 * UNITS.GeV
        UNITS.m = 1e15 / 0.197 / UNITS.GeV
        UNITS.N = 1e-5 / 8.19 * UNITS.GeV**2

UNITS.reset_units()


class CONST(object):

    __slots__ = ['G', 'M_p', 'G_F', 'sin_theta_w_2', 'g_R', 'g_L', 'f_pi',
                 'MeV_to_s_1', 'MeV_to_10_9K', 'MeV4_to_g_cm_3', 'rate_normalization']

    """ ### Physical constants """

    # Gravitational constant
    G = 6.67 * 1e-11 * (UNITS.N / UNITS.kg**2 * UNITS.m**2)
    # Planck mass
    M_p = 1.2209 * 1e19 * UNITS.GeV
    # Fermi constant
    G_F = 1.166 * 1e-5 / UNITS.GeV**2
    # Weinberg angle
    sin_theta_w_2 = 0.2312
    g_R = sin_theta_w_2
    g_L = sin_theta_w_2 + 0.5
    # Pion decay constant
    f_pi = 130. * UNITS.MeV

    MeV_to_s_1 = 1.51926758e21
    MeV_to_10_9K = 11.6045
    MeV4_to_g_cm_3 = 2.32011575e5

    rate_normalization = 17.54459


class Params(PicklableObject):

    __slots__ = ['T_initial', 'T_final', 'm', 'dy', 't', 'H', 'rho',
                 'a_initial', 'a', 'x', 'y', 'dx', 'T', 'aT', '_saveable_fields']

    _saveable_fields = __slots__

    """ ## Parameters
        Master object carrying the cosmological state of the system and initial conditions """

    # Temperature bounds define the simulations boundaries of the system
    T_initial = 10 * UNITS.MeV
    T_final = 10 * UNITS.keV

    # Arbitrary normalization of the conformal scale factor
    m = 1. * UNITS.MeV
    # Conformal scale factor step size during computations
    dy = 0.05
    # Initial time
    t = 0. * UNITS.s
    # Hubble rate
    H = 0.
    # Total energy density
    rho = 0.

    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])

        self.infer()

    def infer(self):
        """ Set initial cosmological parameters based on the value of `T_initial` """
        # As the initial scale factor is arbitrary, it can be use to ensure the initial $aT$ value\
        # equal to 1
        self.a_initial = self.m / self.T_initial

        # Compute present-state parameters that can be inferred from the base ones
        self.a = self.a_initial
        self.x = self.a * self.m
        self.y = numpy.log(self.x)
        self.dx = self.x * (numpy.exp(self.dy) - 1.)
        self.T = self.T_initial
        self.aT = self.a * self.T

    def update(self, rho):
        """ Hubble expansion parameter defined by a Friedmann equation:

            \begin{equation}
                H = \sqrt{\frac{8 \pi}{3} G \rho}
            \end{equation}
        """
        self.rho = rho
        self.H = numpy.sqrt(8./3.*numpy.pi * rho) / CONST.M_p

        old_a = self.a
        """ Physical scale factor and temperature for convenience """
        self.a = self.x / self.m
        self.T = self.aT / self.a

        """ Time step size is inferred from the approximation of the scale factor `a` \
            derivative and a definition of the Hubble parameter `H`:

            \begin{equation}
                H = \frac{\dot{a}}{a} = \frac{1}{a_{i-1}} \frac{a_i - a_{i-1}}{\Delta t} \
                  = \frac{1}{\Delta t}(\frac{a_i}{a_{i-1}} -1)
            \end{equation}

            \begin{equation}
                \Delta t = \frac{1}{H} (\frac{a_i}{a_{i-1}} - 1)
            \end{equation}
        """
        dt = (self.a / old_a - 1) / self.H
        self.t += dt


class Grid(PicklableObject):

    """ ### Distribution functions grid

        TODO: try a log-spaced grid instead of equally-spaced
        TODO: try an irregular grid based on the Gauss-Legendre quadrature roots to minimize \
              the interpolation for the massless particles.

        To capture non-equilibrium effects in the Early Universe, we work with particle species \
        distribution functions $f(\vec{p}, \vec{r}, t)$. Assuming that the Universe is homogeneous\
        and isotropic, we can forget dependency on the position and the momentum direction: \
        $f(p, t)$.

        Resulting functions are sampled across a wide range of momenta. However, momentum range\
        cannot include both 0 momenta and very large momenta (either leads to numerical overflows\
        and errors).
        """

    __slots__ = ['MIN_MOMENTUM', 'MAX_MOMENTUM', 'MOMENTUM_SAMPLES', 'MOMENTUM_STEP', 'TEMPLATE']

    def __init__(self):
        self.MIN_MOMENTUM = 1. * UNITS.eV
        self.MAX_MOMENTUM = self.MIN_MOMENTUM + 20 * UNITS.MeV
        self.MOMENTUM_SAMPLES = 50

        """
        Grid template can be copied when defining a new distribution function and is convenient to\
        calculate any _vectorized_ function over the grid. For example,

        ```python
        numpy.vectorize(particle.conformal_energy)(GRID.TEMPLATE)
        ```

        yields an array of particle conformal energy mapped over the `GRID`
        """
        self.TEMPLATE = numpy.linspace(self.MIN_MOMENTUM, self.MAX_MOMENTUM,
                                       num=self.MOMENTUM_SAMPLES, endpoint=True)

        self.MOMENTUM_STEP = self.TEMPLATE[1] - self.TEMPLATE[0]


GRID = Grid()


def theta(f):
    """ Heaviside $\theta$ function """
    if f > 0:
        return 1.
    if f == 0:
        return 0.5
    return 0.


def linear_interpolation(function, grid):
    """ ### Linear interpolation """

    def interpolation(x):
        index = numpy.searchsorted(grid, x)

        if index >= len(grid) - 1:
            return function[len(grid) - 1]

        if x == grid[index]:
            return function[index]

        # Determine the closest grid points
        i_low = index
        x_low = grid[i_low]

        i_high = index + 1
        x_high = grid[i_high]

        return (
            function[i_low] * (x_high - x) + function[i_high] * (x - x_low)
        ) / (x_high - x_low)

    return interpolation
