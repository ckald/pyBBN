# -*- coding: utf-8 -*-

"""
= Common =

This file contains constants and utilities shared by all other modules in the project.
"""
import numpy
import numericalunits as nu


# Numpy error handling behavior. Uncomment to see a notice each time you get an overflow.
# numpy.seterr(all='print')


class UNITS:

    """ == Units ==
        As we use natural units in the project, all units from `numericalunits` except energy units\
        are useless. Here some useful units are defined in terms of `GeV`s. """

    # eV = nu.eV
    eV = 1e-9

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


class CONST:
    """ === Physical constants === """

    G = 6.67 * 1e-11 * (UNITS.N / UNITS.kg**2 * UNITS.m**2)
    G_F = 1.166 * 1e-5 / UNITS.GeV**2
    sin_theta_w_2 = 0.23
    g_R = sin_theta_w_2
    g_L = sin_theta_w_2 + 0.5


class PARAMS:

    """ == Parameters ==
        Master object carrying the cosmological state of the system and initial conditions """

    # Temperature bounds define the simulations boundaries of the system
    T_initial = 10 * UNITS.MeV
    T_final = 10 * UNITS.keV

    # Arbitrary normalization of the conformal scale factor
    m = 1. * UNITS.MeV
    # Conformal scale factor step size during computations
    dx = 1e-2 * UNITS.MeV
    # Initial time
    t = 0. * UNITS.s
    # Hubble rate
    H = 0.
    # Total energy density
    rho = 0.

    @classmethod
    def infer(cls):
        # Initial scale factor - arbitrary
        PARAMS.a_initial = 1. / (PARAMS.T_initial / UNITS.MeV)
        # PARAMS.a_initial = 1.

        # Compute present-state parameters that can be inferred from the base ones
        PARAMS.a = PARAMS.a_initial
        PARAMS.x = PARAMS.a * PARAMS.m
        PARAMS.T = PARAMS.T_initial
        PARAMS.aT = PARAMS.a * PARAMS.T

        GRID.init()

    @classmethod
    def update(cls, rho):
        """ Hubble expansion parameter defined by a Friedmann equation:

            \begin{equation}
                H = \sqrt{\frac{8 \pi}{3} G \rho}
            \end{equation}
        """
        PARAMS.rho = rho
        PARAMS.H = numpy.sqrt(8./3.*numpy.pi * CONST.G * rho)

        old_a = PARAMS.a
        """ Physical scale factor and temperature for convenience """
        PARAMS.a = PARAMS.x / PARAMS.m
        PARAMS.T = PARAMS.aT / PARAMS.a

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
        dt = (PARAMS.a / old_a - 1) / PARAMS.H
        PARAMS.t += dt


class GRID:

    """ === Distribution functions grid ===

        TODO: try a log-spaced grid instead of equally-spaced

        To capture non-equilibrium effects in the Early Universe, we work with particle species \
        distribution functions $f(\vec{p}, \vec{r}, t)$. Assuming that the Universe is homogeneous\
        and isotropic, we can forget dependency on the position and the momentum direction: \
        $f(p, t)$.

        Resulting functions are sampled across a wide range of momenta. However, momentum range\
        cannot include both 0 momenta and very large momenta (either leads to numerical overflows\
        and errors).
        """

    @classmethod
    def init(cls):
        GRID.MIN_MOMENTUM = 1. * UNITS.eV
        GRID.MAX_MOMENTUM = GRID.MIN_MOMENTUM + 20 * UNITS.MeV
        GRID.MOMENTUM_SAMPLES = 50

        # Grid template can be copied when defining a new distribution function
        GRID.TEMPLATE = numpy.linspace(GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
                                       num=GRID.MOMENTUM_SAMPLES, endpoint=True)

        GRID.MOMENTUM_STEP = GRID.TEMPLATE[1] - GRID.TEMPLATE[0]

PARAMS.infer()


def theta(f):
    """ Theta $\theta$ function """
    if f > 0:
        return 1.
    if f == 0:
        return 0.5
    return 0.
