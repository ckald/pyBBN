# -*- coding: utf-8 -*-

"""
# Common

This file contains constants and utilities shared by all other modules in the project.
"""
import numpy
import numericalunits as nu

import environment
from common import statistics as STATISTICS

class UNITS(object):

    """ ## Units
        As we use natural units in the project, all units from `numericalunits` except energy units
        are useless. Here some useful units are defined in terms of `GeV`s. """

    use_numericalunits = False
    eV = nu.eV if use_numericalunits else 1e-9

    @classmethod
    def reset_units(cls):
        UNITS.keV = UNITS.eV * 1e3
        UNITS.MeV = UNITS.keV * 1e3
        UNITS.GeV = UNITS.MeV * 1e3
        UNITS.TeV = UNITS.GeV * 1e3
        UNITS.s = 1e22 / 6.582119 / UNITS.MeV
        UNITS.kg = 1e27 / 1.782662 * UNITS.GeV
        UNITS.m = 1e15 / 1.239842 / UNITS.GeV
        UNITS.N = 1e-5 / 8.19 * UNITS.GeV**2

        # Temperature: $10^9 K$
        UNITS.K9 = UNITS.MeV / 11.6045
        UNITS.g_cm3 = UNITS.MeV**4 / 2.32011575e5

UNITS.reset_units()


class CONST(object):

    """ ### Physical constants """

    # Planck mass
    M_p = 1.2209 * 1e22 * UNITS.MeV
    # Gravitational constant
    # G = 6.67 * 1e-11 * (UNITS.N / UNITS.kg**2 * UNITS.m**2)
    G = 1 / M_p**2
    # Fermi constant
    G_F = 1.166 * 1e-5 / UNITS.GeV**2
    # Hubble constant
    H = 1. / (4.55e17 * UNITS.s)
    # Weinberg angle
    sin_theta_w_2 = 0.2312
    g_R = sin_theta_w_2
    g_L = sin_theta_w_2 + 0.5
    # Decay constants
    f_K = 155.6 * UNITS.MeV # kaon
    f_eta_c = 335. * UNITS.MeV
    f_D = 212. * UNITS.MeV
    f_Ds = 249. * UNITS.MeV
    f_Ds_star = 315. * UNITS.MeV
    # Baryon-to-photon ratio
    eta = 6.1e-10
    # Fine-structure constant
    alpha = 1 / 137.036
    # QCD transition temperature
    lambda_QCD = 150 * UNITS.MeV
    # W_boson mass
    M_W = 80.379 * UNITS.GeV
    # Neutrino squared mass differences
    delta_m12_sq = 7.37e-5 * UNITS.eV**2
    delta_m13_sq = 2.56e-3 * UNITS.eV**2

    rate_normalization = 17.5426 / UNITS.MeV**5 #17.54459

    MSW_constant = 16. * 1.2021 * numpy.sqrt(2) * G_F / numpy.pi / M_W**2


class Params(object):

    a = None
    T = None
    m = None
    t = None
    H = None
    rho = None
    N_eff = None

    x = None
    dx = None
    y = None
    dy = None

    h = None

    def __init__(self, **kwargs):
        """ ## Parameters
            Master object carrying the cosmological state of the system and initial conditions """

        # Temperature bounds define the simulations boundaries of the system
        self.T = None

        # Arbitrary normalization of the conformal scale factor
        self.m = 1. * UNITS.MeV
        # Initial time
        self.t = 0. * UNITS.s
        # Hubble rate
        self.H = 0.
        # Total energy density
        self.rho = None
        self.N_eff = 0.

        self.dy = 0.05
        self.dx = 0.005 * UNITS.MeV

        for key in kwargs:
            setattr(self, key, kwargs[key])

        # As the initial scale factor is arbitrary, it can be use to ensure the initial $aT$ value\
        # equal to 1
        self.a = 1 * self.m / self.T
        self.a_ini = self.a

        self.infer()

        if environment.get('LOGARITHMIC_TIMESTEP') and not kwargs.get('dy'):
            raise Exception("Using logarithmic timestep, but no Params.dy was specified")
        if not environment.get('LOGARITHMIC_TIMESTEP') and not kwargs.get('dx'):
            raise Exception("Using linear timestep, but no Params.dx was specified")

    def infer(self):
        """ Set initial cosmological parameters based on the value of `T` """

        # Compute present-state parameters that can be inferred from the base ones
        self.x = self.a * self.m
        if environment.get('LOGARITHMIC_TIMESTEP'):
            self.y = numpy.log(self.a)
        self.aT = self.a * self.T

        # Conformal scale factor step size during computations
        if environment.get('LOGARITHMIC_TIMESTEP'):
            self.dx = self.x * self.dy
            self.h = self.dy
        else:
            self.dy = None
            self.h = self.dx

    def init_time(self, rho):
        self.t = numpy.sqrt(3. / (32. * numpy.pi * CONST.G * rho))
        return self.t

    @staticmethod
    def entropic_dof_eq(universe):
        dof = 0
        for particle in universe.particles:
            if universe.params.T > particle.mass and particle.in_equilibrium:
                if particle.statistics == STATISTICS.FERMION:
                    dof += particle.dof * 7. / 8.
                else:
                    dof += particle.dof
        return dof

    def update(self, rho, S):
        """ Hubble expansion parameter defined by a Friedmann equation:

            \begin{equation}
                H = \sqrt{\frac{8 \pi}{3} G \rho}
            \end{equation}
        """
        if environment.get('LOGARITHMIC_TIMESTEP'):
            self.dx = self.x * self.dy
        self.rho = rho
        self.S = S
        self.H = numpy.sqrt(8./3. * numpy.pi * rho) / CONST.M_p

        self.N_eff = (
            (rho - (numpy.pi**2 / 15. * self.T**4))
            / (7./8. * numpy.pi**2 / 15. * (self.T / 1.401)**4)
        )

        old_a = self.a
        """ Physical scale factor and temperature for convenience """
        self.a = self.x / self.m
        self.T = self.aT / self.a

        """ Time step size is inferred from the approximation of the scale factor `a`
            derivative and a definition of the Hubble parameter `H`:

            \begin{equation}
                H = \frac{\dot{a}}{a} = \frac{1}{a_{i-1}} \frac{a_i - a_{i-1}}{\Delta t}
                  = \frac{1}{\Delta t}(\frac{a_i}{a_{i-1}} -1)
            \end{equation}

            \begin{equation}
                \Delta t = \frac{1}{H} (\frac{a_i}{a_{i-1}} - 1)
            \end{equation}
        """
        # dt = (self.a / old_a - 1) / self.H
        dt = (1 - old_a / self.a) / self.H
        # dt = self.dx / self.x / self.H
        self.t += dt


"""
### Distribution functions grid

To capture non-equilibrium effects in the Early Universe, we work with particle species
distribution functions $f(\vec{p}, \vec{r}, t)$. Assuming that the Universe is homogeneous
and isotropic, we can forget dependency on the position and the momentum direction:
$f(p, t)$.

Resulting functions are sampled across a wide range of momenta. However, momentum range
cannot include both 0 momenta and very large momenta (either leads to numerical overflows
and errors).
"""


class LinearSpacedGrid(object):

    def __init__(self, MOMENTUM_SAMPLES=None, MAX_MOMENTUM=None):
        if not MAX_MOMENTUM:
            MAX_MOMENTUM = environment.get('MAX_MOMENTUM_MEV') * UNITS.MeV
        if not MOMENTUM_SAMPLES:
            MOMENTUM_SAMPLES = environment.get('MOMENTUM_SAMPLES')

        self.MIN_MOMENTUM = 0
        self.MAX_MOMENTUM = MAX_MOMENTUM
        self.BOUNDS = (self.MIN_MOMENTUM, self.MAX_MOMENTUM)
        self.MOMENTUM_SAMPLES = MOMENTUM_SAMPLES

        """
        Grid template can be copied when defining a new distribution function and is convenient to\
        calculate any _vectorized_ function over the grid. For example,

        ```python
            particle.conformal_energy(GRID.TEMPLATE)
        ```

        yields an array of particle conformal energy mapped over the `GRID`
        """
        self.TEMPLATE = numpy.linspace(self.MIN_MOMENTUM, self.MAX_MOMENTUM,
                                       num=self.MOMENTUM_SAMPLES, endpoint=True)


class LogSpacedGrid(object):

    def __init__(self, MOMENTUM_SAMPLES=None, MAX_MOMENTUM=None):
        if not MAX_MOMENTUM:
            MAX_MOMENTUM = environment.get('MAX_MOMENTUM_MEV') * UNITS.MeV
        if not MOMENTUM_SAMPLES:
            MOMENTUM_SAMPLES = environment.get('MOMENTUM_SAMPLES')

        self.MIN_MOMENTUM = 0
        self.MAX_MOMENTUM = self.MIN_MOMENTUM + MAX_MOMENTUM
        self.BOUNDS = (self.MIN_MOMENTUM, self.MAX_MOMENTUM)
        self.MOMENTUM_SAMPLES = MOMENTUM_SAMPLES

        self.TEMPLATE = self.generate_template()

    def generate_template(self):
        base = 1.2
        return (
            self.MIN_MOMENTUM
            + (self.MAX_MOMENTUM - self.MIN_MOMENTUM)
            * (base ** numpy.arange(0, self.MOMENTUM_SAMPLES, 1) - 1.)
            / (base ** (self.MOMENTUM_SAMPLES - 1.) - 1.)
        )


class HeuristicGrid(object):

    def __init__(self, M, tau, aT=1*UNITS.MeV, b=0.8, c=200):
        H = 0.5 / UNITS.s  # such that at T=1 <=> t=1
        a_max = numpy.sqrt(2 * H * b * tau)
        T_max = aT / a_max

        T = T_max
        grid = [a_max * (M + numpy.sqrt(M*T))]

        while grid[-1] > 0:
            g = grid[-1]
            T = aT * M * (aT + g*c**2 - numpy.sqrt(aT * (aT + 2*g*c**2))) / (2 * g**2 * c**2)
            grid.append(max(g - numpy.sqrt(M*T/c) * aT/T, 0))

        self.TEMPLATE = numpy.array(grid[::-1])
        self.MOMENTUM_SAMPLES = len(self.TEMPLATE)

        self.MIN_MOMENTUM = 0
        self.MAX_MOMENTUM = grid[0]
        self.BOUNDS = (self.MIN_MOMENTUM, self.MAX_MOMENTUM)


GRID = LinearSpacedGrid()


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


def cubic_interpolation(function, grid):
    """ # Cubic interpolation """

    def interpolation(x):
        index = numpy.searchsorted(grid, x)

        if index >= len(grid) - 1:
            return function[len(grid) - 1]

        if x == grid[index]:
            return function[index]

        if index == 0 or index == len(grid) - 2:
            # Determine the closest grid points
            i_low = index
            x_low = grid[i_low]

            i_high = index + 1
            x_high = grid[i_high]

            return (
                function[i_low] * (x_high - x) + function[i_high] * (x - x_low)
            ) / (x_high - x_low)

        else:
            # Otherwise, use cubic interpolation
            y0 = function[index - 1]
            y1 = function[index]
            y2 = function[index + 1]
            y3 = function[index + 2]

            a0 = -0.5*y0 + 1.5*y1 - 1.5*y2 + 0.5*y3
            a1 = y0 - 2.5*y1 + 2*y2 - 0.5*y3
            a2 = -0.5*y0 + 0.5*y2
            a3 = y1

            mu = (x - grid[index]) / (grid[index + 1] - grid[index])

            return a0*mu*mu*mu + a1*mu*mu + a2*mu + a3

    return interpolation


def binary_find(grid, x):
    temp = grid.TEMPLATE
    head = 0
    tail = grid.MOMENTUM_SAMPLES - 1

    if temp[-1] < x:
        return tail, -1

    if temp[0] > x:
        return -1, head

    if x in temp:
        pos = numpy.where(temp == x)
        return pos, pos

    return numpy.searchsorted(temp, x, side='right') - 1, numpy.searchsorted(temp, x, side='right')


def interp(particle, p, conformal_mass):
    if particle.in_equilibrium:
        return 1. / (
        numpy.exp(particle.conformal_energy(p, conformal_mass) / particle.aT)
        + particle.eta
    )

    i_lo, i_hi = binary_find(particle.grid, p)
    if i_lo == -1:
        raise ValueError("Input momentum is too small for the given grid")

    if i_hi == -1:
        return particle._distribution[-1] \
            * numpy.exp((particle.conformal_energy(particle.grid.MAX_MOMENTUM, conformal_mass) - particle.conformal_energy(p, conformal_mass)) / particle.aT)

    if i_lo == i_hi:
        return particle._distribution[i_lo]

    p_lo = particle.grid.TEMPLATE[i_lo]
    p_hi = particle.grid.TEMPLATE[i_hi]

    E_p  = particle.conformal_energy(p, conformal_mass)
    E_lo = particle.conformal_energy(p_lo, conformal_mass)
    E_hi = particle.conformal_energy(p_hi, conformal_mass)

    g_hi = particle._distribution[i_hi]
    g_lo = particle._distribution[i_lo]

    g_hi = (1. / g_hi - particle.eta)
    if g_hi > 0:
        g_hi = numpy.log(g_hi)
    else:
        return 0.

    g_lo = (1. / g_lo - particle.eta)
    if g_lo > 0:
        g_lo = numpy.log(g_lo)
    else:
        return 0.

    g = ((E_p - E_lo) * g_hi + (E_hi - E_p) * g_lo) / (E_hi - E_lo)

    g = 1. / (numpy.exp(g) + particle.eta)
    if not numpy.isfinite(g):
        return 0.

    return g
