"""
For intermediate regime equilibrium particles, density, energy density and pressure\
are obtained through integration of distribution function
"""
from scipy import integrate
import numpy
from common import GRID, PARAMS


name = 'intermediate'


def density(particle):
    """ == Particle density ==

        \begin{equation}
            g \int \frac{p^2 dp}{2 \pi^2} f\left( \frac{p}{T} \right)
        \end{equation}
    """
    density, _ = integrate.quad(
        lambda p: (
            particle.distribution_function(
                particle.energy(p) / particle.T
            ) * p**2 * particle.dof / 2. / numpy.pi**2
        ), GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
        epsrel=1e-8, epsabs=0
    )
    return density


# == Particle energy density ==

def energy_density_integrand(p, particle):
    """ \begin{equation}
            I_\rho = \frac{g}{2 \pi^2} p^2 E(p) \, f\left( \frac{E(p)}{T} \right)
        \end{equation}
    """
    E = particle.energy(p)
    return (
        particle.distribution_function(E / particle.T)
        * p**2 * E
        * particle.dof / 2. / numpy.pi**2
    )


def energy_density(particle):
    """ \begin{equation}
            \rho = \int dp I_\rho
        \end{equation}
    """
    energy_density, _ = integrate.quad(
        lambda p: energy_density_integrand(p, particle),
        GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
        epsrel=1e-8, epsabs=0
    )
    return energy_density


# == Particle pressure ==

def pressure_integrand(p, particle):
    """ \begin{equation}
            I_P = \frac{g}{6 \pi^2} \frac{p^4}{E(p)} \, f\left(\frac{E(p)}{T}\right)
        \end{equation}
    """
    E = particle.energy(p)
    return (
        particle.distribution_function(E / particle.T)
        * p**4 / E
        * particle.dof / 6. / numpy.pi**2
    )


def pressure(particle):
    """ \begin{equation}
            P = \int dp I_P
        \end{equation}
    """
    pressure, _ = integrate.quad(
        lambda p: pressure_integrand(p, particle),
        GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
        epsrel=1e-8, epsabs=0
    )
    return pressure


# == Master equation terms ==

def numerator(particle):
    """ \begin{equation}
            \frac{M^2 x}{m^2 a T} * I(2)
        \end{equation}
    """
    return particle.mass**2 * PARAMS.x / PARAMS.m**2 / PARAMS.aT * I(particle, 2)


def denominator(particle):
    """ \begin{equation}
            \frac{(I(4) + M_N^2 I(2))}{(a T)^2}
        \end{equation}
    """
    return (I(particle, 4) + particle.conformal_mass**2 * I(particle)) / PARAMS.aT**2


def I(particle, y_power=2):
    """ \begin{equation}
            I(n) = \frac{g}{2 \pi^2} \int y^n dy \frac{ e^{-\frac{E_N(y)}{a T}} } \
            { \left(e^{-\frac{E_N(y)}{a T}} + \eta \right)^2 }
        \end{equation}
    """
    return particle.dof / 2. / numpy.pi**2 * integrate.quad(
        lambda y: (
            y**y_power * numpy.exp(-particle.conformal_energy(y) / PARAMS.aT)
            / (numpy.exp(-particle.conformal_energy(y) / PARAMS.aT) + particle.eta) ** 2
        ), GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM, epsrel=1e-8, epsabs=0
    )[0]
