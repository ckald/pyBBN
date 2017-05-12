"""
For intermediate regime equilibrium particles, density, energy density and pressure\
are obtained through integration of distribution function
"""

import numpy

import environment
from common import integrators
from common.integrators import gauss_laguerre


name = 'intermediate'


def density(particle):
    """ ## Particle density

        \begin{equation}
            g \int \frac{p^2 dp}{2 \pi^2} f\left( \frac{p}{T} \right)
        \end{equation}
    """
    density, _ = integrators.integrate_1D(
        lambda p: (
            particle.equilibrium_distribution_function(particle.energy(p) / particle.T)
            * p**2 * particle.dof / 2. / numpy.pi**2
        ), (0, 20 * particle.T)
    )
    return density


# ## Particle energy density

def energy_density_integrand(p, particle):
    """ \begin{equation}
            I_\rho = \frac{g}{2 \pi^2} p^2 E(p) \, f\left( \frac{E(p)}{T} \right)
        \end{equation}
    """
    E = particle.energy(p)
    return (
        particle.equilibrium_distribution_function(E / particle.T)
        * p**2 * E
        * particle.dof / 2. / numpy.pi**2
    )


def energy_density(particle):
    """ \begin{equation}
            \rho = \int dp I_\rho
        \end{equation}
    """
    energy_density, _ = integrators.integrate_1D(
        lambda p: energy_density_integrand(p, particle),
        (0, 20 * particle.T)
    )
    return energy_density


# ## Particle pressure

def pressure_integrand(p, particle):
    """ \begin{equation}
            I_P = \frac{g}{6 \pi^2} \frac{p^4}{E(p)} \, f\left(\frac{E(p)}{T}\right)
        \end{equation}
    """
    E = particle.energy(p)
    return (
        particle.equilibrium_distribution_function(E / particle.T)
        * p**4 / E
        * particle.dof / 6. / numpy.pi**2
    )


def pressure(particle):
    """ \begin{equation}
            P = \int dp I_P
        \end{equation}
    """
    pressure, _ = integrators.integrate_1D(
        lambda p: pressure_integrand(p, particle),
        (0, 20 * particle.T)
    )
    return pressure


def entropy(particle):
    """ ## Entropy

        \begin{equation}
            s = \frac{\rho + P}{T}
        \end{equation}
    """

    return (energy_density(particle) + pressure(particle)) / particle.params.T


# ## Master equation terms

def numerator(particle):
    """ \begin{equation}
            \frac{M^2 x}{m^2 a T} * I(2)
        \end{equation}
    """
    return (particle.mass**2 * particle.params.x
            / particle.params.m**2 / particle.params.aT
            * Int(particle, 2))


def denominator(particle):
    """ \begin{equation}
            \frac{I(4) + M_N^2 I(2)}{(a T)^2}
        \end{equation}
    """
    return (Int(particle, 4)
            + particle.conformal_mass**2 * Int(particle, 2)) / particle.params.aT**2


def Int(particle, y_power=2):
    if environment.get('LAGUERRE_GAUSS_FOR_MASSIVE_EQUILIBRIUM_PARTICLES'):
        aT = particle.params.aT
        mat = particle.conformal_mass / aT

        laguerre = (
            particle.dof / 2. / numpy.pi**2 * aT**(y_power+1) * numpy.exp(-mat)
            * gauss_laguerre.integrate_1D(lambda eps: (
                (eps+mat) * (eps * (eps + 2. * mat))**((y_power-1.)/2.)
                / (numpy.exp(-eps-mat) + particle.eta)**2
            ))[0]
        )
        return laguerre
    else:
        legendre = particle.dof / 2. / numpy.pi**2 * integrators.integrate_1D(
            lambda y: (
                y**y_power * numpy.exp(-particle.conformal_energy(y) / particle.params.aT)
                / (numpy.exp(-particle.conformal_energy(y) / particle.params.aT) + particle.eta)**2
            ),
            (particle.grid.MIN_MOMENTUM, particle.grid.MAX_MOMENTUM)
        )[0]

        return legendre
