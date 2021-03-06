"""
Ultra-relativistic simplifications of density, energy density and pressure calculations
"""
import numpy
from common import statistics as STATISTICS

__path__ = "../.."


name = 'radiation'


def density(particle):
    """ ## Density
        \begin{equation}
            n = \frac{g}{\pi^2} (aT)^3 \zeta(3)
        \end{equation}
    """

    density = particle.T**3 * particle.dof / numpy.pi**2 * 1.2
    if particle.statistics == STATISTICS.FERMION:
        # Multiplied by $\frac34$ for fermions
        density *= 3./4.
    return density


def energy_density(particle):
    """ ## Energy density

        \begin{equation}
            \rho = \frac{\pi^2}{30} g T^4
        \end{equation}
    """

    density = particle.dof * numpy.pi**2 / 30. * particle.T**4
    if particle.statistics == STATISTICS.FERMION:
        # Multiplied by $\frac78$ for fermions
        density *= 7./8.
    return density


def pressure(particle):
    """ ## Pressure

        \begin{equation}
            P = \frac{\rho}{3}
        \end{equation}
    """
    return 1. * energy_density(particle) / 3.


def entropy(particle):
    """ ## Entropy

        \begin{equation}
            s = \frac{\rho + P}{T} = \frac43 \frac{\rho}{T}
        \end{equation}
    """

    return 4. / 3. * energy_density(particle) / particle.params.T


# ## Master equation terms

def numerator(particle):
    """ \begin{equation}0\end{equation} """
    return 0.


def denominator(particle):
    """ \begin{equation}
            \frac{2 \pi^2}{15} g (a T)^3
        \end{equation}
    """
    density = 2. * numpy.pi**2 / 15. * particle.dof * (particle.aT)**3
    if particle.statistics == STATISTICS.FERMION:
        # Multiplied by $\frac78$ for fermions
        density *= 7./8.
    return density
