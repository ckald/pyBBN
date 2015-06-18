"""
Non-relativistic simplifications of density, energy density and pressure calculations
"""
import math
import IntermediateParticle


name = 'dust'


def density(particle):
    """ ## Density
        \begin{equation}
            n = g \left(\frac{M T}{2 \pi}\right)^{3/2} e^{-\frac{M}{T}}
        \end{equation}
    """
    return (
        particle.dof
        * math.sqrt(particle.mass * particle.T / 2. / math.pi)**3
        * math.exp(- particle.mass / particle.T)
    )


def energy_density(particle):
    """ ## Energy density
        \begin{equation}
            \rho = n \left(M + \frac32 T\right)
        \end{equation}
    """
    return (particle.mass + 3./2. * particle.T) * density(particle)


def pressure(particle):
    """ ## Pressure
        \begin{equation}
            p = n T
        \end{equation}
    """
    return density(particle) * particle.T


# ## Master equation terms
# Dust regime does not differ from intermediate regime here.
numerator = IntermediateParticle.numerator
denominator = IntermediateParticle.denominator
