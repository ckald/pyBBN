from interaction import Interaction
from particles import Particle, library
from evolution import Universe
from common import CONST, UNITS, PARAMS, GRID


PARAMS.T_initial = 3 * UNITS.MeV
PARAMS.T_final = 0.075 * UNITS.MeV
PARAMS.dx = 1e-3 * UNITS.MeV
PARAMS.infer()


Particles = []
Interactions = []

photon = Particle(**library.StandardModelParticles.photon)
electron = Particle(**library.StandardModelParticles.electron)
neutrino_e = Particle(**library.StandardModelParticles.neutrino_e)
neutrino_mu = Particle(**library.StandardModelParticles.neutrino_mu)
neutrino_tau = Particle(**library.StandardModelParticles.neutrino_tau)

Particles += [photon, electron, neutrino_e, neutrino_mu, neutrino_tau]

""" \begin{align}
        \nu_e + \nu_e &\to \nu_e + \nu_e
        \\\\ \nu_e + \overline{\nu_e} &\to \nu_e + \overline{\nu_e}
    \end{align} """
neutrino_e_scattering = Interaction(
    in_particles=[neutrino_e, neutrino_e],
    out_particles=[neutrino_e, neutrino_e],
    decoupling_temperature=0 * UNITS.MeV,
    K1=(64 + 128) * CONST.G_F**2,
    K2=0.,
    order=(0, 1, 2, 3)
)
Interactions.append(neutrino_e_scattering)

""" \begin{align}
        \nu_\mu + \nu_\mu &\to \nu_\mu + \nu_\mu
        \\\\ \nu_\mu + \overline{\nu_\mu} &\to \nu_\mu + \overline{\nu_\mu}
    \end{align} """
neutrino_mu_scattering = Interaction(
    in_particles=[neutrino_mu, neutrino_mu],
    out_particles=[neutrino_mu, neutrino_mu],
    decoupling_temperature=0 * UNITS.MeV,
    K1=(64 + 128) * CONST.G_F**2,
    K2=0.,
    order=(0, 1, 2, 3)
)
Interactions.append(neutrino_mu_scattering)

neutrino_tau_scattering = Interaction(
    in_particles=[neutrino_tau, neutrino_tau],
    out_particles=[neutrino_tau, neutrino_tau],
    decoupling_temperature=0 * UNITS.MeV,
    K1=(64 + 128) * CONST.G_F**2,
    K2=0.,
    order=(0, 1, 2, 3)
)
Interactions.append(neutrino_tau_scattering)

""" \begin{align}
        \nu_e + \nu_\mu &\to \nu_e + \nu_\mu
        \\\\ \nu_e + \overline{\nu_\mu} &\to \nu_e + \overline{\nu_\mu}
    \end{align} """
neutrino_e_mu_scattering = Interaction(
    in_particles=[neutrino_e, neutrino_mu],
    out_particles=[neutrino_e, neutrino_mu],
    decoupling_temperature=0 * UNITS.MeV,
    K1=32 * CONST.G_F**2 * 2,
    K2=0.,
    order=(0, 1, 2, 3)
)
Interactions.append(neutrino_e_mu_scattering)

neutrino_mu_tau_scattering = Interaction(
    in_particles=[neutrino_tau, neutrino_mu],
    out_particles=[neutrino_tau, neutrino_mu],
    decoupling_temperature=0 * UNITS.MeV,
    K1=32 * CONST.G_F**2 * 2,
    K2=0.,
    order=(0, 1, 2, 3)
)
Interactions.append(neutrino_mu_tau_scattering)

neutrino_tau_e_scattering = Interaction(
    in_particles=[neutrino_tau, neutrino_e],
    out_particles=[neutrino_tau, neutrino_e],
    decoupling_temperature=0 * UNITS.MeV,
    K1=32 * CONST.G_F**2 * 2,
    K2=0.,
    order=(0, 1, 2, 3)
)
Interactions.append(neutrino_tau_e_scattering)

""" Add non-equilibrium parts to the distribution functions """
import numpy
neutrino_e._distribution += \
    numpy.vectorize(lambda x: 0.01 * numpy.exp(-(x/UNITS.MeV-5)**2))(GRID.TEMPLATE)
neutrino_mu._distribution += \
    numpy.vectorize(lambda x: 0.01 * numpy.exp(-(x/UNITS.MeV-5)**2))(GRID.TEMPLATE)
neutrino_tau._distribution += \
    numpy.vectorize(lambda x: 0.01 * numpy.exp(-(x/UNITS.MeV-5)**2))(GRID.TEMPLATE)

universe = Universe(Particles, Interactions)
universe.graphics.monitor(particles=[neutrino_e, neutrino_mu, neutrino_tau])
universe.evolve()

raw_input("...")
