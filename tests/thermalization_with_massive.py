from particles import Particle
from library import StandardModelParticles as SMP, StandardModelInteractions as SMI
from evolution import Universe
from common import CONST, UNITS, PARAMS


PARAMS.T_initial = 10. * UNITS.MeV
PARAMS.T_final = 0.050 * UNITS.MeV
PARAMS.dx = 1e-6 * UNITS.MeV
PARAMS.infer()


Particles = []
Interactions = []

photon = Particle(**SMP.photon)
electron = Particle(**SMP.electron)
neutrino_e = Particle(**SMP.neutrino_e)
neutrino_mu = Particle(**SMP.neutrino_mu)

neutrino_e.decoupling_temperature = PARAMS.T_initial

Particles += [
    photon,
    electron,
    neutrino_e,
    neutrino_mu,
]

Interactions += [
    # SMI.neutrino_self_scattering(neutrino_e),
    # SMI.neutrino_self_scattering(neutrino_mu),
    # SMI.neutrino_inter_scattering(neutrino_e, neutrino_mu),
    SMI.neutrinos_to_electrons(neutrino=neutrino_e, electron=electron, g_L=CONST.g_R+0.5),
    # SMI.neutrinos_to_electrons(neutrino=neutrino_mu, electron=electron, g_L=CONST.g_R-0.5),
    # SMI.neutrino_electron_scattering(neutrino=neutrino_e, electron=electron, g_L=CONST.g_R+0.5),
    # SMI.neutrino_electron_scattering(neutrino=neutrino_mu, electron=electron, g_L=CONST.g_R-0.5),
]

universe = Universe(Particles, Interactions)
universe.graphics.monitor(particles=[
    neutrino_e,
    neutrino_mu,
    electron,
    photon
])


universe.evolve()

raw_input("...")
