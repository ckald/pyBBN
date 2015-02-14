from particles import Particle
from library.SM import particles as SMP, interactions as SMI
from evolution import Universe
from common import CONST, UNITS, PARAMS


PARAMS.T_initial = 5. * UNITS.MeV
PARAMS.T_final = 0.075 * UNITS.MeV
PARAMS.dx = 1e-4 * UNITS.MeV
PARAMS.infer()


Particles = []
Interactions = []

photon = Particle(**SMP.photon)
electron = Particle(**SMP.electron)
neutrino_e = Particle(**SMP.neutrino_e)
neutrino_mu = Particle(**SMP.neutrino_mu)

# Massive tau neutrino model
massive_tau = SMP.neutrino_tau
massive_tau['mass'] = 20 * UNITS.MeV
neutrino_tau = Particle(**massive_tau)

Particles += [
    photon,
    electron,
    neutrino_e,
    neutrino_mu,
    neutrino_tau,
]

Interactions += [
    SMI.neutrino_self_scattering(neutrino_e),
    SMI.neutrino_self_scattering(neutrino_mu),
    SMI.neutrino_self_scattering(neutrino_tau),
    SMI.neutrino_inter_scattering(neutrino_e, neutrino_mu),
    SMI.neutrino_inter_scattering(neutrino_mu, neutrino_tau),
    SMI.neutrino_inter_scattering(neutrino_tau, neutrino_e),
    SMI.neutrinos_to_electrons(neutrino=neutrino_e, electron=electron, g_L=CONST.g_R+0.5),
    SMI.neutrinos_to_electrons(neutrino=neutrino_mu, electron=electron, g_L=CONST.g_R-0.5),
    SMI.neutrinos_to_electrons(neutrino=neutrino_tau, electron=electron, g_L=CONST.g_R-0.5),
    SMI.neutrino_electron_scattering(neutrino=neutrino_e, electron=electron, g_L=CONST.g_R+0.5),
    SMI.neutrino_electron_scattering(neutrino=neutrino_mu, electron=electron, g_L=CONST.g_R-0.5),
    SMI.neutrino_electron_scattering(neutrino=neutrino_tau, electron=electron, g_L=CONST.g_R-0.5),
]

universe = Universe(Particles, Interactions)
universe.graphics.monitor(particles=[
    neutrino_e,
    neutrino_mu,
    neutrino_tau
])


universe.evolve()

raw_input("...")
