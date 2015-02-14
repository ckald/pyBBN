from particles import Particle
from evolution import Universe
from library.SM import particles as SMP, interactions as SMI


universe = Universe()
params = universe.params

Interactions = []
photon = Particle(params=params, **SMP.photon)
neutrino = Particle(params=params, **SMP.neutrino_e)
electron = Particle(params=params, **SMP.electron)

universe.particles += [photon, neutrino, electron]
universe.interactions += [
    SMI.neutrino_self_scattering(neutrino),
    SMI.neutrinos_to_electrons(neutrino=neutrino, electron=electron),
    SMI.neutrino_electron_scattering(neutrino=neutrino, electron=electron)
]

universe.graphics.monitor(particles=[neutrino])
universe.evolve()

raw_input("...")
