from interaction import Interaction
from particle import Particle
from evolution import Universe
from common import STATISTICS, CONST, UNITS


Particles = []
Interactions = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON)
Particles.append(photon)

neutron = Particle(name='Neutron',
                   statistics=STATISTICS.FERMION,
                   mass=0.939 * UNITS.GeV,
                   decoupling_temperature=1.3 * UNITS.MeV)
Particles.append(neutron)

proton = Particle(name='Proton',
                  statistics=STATISTICS.FERMION,
                  mass=0.938 * UNITS.GeV,
                  )
Particles.append(proton)

neutrino = Particle(name='Neutrino',
                    statistics=STATISTICS.FERMION,
                    dof=4,
                    )
Particles.append(neutrino)

electron = Particle(name='Electron',
                    mass=0.511 * UNITS.MeV,
                    statistics=STATISTICS.FERMION,
                    dof=4)
Particles.append(electron)

neutron_decay = Interaction(
    in_particles=[neutron],
    out_particles=[proton, neutrino, electron],
    decoupling_temperature=0 * UNITS.MeV,

    K1=64 * CONST.G_F**2,
    K2=0
)
Interactions.append(neutron_decay)

universe = Universe(Particles, Interactions)
universe.graphics.monitor(particles=[neutron])
universe.evolve()

for particle in Particles:
    print particle

raw_input("...")
