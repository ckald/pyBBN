from interaction import Interaction
from particle import Particle
from evolution import Universe
from common import STATISTICS, CONST, UNITS


Particles = []
Interactions = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON)
Particles.append(photon)

# neutron = Particle(name='Neutron',
#                    statistics=STATISTICS.FERMION,
#                    mass=0.939 * UNITS.GeV,
#                    decoupling_temperature=1.3 * UNITS.MeV)
# Particles.append(neutron)

# proton = Particle(name='Proton',
#                   statistics=STATISTICS.FERMION,
#                   mass=0.938 * UNITS.GeV,
#                   )
# Particles.append(proton)

neutrino = Particle(name='Neutrino',
                    statistics=STATISTICS.FERMION,
                    dof=4,
                    decoupling_temperature=2 * UNITS.MeV
                    )
Particles.append(neutrino)

electron = Particle(name='Electron',
                    mass=0.511 * UNITS.MeV,
                    statistics=STATISTICS.FERMION,
                    dof=4)
Particles.append(electron)

neutrino_scattering = Interaction(
    in_particles=[neutrino, neutrino],
    out_particles=[neutrino, neutrino],
    decoupling_temperature=0 * UNITS.MeV,
    K1=192 * CONST.G_F**2,
    K2=0
)
Interactions.append(neutrino_scattering)

neutrino_annihilation = Interaction(
    in_particles=[neutrino, neutrino],
    out_particles=[electron, electron],
    decoupling_temperature=0 * UNITS.MeV,

    K1=128 * CONST.G_F**2 * (CONST.g_L**2 + CONST.g_L**2),
    K2=(-1) * 128 * CONST.G_F**2 * CONST.g_L * CONST.g_R * electron.mass**2
)
Interactions.append(neutrino_annihilation)

neutrino_electron_scattering = Interaction(
    in_particles=[neutrino, electron],
    out_particles=[neutrino, electron],
    decoupling_temperature=0 * UNITS.MeV,

    K1=128 * CONST.G_F**2 * 0.5,
    K2=128 * CONST.G_F**2 * CONST.g_L * CONST.g_R * electron.mass**2
)
Interactions.append(neutrino_electron_scattering)

universe = Universe(Particles, Interactions)
universe.graphics.monitor(particles=[neutrino])
universe.evolve()

for particle in Particles:
    print particle

raw_input("...")
