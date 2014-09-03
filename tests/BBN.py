from particle import Particle
from evolution import Universe
from common import STATISTICS
import numericalunits as nu


Particles = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON)
Particles.append(photon)

# neutron = Particle(name='Neutron',
#                    statistics=STATISTICS.FERMION,
#                    mass=0.939 * nu.GeV)
# Particles.append(neutron)

# proton = Particle(name='Proton',
#                   statistics=STATISTICS.FERMION,
#                   mass=0.938 * nu.GeV)
# Particles.append(proton)

# neutrino = Particle(name='Neutrino',
#                     statistics=STATISTICS.FERMION,
#                     dof=4,
#                     )
# Particles.append(neutrino)

# electron = Particle(name='Electron',
#                     mass=0.511 * nu.MeV,
#                     statistics=STATISTICS.FERMION,
#                     dof=4)
# Particles.append(electron)

universe = Universe(Particles)
universe.evolve()

for particle in Particles:
    print particle

raw_input("...")
