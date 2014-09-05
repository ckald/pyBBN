from particle import Particle
from evolution import Universe
from common import STATISTICS, UNITS


Particles = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON)
Particles.append(photon)

neutron = Particle(name='Neutron',
                   statistics=STATISTICS.FERMION,
                   mass=0.939 * UNITS.GeV)
Particles.append(neutron)

proton = Particle(name='Proton',
                  statistics=STATISTICS.FERMION,
                  mass=0.938 * UNITS.GeV)
Particles.append(proton)

neutrino = Particle(name='Neutrino',
                    statistics=STATISTICS.FERMION,
                    dof=4,
                    decoupling_temperature=3 * UNITS.MeV
                    )
Particles.append(neutrino)

electron = Particle(name='Electron',
                    mass=0.511 * UNITS.MeV,
                    statistics=STATISTICS.FERMION,
                    dof=4)
Particles.append(electron)

universe = Universe(Particles)
universe.graphics.monitor([neutrino])
universe.evolve()

initial_aT = universe.data['aT'][0]
print "a * T is conserved: {}".format(all([initial_aT == value for value in universe.data['aT']]))
initial_a = universe.data['a'][len(universe.data['a'])/2]
initial_t = universe.data['t'][len(universe.data['a'])/2] / UNITS.s
last_a = universe.data['a'][-1]
last_t = universe.data['t'][-1] / UNITS.s

print "a scaling discrepancy is: {:.2f}%"\
    .format(100 * (last_a / initial_a - (last_t / initial_t) ** (2./3.))
            / (last_t / initial_t) ** (2./3.))

universe.graphics.save(__file__)

raw_input("...")
