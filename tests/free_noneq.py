from particles import Particle
from library import StandardModelParticles as SMP
from evolution import Universe
from common import UNITS, PARAMS


PARAMS.T_initial = 100 * UNITS.MeV
PARAMS.T_final = 10 * UNITS.keV
PARAMS.dx = 1e-3 * UNITS.MeV
PARAMS.infer()


photon = Particle(**SMP.photon)
neutron = Particle(**SMP.neutron)
proton = Particle(**SMP.proton)
neutrino = Particle(**SMP.neutrino_e)
electron = Particle(**SMP.electron)

Particles = [
    photon,
    neutron,
    proton,
    neutrino,
    electron
]

universe = Universe(Particles)
# universe.graphics.monitor(particles=[neutrino])
universe.evolve()

initial_aT = universe.data['aT'][0]
print "a * T is conserved: {}".format(all([initial_aT == value for value in universe.data['aT']]))
initial_a = universe.data['a'][len(universe.data['a'])/2]
initial_t = universe.data['t'][len(universe.data['a'])/2] / UNITS.s
last_a = universe.data['a'][-1]
last_t = universe.data['t'][-1] / UNITS.s

print "a scaling discrepancy is: {:.2f}%"\
    .format(100 * (last_a / initial_a / (last_t / initial_t) ** (1./2.) - 1))

universe.graphics.save(__file__)

raw_input("...")
