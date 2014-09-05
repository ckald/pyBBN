"""
== Neutron universe test ==

<img src="plots.png" width=100% />

This test checks that in the photon universe:

  * $a \propto t^{2/3}$

[Log file](log.txt)

TODO: why is there a discrepancy $\sim 15 \%$?
TODO: why is `aT` not conserved?
TODO: why scale factor law changes over time?

"""

from particles import Particle, STATISTICS
from evolution import Universe
from common import PARAMS, UNITS


PARAMS.T_initial = 1000 * UNITS.MeV
PARAMS.T_final = 100 * UNITS.MeV
PARAMS.dx = 1e-4 * UNITS.MeV
PARAMS.infer()


Particles = []
neutron = Particle(name='Neutron',
                   statistics=STATISTICS.FERMION,
                   mass=0.939 * UNITS.GeV,
                   dof=4)
Particles.append(neutron)

universe = Universe(Particles, logfile="tests/neutron_universe/log.txt")
universe.evolve(dx=PARAMS.dx, T_final=PARAMS.T_final)

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
