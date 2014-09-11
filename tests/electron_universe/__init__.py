"""
== Electron universe test ==

<img src="plots.png" width=100% />

This test checks that in the photon universe:

  * $a \propto t^{2/3}$

[Log file](log.txt)

TODO: why is there a discrepancy $\sim 13\%$? Test if it is because of the electrons regime
TODO: why is `aT` not conserved? Does it mean that fewer particles are moving faster?

"""

from particles import Particle, STATISTICS
from evolution import Universe
from common import PARAMS, UNITS


PARAMS.T_initial = 100 * UNITS.MeV
PARAMS.T_final = 0.1 * UNITS.MeV
PARAMS.dx = 1e-2 * UNITS.MeV
PARAMS.infer()


Particles = []
electron = Particle(name='Electron',
                    mass=0.511 * UNITS.MeV,
                    statistics=STATISTICS.FERMION,
                    dof=4)
Particles.append(electron)

universe = Universe(Particles, logfile="tests/electron_universe/log.txt")
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
