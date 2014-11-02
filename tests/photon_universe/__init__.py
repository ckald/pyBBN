"""
== Photon universe test ==

<img src="plots.png" width=100% />

This test checks that in the photon universe:

  * $a * T$ is conserved exactly.
  * $a \propto t^{1/2}$

[Log file](log.txt)

"""

import numpy

from particles import Particle, STATISTICS
from evolution import Universe
from common import PARAMS, UNITS


PARAMS.T_initial = 100 * UNITS.MeV
PARAMS.dx = 1e-2 * UNITS.MeV
PARAMS.infer()


Particles = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON,
                  dof=2)
Particles.append(photon)

universe = Universe(Particles, logfile="tests/photon_universe/log.txt")
universe.evolve()

initial_aT = universe.data['aT'][0]
print "a * T is conserved: {}".format(all([initial_aT == value for value in universe.data['aT']]))
initial_a = universe.data['a'][len(universe.data['a'])/2]
initial_t = universe.data['t'][len(universe.data['a'])/2] / UNITS.s
last_a = universe.data['a'][-1]
last_t = universe.data['t'][-1] / UNITS.s

print "a scaling discrepancy is: {:.2f}%"\
    .format(100 * (last_a / initial_a - numpy.sqrt(last_t / initial_t))
            / numpy.sqrt(last_t / initial_t))

universe.graphics.save(__file__)

raw_input("...")
