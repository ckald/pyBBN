"""
== Photon universe test ==

<img src="plots.png" width=100% />

This test checks that in the photon universe:

  * $a * T$ is conserved exactly.
  * $a \propto t^{1/2}$

[Log file](log.txt)

"""

import numpy

from particles import Particle
from evolution import Universe
from library import StandardModelParticles as SMP
from common import Params, UNITS


params = Params(T_initial=100 * UNITS.MeV,
                T_final=100 * UNITS.keV,
                dx=1e-2 * UNITS.MeV)

photon = Particle(params=params, **SMP.photon)

universe = Universe(params=params, logfile="tests/photon_universe/log.txt")
universe.particles.append(photon)
universe.evolve()

initial_aT = universe.data['aT'][0]
print "a * T is conserved: {}".format(all([initial_aT == value for value in universe.data['aT']]))
initial_a = universe.data['a'][len(universe.data['a'])/2]
initial_t = universe.data['t'][len(universe.data['a'])/2] / UNITS.s
last_a = universe.data['a'][-1]
last_t = universe.data['t'][-1] / UNITS.s

print "a scaling discrepancy is: {:.2f}%"\
    .format(100 * (last_a / initial_a / numpy.sqrt(last_t / initial_t) - 1))

universe.graphics.save(__file__)

raw_input("...")
