"""
## Photon universe test

<img src="plots.svg" width=100% />

This test checks that in the photon universe:

  * $a * T$ is conserved exactly.
  * $a \propto t^{1/2}$

[Log file](log.txt)

"""

import os
import numpy

from particles import Particle
from evolution import Universe
from library.SM import particles as SMP
from common import Params, UNITS


folder = os.path.join(os.path.split(__file__)[0], 'output')

T_final = 100 * UNITS.keV
params = Params(T=100 * UNITS.MeV,
                dy=0.05)

universe = Universe(params=params, folder=folder)

photon = Particle(**SMP.photon)

universe.add_particles([photon])
universe.evolve(T_final)

initial_aT = universe.data['aT'][0]
print("a * T is conserved: {}".format(all([initial_aT == value for value in universe.data['aT']])))
initial_a = universe.data['a'][int(len(universe.data)/2)]
initial_t = universe.data['t'][int(len(universe.data)/2)] / UNITS.s
last_a = universe.data['a'][-1]
last_t = universe.data['t'][-1] / UNITS.s

print("a scaling discrepancy is: {:.2f}%"
      .format(100 * (last_a / initial_a / numpy.sqrt(last_t / initial_t) - 1)))
