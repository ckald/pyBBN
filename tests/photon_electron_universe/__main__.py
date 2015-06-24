"""
## Photon-electron universe test

<img src="plots.svg" width=100% />

This test checks that in the photon-electron universe:

  * $a \propto t^{1/2}$
  * $a * T$ is not conserved by a factor around `1.401`

[Log file](log.txt)

"""

import numpy

from particles import Particle
from evolution import Universe
from library.SM import particles as SMP
from common import Params, UNITS


params = Params(T=100 * UNITS.MeV,
                dx=1e-1 * UNITS.MeV)

universe = Universe(params=params, logfile="tests/photon_electron_universe/log.txt")

photon = Particle(**SMP.photon)
electron = Particle(**SMP.leptons.electron)

universe.add_particles([photon, electron])
universe.evolve()

initial_aT = universe.data['aT'][0]
print "a * T is not conserved: {}"\
    .format(any([initial_aT != value for value in universe.data['aT']]))
initial_a = universe.data['a'][len(universe.data['a'])/2]
initial_t = universe.data['t'][len(universe.data['a'])/2] / UNITS.s
last_a = universe.data['a'][-1]
last_t = universe.data['t'][-1] / UNITS.s

print """
    a scaling discrepancy is: {:.2f}%"""\
    .format(100 * (last_a / initial_a - numpy.sqrt(last_t / initial_t))
            / numpy.sqrt(last_t / initial_t))

print """
    Cosmic photon background temperature is {:.3f} times bigger than cosmic neutrinos temperature.
    Relative error is {:.3f} %""".format(universe.params.aT / UNITS.MeV,
                                         (universe.params.aT / UNITS.MeV - 1.401) / 1.401 * 100)

universe.graphics.save(__file__)

raw_input("...")
