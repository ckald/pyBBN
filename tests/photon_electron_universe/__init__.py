"""
== Photon-electron universe test ==

<img src="plots.png" width=100% />

This test checks that in the photon-electron universe:

  * $a \propto t^{1/2}$
  * $a * T$ is not conserved by a factor around `1.401`

[Log file](log.txt)

"""

import numpy

from particles import Particle, STATISTICS
from evolution import Universe
from common import PARAMS, UNITS


PARAMS.T_initial = 100 * UNITS.MeV
PARAMS.dx = 1e-1 * UNITS.MeV
PARAMS.infer()


Particles = []
photon = Particle(name='Photon',
                  statistics=STATISTICS.BOSON,
                  dof=2)
Particles.append(photon)

electron = Particle(name='Electron',
                    mass=0.511 * UNITS.MeV,
                    statistics=STATISTICS.FERMION,
                    dof=4)
Particles.append(electron)

universe = Universe(Particles, logfile="tests/photon_electron_universe/log.txt")
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
    Relative error is {:.3f} %""".format(PARAMS.aT / UNITS.MeV,
                                         (PARAMS.aT / UNITS.MeV - 1.401) / 1.401 * 100)

universe.graphics.save(__file__)

raw_input("...")
