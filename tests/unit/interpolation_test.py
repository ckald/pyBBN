import numpy
from nose import with_setup

from common import GRID
from particles import Particle
from library import StandardModelParticles as SMP

from . import eps, setup


@with_setup(setup)
def init_distribution_test():

    photon = Particle(**SMP.photon)
    neutrino = Particle(**SMP.neutrino_e)

    print neutrino._distribution - numpy.vectorize(neutrino.distribution)(GRID.TEMPLATE)

    assert all(photon._distribution == numpy.vectorize(photon.distribution)(GRID.TEMPLATE))
    assert all(neutrino._distribution == numpy.vectorize(neutrino.distribution)(GRID.TEMPLATE))


@with_setup(setup)
def distribution_interpolation_accuracy_test():

    neutrino = Particle(**SMP.neutrino_e)

    detailed_grid = numpy.linspace(GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM*2,
                                   num=GRID.MOMENTUM_SAMPLES*10, endpoint=True)

    print numpy.abs(
        neutrino.distribution_function(detailed_grid / neutrino.aT)
        - numpy.vectorize(neutrino.distribution)(detailed_grid)
    )

    assert all(numpy.abs(
        neutrino.distribution_function(detailed_grid / neutrino.aT)
        - numpy.vectorize(neutrino.distribution)(detailed_grid)
    ) < 1e-9)
