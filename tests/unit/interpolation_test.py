import numpy
from nose import with_setup_args

from common import GRID
from particles import Particle
from library.SM import particles as SMP

from . import setup


@with_setup_args(setup)
def init_distribution_test(params):

    photon = Particle(params=params, **SMP.photon)
    neutrino = Particle(params=params, **SMP.leptons.neutrino_e)

    print neutrino._distribution - numpy.vectorize(neutrino.distribution)(GRID.TEMPLATE)

    assert all(photon._distribution == numpy.vectorize(photon.distribution)(GRID.TEMPLATE))
    assert all(neutrino._distribution == numpy.vectorize(neutrino.distribution)(GRID.TEMPLATE))


@with_setup_args(setup)
def distribution_interpolation_accuracy_test(params):

    neutrino = Particle(params=params, **SMP.leptons.neutrino_e)

    detailed_grid = numpy.linspace(GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM*2,
                                   num=GRID.MOMENTUM_SAMPLES*10, endpoint=True)

    print numpy.abs(
        neutrino.equilibrium_distribution_function(detailed_grid / neutrino.aT)
        - numpy.vectorize(neutrino.distribution)(detailed_grid)
    )

    assert all(numpy.abs(
        neutrino.equilibrium_distribution_function(detailed_grid / neutrino.aT)
        - numpy.vectorize(neutrino.distribution)(detailed_grid)
    ) < 1e-9)
