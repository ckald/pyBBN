import numpy

from common import UNITS
from particles import Particle
from library.SM import particles as SMP
from interactions.four_particle.integral import binary_search

from . import setup, with_setup_args


@with_setup_args(setup)
def init_distribution_test(params):

    photon = Particle(params=params, **SMP.photon)
    neutrino = Particle(params=params, **SMP.leptons.neutrino_e)

    neutrino.update()

    print neutrino._distribution - numpy.vectorize(neutrino.distribution)(neutrino.grid.TEMPLATE)

    assert all(photon._distribution == numpy.vectorize(photon.distribution)(photon.grid.TEMPLATE))
    assert all(neutrino._distribution == numpy.vectorize(neutrino.distribution)(neutrino.grid.TEMPLATE))


@with_setup_args(setup)
def distribution_interpolation_accuracy_test(params):
    eps = 1e-8

    neutrino = Particle(params=params, **SMP.leptons.neutrino_e)
    neutrino.update()

    detailed_grid = numpy.linspace(neutrino.grid.MIN_MOMENTUM, neutrino.grid.MAX_MOMENTUM*2,
                                   num=neutrino.grid.MOMENTUM_SAMPLES*10, endpoint=True)

    x = numpy.abs(
        neutrino.equilibrium_distribution_function(detailed_grid / neutrino.aT)
        / numpy.vectorize(neutrino.distribution)(detailed_grid)
        - 1
    )

    assert all(x < eps), ("Precision is worse than 1e-8 for momenta {}"
                          .format(detailed_grid[x >= eps] / UNITS.MeV))


def binary_search_test():

    # Test a basic case
    haystack = numpy.array(range(100), dtype=numpy.float_)
    length = len(haystack)
    assert binary_search(haystack, length, 50.) == 50  # (51, 51)

    # Test a case with odd len
    haystack = numpy.array(range(99), dtype=numpy.float_)
    assert binary_search(haystack, length, 50.) == 50  # (51, 51)

    # Range test case
    assert binary_search(haystack, length, 3.1) == 3  # (4, 5)
    assert binary_search(haystack, length, 97.9) == 97  # (98, 99)
    assert binary_search(haystack, length, 0.1) == 0  # (1, 2)

    # Test the case with needle outside the haystack range
    assert binary_search(haystack, length, -1.) == 0  # (1, 1)
    assert binary_search(haystack, length, 99.) == 98  # (99, 99)

    # Corner cases
    assert binary_search(haystack, length, 0.) == 0  # (1, 1)
    assert binary_search(haystack, length, 98.) == 98  # (99, 99)
