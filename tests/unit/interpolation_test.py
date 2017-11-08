import numpy

from common import UNITS
from particles import Particle
from library.SM import particles as SMP
from interactions.four_particle.cpp.integral import binary_find


from . import setup, with_setup_args


@with_setup_args(setup)
def init_distribution_test(params):

    photon = Particle(params=params, **SMP.photon)
    neutrino = Particle(params=params, **SMP.leptons.neutrino_e)

    neutrino.update()

    assert numpy.allclose(photon._distribution,
                          numpy.vectorize(photon.distribution)(photon.grid.TEMPLATE))
    assert numpy.allclose(neutrino._distribution,
                          numpy.vectorize(neutrino.distribution)(neutrino.grid.TEMPLATE))


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
    assert binary_find(haystack, 50.) == (50, 50)  # (head=50, tail=50)
    assert binary_find(haystack, 49.) == (49, 49)

     # Range test case
    assert binary_find(haystack, 3.1) == (3, 4)
    assert binary_find(haystack, 4.2) == (4, 5)
    assert binary_find(haystack, 96.5) == (96, 97)
    assert binary_find(haystack, 97.9) == (97, 98)
    assert binary_find(haystack, 0.1) == (0, 1)

    # Test the case with needle outside the haystack range
    assert binary_find(haystack, -1.) == (0, 0)
    assert binary_find(haystack, 990.) == (99, 99)

    # Corner cases
    assert binary_find(haystack, 0.) == (0, 0)
    assert binary_find(haystack, 98.) == (98, 98)


    # Test a case with odd len
    haystack = numpy.array(range(99), dtype=numpy.float_)
    assert binary_find(haystack, 50.) == (50, 50)
    assert binary_find(haystack, 49.) == (49, 49)

    # Range test case
    assert binary_find(haystack, 3.1) == (3, 4)
    assert binary_find(haystack, 4.2) == (4, 5)
    assert binary_find(haystack, 96.5) == (96, 97)
    assert binary_find(haystack, 97.9) == (97, 98)
    assert binary_find(haystack, 0.1) == (0, 1)

    # Test the case with needle outside the haystack range
    assert binary_find(haystack, -1.) == (0, 0)
    assert binary_find(haystack, 990.) == (98, 98)

    # Corner cases
    assert binary_find(haystack, 0.) == (0, 0)
    assert binary_find(haystack, 98.) == (98, 98)
