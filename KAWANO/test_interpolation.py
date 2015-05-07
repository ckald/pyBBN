""" interpolation.f90 test suite

    WARNING: don't forget that Fortran's default array indexing is (1...size) unlike Python's
             (0...size-1)
"""

import math
from sys import float_info


from interpolation import binary_search, interp_values, log_interp_values, dist_interp_values


def test_binary_search():

    # Test a basic case
    haystack = list(range(100))
    assert binary_search(50, haystack) == (51, 51)

    # Test a case with odd len
    haystack = list(range(99))
    assert binary_search(50, haystack) == (51, 51)

    # Range test case
    assert binary_search(3.1, haystack) == (4, 5)
    assert binary_search(97.9, haystack) == (98, 99)
    assert binary_search(0.1, haystack) == (1, 2)

    # Test the case with needle outside the haystack range
    assert binary_search(-1, haystack) == (1, 1)
    assert binary_search(99, haystack) == (99, 99)

    # Corner cases
    assert binary_search(0, haystack) == (1, 1)
    assert binary_search(98, haystack) == (99, 99)


def test_binary_search_reversed():

    # Test a basic case
    haystack = list(reversed(range(100)))
    assert binary_search(50, haystack, decreasing=True) == (50, 50)

    # Test a case with odd len
    haystack = list(reversed(range(99)))
    assert binary_search(50, haystack, decreasing=True) == (49, 49)

    # Range test case
    assert binary_search(3.1, haystack, decreasing=True) == (95, 96)
    assert binary_search(97.9, haystack, decreasing=True) == (1, 2)
    assert binary_search(0.1, haystack, decreasing=True) == (98, 99)

    # Test the case with needle outside the haystack range
    assert binary_search(-1, haystack, decreasing=True) == (99, 99)
    assert binary_search(99, haystack, decreasing=True) == (1, 1)

    # Corner cases
    assert binary_search(0, haystack, decreasing=True) == (99, 99)
    assert binary_search(98, haystack, decreasing=True) == (1, 1)


def test_linear_interpolation():

    xs = [0, 1, 2, 3, 4]
    ys = [0, 0, 1, 0, 5]

    assert ys == [interp_values(x, xs, ys) for x in xs]
    assert 0 == interp_values(0.5, xs, ys)
    assert 0.5 == interp_values(1.5, xs, ys) == interp_values(2.5, xs, ys)

    assert abs(4 - interp_values(3.8, xs, ys)) < 10 * float_info.epsilon


def test_logarithmic_interpolation():

    def f(x):
        return math.exp(-math.pi*x)

    xs = range(-5, 5)
    ys = map(f, xs)

    xxs = [x / 5. for x in range(-25, 21)]
    yys = map(f, xxs)

    assert [abs(log_interp_values(x, xs, ys) - y) < 10 * float_info.epsilon
            for x, y in zip(xs, ys)]
    assert [abs(log_interp_values(x, xs, ys) - y) < 10 * float_info.epsilon
            for x, y in zip(xxs, yys)]


def test_distribution_interpolation():

    def f(x):
        return 1. / (math.exp(x / math.pi) + 1)

    xs = range(-5, 5)
    ys = map(f, xs)

    xxs = [x / 5. for x in range(-25, 21)]
    yys = map(f, xxs)

    assert [abs(dist_interp_values(x, xs, ys) - y) < 10 * float_info.epsilon
            for x, y in zip(xs, ys)]
    assert [abs(dist_interp_values(x, xs, ys) - y) < 10 * float_info.epsilon
            for x, y in zip(xxs, yys)]
