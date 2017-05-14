import numpy

import environment


def init_quadrature():
    global points, weights, grid
    points, weights = numpy.polynomial.legendre.leggauss(environment.get('GAUSS_LEGENDRE_ORDER'))
    grid = numpy.meshgrid(points, points)


init_quadrature()
