import numpy

import environment


def init_quadrature():
    global points, weights, grid
    points, weights = numpy.polynomial.laguerre.laggauss(environment.get('GAUSS_LAGUERRE_ORDER'))
    grid = numpy.meshgrid(points, points)


init_quadrature()


def gaussian(f):
    return numpy.dot(f(points), weights)


def integrate_1D(integrand):
    integral = gaussian(integrand)
    error = numpy.nan

    return integral, error
