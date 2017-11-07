import numpy
from scipy.integrate import simps

from common.integrators.gauss_legendre import weights, grid


def integrate_2D_simpsons(integrands, bounds, grids):
    integral_1 = simps(simps(integrands[0], grids[0]), grids[1])
    integral_f = simps(simps(integrands[1], grids[0]), grids[1])
    return integral_1, integral_f


def integrate_2D(integrands, bounds):
    integrals = double_gaussian(
        integrands,
        bounds[0][0], bounds[0][1],
        bounds[1][0], bounds[1][1]
    )

    return integrals


def remap_interval(fs, x, y, bounds):
    a, b, g, h = bounds

    sub_x = (b - a) / 2.
    add_x = (b + a) / 2.
    norm_x = sub_x * x + add_x

    h_x = h(norm_x)
    g_x = g(norm_x)
    sub_y = (h_x - g_x) / 2.
    add_y = (h_x + g_x) / 2.
    norm_y = sub_y * y + add_y

    f_1, f_f = fs(norm_x, norm_y)
    return sub_x * sub_y * f_1, sub_x * sub_y * f_f


def double_gaussian(fs, a, b, g, h):

    x, y = grid
    mesh_1, mesh_f = remap_interval(fs, x, y, bounds=(a, b, g, h))
    integral_1 = numpy.dot(numpy.transpose(weights), numpy.dot(mesh_1, weights))
    integral_f = numpy.dot(numpy.transpose(weights), numpy.dot(mesh_f, weights))

    return integral_1, integral_f
