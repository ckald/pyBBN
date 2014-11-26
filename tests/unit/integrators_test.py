from scipy import special, integrate

from common import integrators


def gaussian_test():
    func = lambda z: special.jn(3, z)
    adaptive_result, error = integrate.quad(func, 0, 10)
    fixed_result, _ = integrate.fixed_quad(func, 0, 10, n=integrators.GAUSS_LEGENDRE_ORDER)
    own_result = integrators.gaussian(func, 0, 10)

    print adaptive_result, error
    print fixed_result
    print own_result

    assert adaptive_result - fixed_result < error, "Gauss-Legendre quadrature order is insufficient"
    assert adaptive_result - own_result < error, "Own integrator is inaccurate"


def double_gaussian_test():
    func = lambda x, y: special.jn(x, y)
    adaptive_result, error = integrate.dblquad(func, 0, 10, lambda x: 1, lambda y: 3)
    own_result = integrators.double_gaussian(func, 0, 10, lambda x: 1, lambda y: 3)

    print adaptive_result, error
    print own_result

    assert adaptive_result - own_result < error, "Own integrator is inaccurate"
