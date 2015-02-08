import numpy

from array import array
from scipy import special, integrate

from common import integrators

from . import eps


def solve_explicit(f, y_0, h, t_i=0., t_f=5):
    y = array('f', [y_0])
    t = t_i

    while t <= t_f:
        y.append(y[-1] + integrators.euler_correction(y=y[-1], t=t, f=f, h=h))
        t += h

    return y


def solve_implicit(A, B, y_0, h, t_i=0., t_f=5):
    y = array('f', [y_0])
    t = t_i

    while t <= t_f:
        y.append(integrators.implicit_euler(y=y[-1], t=t, A=A, B=B, h=h))
        t += h

    return y


def solve_heun(f, y_0, h, t_i=0., t_f=5):
    y = array('f', [y_0])
    t = t_i

    while t <= t_f:
        y.append(y[-1] + integrators.heun_correction(y=y[-1], t=t, f=f, h=h))
        t += h

    return y


def explicit_euler_test():

    f = lambda t, y: -2.3 * y

    y_unstable = solve_explicit(f, 1., 1., 0., 5.)
    assert all(numpy.diff(numpy.abs(y_unstable)) >= 0),\
        "Explicit Euler method should be unstable here"
    y_stable = solve_explicit(f, 1., 0.5, 0., 5.)
    assert all(numpy.diff(numpy.abs(y_stable)) <= 0),\
        "Explicit Euler method should be stable here"

    assert y_stable[-1] < eps


def implicit_euler_test():

    y_implicit = solve_implicit(0, -2.3, 1, 1, 0., 5.)
    assert y_implicit >= 0, "Implicit Euler method be positive"
    assert all(numpy.diff(numpy.abs(y_implicit)) <= 0), "Implicit Euler method should be stable"

    y_implicit_coarse = solve_implicit(0, -2.3, 1, 2, 0., 5.)
    print y_implicit_coarse
    assert y_implicit_coarse >= 0, "Implicit Euler method solution be positive"
    assert all(numpy.diff(numpy.abs(y_implicit_coarse)) <= 0), \
        "Implicit Euler method should be stable"


def heun_test():

    f = lambda t, y: -2.3 * y

    exact = numpy.vectorize(lambda t: numpy.exp(-2.3 * t))

    y_euler = solve_explicit(f, 1., 0.5, 0., 5.)
    assert not all(numpy.array(y_euler) >= 0), "Explicit Euler method should be unstable here"
    y_heun = solve_heun(f, 1., 0.5, 0., 5.)
    assert all(numpy.diff(numpy.abs(y_heun)) <= 0), "Heun method should be stable here"

    y_heun_detailed = solve_heun(f, 1., 0.05, 0., 5.)
    assert all((y_heun[1:-1] - exact(numpy.arange(0, 5, 0.5)))[2:] >=
               (y_heun_detailed[1:-1] - exact(numpy.arange(0, 5, 0.05)))[::10][2:]), \
        "Heun method should be more accurate"


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
