import numpy

from array import array
from scipy import special, integrate

import environment
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


def solve_adams_bashforth(f, y_0, h, t_i=0, t_f=5):
    y = array('f', [y_0])
    fs = array('f', [y_0])
    t = t_i

    while t <= t_f:
        fs.append(f(t, y[-1]))
        y.append(y[-1] + integrators.adams_bashforth_correction(fs[-3:], h))
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

    y_implicit = numpy.array(solve_implicit(0, -2.3, 1, 1, 0., 5.))
    assert all(y_implicit >= 0), "Implicit Euler method be positive"
    assert all(numpy.diff(numpy.abs(y_implicit)) <= 0), "Implicit Euler method should be stable"

    y_implicit_coarse = numpy.array(solve_implicit(0, -2.3, 1, 2, 0., 5.))
    print(y_implicit_coarse)
    assert all(y_implicit_coarse >= 0), "Implicit Euler method solution be positive"
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


# def adams_bashforth_test():

#     f = lambda t, y: -15 * y

#     exact = numpy.vectorize(lambda t: numpy.exp(-15 * t))

#     y_euler = solve_explicit(f, 1., 1./32., 0., 5.)
#     # assert not all(numpy.array(y_euler) >= 0), "Explicit Euler method should be unstable here"
#     # y_heun = solve_adams_bashforth(f, 1., 0.5, 0., 5.)
#     # assert all(numpy.diff(numpy.abs(y_heun)) <= 0), "Heun method should be stable here"

#     y_adams_bashforth = solve_adams_bashforth(f, 1., 1./32., 0., 5.)
#     print(y_euler)
#     print(y_adams_bashforth)
#     assert all((y_euler[1:-1] - exact(numpy.arange(0, 5, 1./32.)))[2:] >=
#                (y_adams_bashforth[1:-1] - exact(numpy.arange(0, 5, 1./32.)))[2:]), \
#         "Heun method should be more accurate"



def gaussian_test():
    func = lambda z: special.jn(3, z)
    adaptive_result, error = integrate.quad(func, 0, 10)
    fixed_result, _ = integrate.fixed_quad(func, 0, 10, n=environment.get('GAUSS_LEGENDRE_ORDER'))
    own_result = integrators.gaussian(func, 0, 10)

    print(adaptive_result, error)
    print(fixed_result)
    print(own_result)

    assert numpy.isclose(integrators.gaussian(func, 0, 10), -integrators.gaussian(func, 10, 0))

    assert adaptive_result - fixed_result < error, "Gauss-Legendre quadrature order is insufficient"
    assert adaptive_result - own_result < error, "Own integrator is inaccurate"


def double_gaussian_test():
    def func(x, y):
        return special.jn(x, y)
    adaptive_result, error = integrate.dblquad(func, 0, 20, lambda x: 1, lambda y: 20)
    own_result = integrators.double_gaussian(func, 0, 20, lambda x: 1, lambda y: 20)

    print(adaptive_result, error)
    print(own_result)

    assert adaptive_result - own_result < error, "Own integrator is inaccurate"


def neutron_lifetime_test():
    from kawano import q as Q, m_e
    q = Q / m_e

    def func(e):
        return e * (e - q)**2 * numpy.sqrt(e**2 - 1)

    adaptive_result, error = integrate.quad(func, 1, q)
    fixed_result, _ = integrate.fixed_quad(func, 1, q, n=environment.get('GAUSS_LEGENDRE_ORDER'))
    own_result = integrators.gaussian(func, 1, q)

    print(adaptive_result, error)
    print(fixed_result)
    print(own_result)

    assert adaptive_result - fixed_result < error, "Gauss-Legendre quadrature order is insufficient"
    assert adaptive_result - own_result < error, "Own integrator is inaccurate"
