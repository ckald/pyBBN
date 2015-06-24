import numpy
import functools

from common import GRID
from common import orthopoly


def euler_correction(y, t, f, h):
    """
    Forward Euler integration method is a most basic way to solve an ODE of the kind:

    \begin{equation}
        \frac{d y(t)}{dt} = f(t, y(t))
    \end{equation}

    Derivation:

    \begin{equation}
        \frac{d y(t)}{dt} \approx \frac{y(t) - y(t-h)}{h}
    \end{equation}
    \begin{equation}
        y_{n+1} = y_n + h f(t_n, y_n)
    \end{equation}

    Local Truncation Error (LTE):  $y(t_0+h) - y_1 = \frac{h^2}{2} y''(t_0) + O(h^3) $
    """

    return h * f(t, y)


def heun_correction(y, t, f, h):
    """
    Second order Heun's method
    """
    euler_y = euler_correction(y=y,
                               t=t,
                               f=f,
                               h=h)

    heun_y = euler_correction(y=y + euler_y,
                              t=t + h,
                              f=f,
                              h=h)

    return (euler_y + heun_y) / 2.


def implicit_euler(y, t, A, B, h):
    """
    Implicit Euler solver for ODE with a linear function.

    \begin{equation}
        \frac{d y(t)}{dt} = A(t) + B(t) y(t)
    \end{equation}

    Derivation:
    \begin{equation}
        \frac{d y(t)}{dt} \approx \frac{y(t+h) - y(t)}{h} = A(t+h) + B(t+h) y(t+h)
    \end{equation}
    \begin{equation}
        y(t+h) = \frac{y(t) + A(t+h) h}{1 - B(t+h) h}
    \end{equation}
    """

    return (y + A * h) / (1 - B * h)


ADAMS_BASHFORTH_COEFFICIENTS = {
    1: ([1.], 1.),
    2: ([-1., 3.], 2.),
    3: ([5., -16., 23.], 12.),
    4: ([-9., 37., -59., 55.], 24.),
    5: ([251., -1274., 2616., -2774., 1901.], 720.)
}


def adams_bashforth_correction(fs, h, order=None):
    if order is None:
        order = min(len(ADAMS_BASHFORTH_COEFFICIENTS), len(fs))

    bs, divider = ADAMS_BASHFORTH_COEFFICIENTS[order]

    return h * sum(b * f for b, f in zip(bs, fs[-order:])) / divider


ADAMS_MOULTON_COEFFICIENTS = {
    1: ([1.], 1.),
    2: ([1., 1.], 2.),
    3: ([-1., 8., 5.], 12.),
    4: ([1., -5., 19., 9.], 24.),
    5: ([-19., 106., -264., 646., 251.], 720.)
}


def adams_moulton_solver(y, fs, A, B, h, order=None):
    if order is None:
        order = min(len(ADAMS_MOULTON_COEFFICIENTS), len(fs))

    bs, divider = ADAMS_MOULTON_COEFFICIENTS[order]

    return (y + h * sum(b * f for b, f in zip(bs, fs[-order+1:] + [A])) / divider) /\
        (1 - h * B * bs[-1] / divider)


def integrate_1D(integrand, bounds):
    integral = gaussian(
        integrand,
        bounds[0], bounds[1]
    )
    error = numpy.nan

    return integral, error


def integrate_2D(integrand, bounds):
    integral = double_gaussian(
        integrand,
        bounds[0][0], bounds[0][1],
        bounds[1][0], bounds[1][1]
    )
    error = numpy.nan

    return integral, error


GAUSS_LEGENDRE_ORDER = 30
# points, weights = polynomial.legendre.leggauss(GAUSS_LEGENDRE_ORDER)
alpha, beta = orthopoly.rec_jacobi(GAUSS_LEGENDRE_ORDER, 0, 0)
points, weights = orthopoly.lobatto(alpha, beta, -1, 1)
grid = numpy.meshgrid(points, points)


def gaussian(f, a, b):

    sub = (b - a) / 2.
    add = (b + a) / 2.

    if sub == 0:
        return 0.

    return sub * numpy.dot(f(sub * points + add), weights)


def remap_interval(f, x, y, bounds):
    a, b, g, h = bounds

    sub_x = (b - a) / 2.
    add_x = (b + a) / 2.
    norm_x = sub_x * x + add_x

    h_x = numpy.vectorize(h)(norm_x)
    g_x = numpy.vectorize(g)(norm_x)
    sub_y = (h_x - g_x) / 2.
    add_y = (h_x + g_x) / 2.
    norm_y = sub_y * y + add_y

    return sub_x * sub_y * f(norm_x, norm_y)

remap_interval = numpy.vectorize(remap_interval, otypes=[numpy.float_], excluded=[0, 3, 'bounds'])


def double_gaussian(f, a, b, g, h):

    x, y = grid
    mesh = remap_interval(f, x, y, bounds=(a, b, g, h))
    integral = numpy.dot(numpy.transpose(weights), numpy.dot(mesh, weights))

    return integral


def lambda_integrate(bounds=(GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,)):
    """ Gaussian integration over the momentum space of the lambda function """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kw):
            fpp = func(*args, **kw)
            result, _ = integrate_1D(fpp, bounds)

            return result
        return wrapper
    return decorator
