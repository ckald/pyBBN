import numpy
from numpy import polynomial

from scipy import integrate
from skmonaco import mcquad

from common import GRID


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

    return (y + A(t + h) * h) / (1 - B(t + h) * h)


def integrate_2D(integrand, bounds, method='mcquad'):
    if method == 'dblquad':
        integral, error = integrate.dblquad(
            lambda p1, p2: integrand((p1, p2)),
            bounds[0][0], bounds[0][1],
            bounds[1][0], bounds[1][1],
            epsrel=1e-1, epsabs=0
        )
    elif method == 'mcquad':
        integral, error = mcquad(
            integrand,
            xl=[GRID.MIN_MOMENTUM, GRID.MIN_MOMENTUM],
            xu=[GRID.MAX_MOMENTUM, GRID.MAX_MOMENTUM],
            npoints=2*1e3
        )
    elif method == 'fixed':
        integral = double_gaussian(
            lambda p1, p2: integrand((p1, p2)),
            bounds[0][0], bounds[0][1],
            bounds[1][0], bounds[1][1]
        )
        error = numpy.nan

    return integral, error


GAUSS_LEGENDRE_ORDER = 40
points, weights = polynomial.legendre.leggauss(GAUSS_LEGENDRE_ORDER)
grid = numpy.meshgrid(points, points)


def gaussian(f, a, b):
    global points, weights

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

    h_x = h(norm_x)
    g_x = g(norm_x)
    sub_y = (h_x - g_x) / 2.
    add_y = (h_x + g_x) / 2.
    norm_y = sub_y * y + add_y

    return sub_x * sub_y * f(norm_x, norm_y)

remap_interval = numpy.vectorize(remap_interval, otypes=[numpy.float_], excluded=['f', 'bounds'])


def double_gaussian(f, a, b, g, h):
    global grid, weights

    x, y = grid
    mesh = remap_interval(f, x, y, bounds=(a, b, g, h))
    integral = numpy.dot(numpy.transpose(weights), numpy.dot(mesh, weights))

    return integral
