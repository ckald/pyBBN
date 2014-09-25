# -*- coding: utf-8 -*-

"""
= Common =

This file contains constants and utilities shared by all other modules in the project.
"""
import sys
import time
import numpy
import codecs
from scipy import integrate
from functools import wraps
import numericalunits as nu


# Numpy error handling behavior. Uncomment to see a notice each time you get an overflow.
# numpy.seterr(all='print')


class UNITS:

    """ == Units ==
        As we use natural units in the project, all units from `numericalunits` except energy units\
        are useless. Here some useful units are defined in terms of `GeV`s. """

    # eV = nu.eV
    eV = 1e-9

    @classmethod
    def reset_units(cls):
        UNITS.keV = UNITS.eV * 1e3
        UNITS.MeV = UNITS.keV * 1e3
        UNITS.GeV = UNITS.MeV * 1e3
        UNITS.TeV = UNITS.GeV * 1e3
        UNITS.s = 1. / 6.58 * 1e25 / UNITS.GeV
        UNITS.kg = 1e27 / 1.8 * UNITS.GeV
        UNITS.m = 1e15 / 0.197 / UNITS.GeV
        UNITS.N = 1e-5 / 8.19 * UNITS.GeV**2

UNITS.reset_units()


class CONST:
    """ === Physical constants === """

    G = 6.67 * 1e-11 * (UNITS.N / UNITS.kg**2 * UNITS.m**2)
    G_F = 1.166 * 1e-5 / UNITS.GeV**2
    sin_theta_w_2 = 0.23
    g_R = sin_theta_w_2
    g_L = sin_theta_w_2 - 0.5


class PARAMS:

    """ == Parameters ==
        Master object carrying the cosmological state of the system and initial conditions """

    # Temperature bounds define the simulations boundaries of the system
    T_initial = 10 * UNITS.MeV
    T_final = 10 * UNITS.keV

    # Arbitrary normalization of the conformal scale factor
    m = 1. * UNITS.MeV
    # Conformal scale factor step size during computations
    dx = 1e-2 * UNITS.MeV
    # Initial time
    t = 0. * UNITS.s
    # Hubble rate
    H = 0.
    # Total energy density
    rho = 0.

    @classmethod
    def infer(cls):
        # Initial scale factor - arbitrary
        PARAMS.a_initial = 1. / (PARAMS.T_initial / UNITS.MeV)
        # PARAMS.a_initial = 1.

        # Compute present-state parameters that can be inferred from the base ones
        PARAMS.a = PARAMS.a_initial
        PARAMS.x = PARAMS.a * PARAMS.m
        PARAMS.T = PARAMS.T_initial
        PARAMS.aT = PARAMS.a * PARAMS.T

        GRID.init()


class GRID:

    """ === Distribution functions grid ===

        TODO: try a log-spaced grid instead of equally-spaced

        To capture non-equilibrium effects in the Early Universe, we work with particle species \
        distribution functions $f(\vec{p}, \vec{r}, t)$. Assuming that the Universe is homogeneous\
        and isotropic, we can forget dependency on the position and the momentum direction: \
        $f(p, t)$.

        Resulting functions are sampled across a wide range of momenta. However, momentum range\
        cannot include both 0 momenta and very large momenta (either leads to numerical overflows\
        and errors).
        """

    @classmethod
    def init(cls):
        # Momentum range `(MIN_MOMENTUM, MAX_MOMENTUM)` must be divisible by `MOMENTUM_STEP`
        GRID.MIN_MOMENTUM = 0.  # 1 * UNITS.eV
        GRID.MAX_MOMENTUM = GRID.MIN_MOMENTUM + 20 * UNITS.MeV
        GRID.MOMENTUM_STEP = 1 * UNITS.MeV
        GRID.MOMENTUM_SAMPLES = int(numpy.round_((GRID.MAX_MOMENTUM - GRID.MIN_MOMENTUM)
                                    / GRID.MOMENTUM_STEP))

        # Grid template can be copied when defining a new distribution function
        GRID.TEMPLATE = numpy.linspace(GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
                                       num=GRID.MOMENTUM_SAMPLES, endpoint=True)

PARAMS.infer()


def memodict(f):
    """ Memoization decorator for a function taking a single argument """
    class memodict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret
    return memodict().__getitem__


class benchmark(object):
    """ Simple benchmarking context manager """
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, ty, val, tb):
        end = time.time()
        print("%s : %0.10f seconds" % (self.name, end-self.start))
        return False


# == Grid utilities ==


@memodict
def momentum_to_index(p):
    """ A map between values of momenta and distribution function grid points """
    index = (
        (GRID.MOMENTUM_SAMPLES * numpy.abs(p) - GRID.MIN_MOMENTUM)
        / (GRID.MAX_MOMENTUM - GRID.MIN_MOMENTUM)
    )
    if index >= GRID.MOMENTUM_SAMPLES:
        index = GRID.MOMENTUM_SAMPLES - 1
    elif index < 0:
        index = 0

    return int(index)


@memodict
def index_to_momentum(i):
    """ Inverse map: grid point index converts to momenta value """
    return GRID.MIN_MOMENTUM + (GRID.MAX_MOMENTUM - GRID.MIN_MOMENTUM) * i / GRID.MOMENTUM_SAMPLES


def theta(f):
    """ Theta $\theta$ function """
    if f > 0:
        return 1.
    if f == 0:
        return 0.5
    return 0.


def lambda_integrate(func):
    """ Scipy Simpson integration over the momentum space of the lambda function applied to the \
        grid template """
    @wraps(func)
    def wrapper(*args, **kw):
        # fpp = func(*args, **kw)(GRID.TEMPLATE)
        # result = integrate.simps(fpp, dx=GRID.MOMENTUM_STEP)

        fpp = func(*args, **kw)
        result, err = integrate.quad(fpp, GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM)
        # print result, err

        # fpp = func(*args, **kw)
        # result = integrate.romberg(fpp, GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM)

        return result
    return wrapper


def echo(func):
    @wraps(func)
    def wrapper(*args, **kw):
        val = func(*args, **kw)
        print val
        return val
    return wrapper


class MemoizeMutable:
    def __init__(self, fn):
        self.fn = fn
        self.memo = {}

    def __call__(self, *args):
        pickle = tuple(args)
        if not pickle in self.memo:
            self.memo[pickle] = self.fn(*args)

        return self.memo[pickle]


import multiprocessing
from multiprocessing import Process, Pipe
from itertools import izip


def spawn(f):
    def fun(pipe, x):
        pipe.send(f(x))
        pipe.close()
    return fun


def parmap(f, X):
    pipe = [Pipe() for x in X]
    processes = [Process(target=spawn(f), args=(c, x)) for x, (p, c) in izip(X, pipe)]
    numProcesses = len(processes)
    processNum = 0
    outputList = []
    while processNum < numProcesses:
        endProcessNum = min(processNum+multiprocessing.cpu_count(), numProcesses)
        for proc in processes[processNum:endProcessNum]:
            proc.start()
        for proc in processes[processNum:endProcessNum]:
            proc.join()
        for proc, c in pipe[processNum:endProcessNum]:
            outputList.append(proc.recv())
        processNum = endProcessNum
    return outputList


class Logger(object):
    """ Convenient double logger that redirects `stdout` and save the output also to the file """
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = codecs.open(filename, "w", encoding="utf8")

    def write(self, message, terminal=True, log=True):
        if terminal:
            self.terminal.write(message)
        if log:
            self.log.write(message.decode('utf8'))

    def __del__(self):
        sys.stdout = self.terminal


def forward_euler_integrator(y, t, f, h):
    """
    Forward Euler integration method is a most basic way to solve an ODE of the kind:

    \begin{equation}
        \frac{d y(t)}{dt} = f(t, y(t))
    \end{equation}

    Derivation:

    \begin{equation}
        \frac{d y(t)}{dt} \approx \frac{y(t) - y(t-1)}{h}
    \end{equation}
    \begin{equation}
        y_{n+1} = y_n + h f(t_n, y_n)
    \end{equation}

    Local Truncation Error (LTE):  $y(t_0+h) - y_1 = \frac{h^2}{2} y''(t_0) + O(h^3) $
    """

    return y[-1] + h * f(t[-1], y[-1])
