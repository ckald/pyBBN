import time
import numpy
from scipy import integrate
from functools import wraps
import numericalunits as nu


class UNITS:
    s = 1. / 6.58 * 1e25 / nu.GeV
    kg = 1e27 / 1.8 * nu.GeV
    m = 1e15 / 0.197 / nu.GeV
    N = 1e-5 / 8.19 * nu.GeV**2


class PARAMS:
    T = 100. * nu.MeV

    a = 1.

    t = 0. * UNITS.s
    dt = 1e-5 * UNITS.s
    t_i = 0. * UNITS.s
    t_f = 1. * UNITS.s


class GRID:
    MAX_MOMENTUM = 1e8 * nu.eV + 1e0 * nu.eV
    MIN_MOMENTUM = 1e0 * nu.eV
    MOMENTUM_STEP = 1e4 * nu.eV
    MOMENTUM_SAMPLES = int(numpy.round_((MAX_MOMENTUM - MIN_MOMENTUM) / MOMENTUM_STEP))

GRID.TEMPLATE = numpy.arange(GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM, GRID.MOMENTUM_STEP)


class CONST:
    G = 6.67 * 1e-11 * (UNITS.N / UNITS.kg**2 * UNITS.m**2)
    G_F = 1.166 * 1e-5 / nu.GeV**2


class STATISTICS:
    BOSON = "Boson"
    FERMION = "Fermion"


class REGIMES:
    RADIATION = 'radiation'
    DUST = 'dust'
    INTERMEDIATE = 'intermediate'
    NONEQ = 'non-equilibrium'


def memodict(f):
    """ Memoization decorator for a function taking a single argument """
    class memodict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret
    return memodict().__getitem__


class benchmark(object):
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, ty, val, tb):
        end = time.time()
        print("%s : %0.10f seconds" % (self.name, end-self.start))
        return False


@memodict
def momentum_to_index(p):
    return int(
        (GRID.MOMENTUM_SAMPLES * numpy.abs(p) - GRID.MIN_MOMENTUM)
        / (GRID.MAX_MOMENTUM - GRID.MIN_MOMENTUM)
    )


@memodict
def index_to_momentum(i):
    return GRID.MIN_MOMENTUM + (GRID.MAX_MOMENTUM - GRID.MIN_MOMENTUM) * i / GRID.MOMENTUM_SAMPLES


class Distributions:
    @staticmethod
    def Fermi(e):
        return 1 / (numpy.exp(e) + 1)

    @staticmethod
    def Bose(e):
        return 1 / (numpy.exp(e) - 1)

Distributions.FermiV = numpy.vectorize(Distributions.Fermi)
Distributions.BoseV = numpy.vectorize(Distributions.Bose)


def lambda_integrate(func):
    @wraps(func)
    def wrapper(*args, **kw):
        # fpp = GRID.TEMPLATE
        # for i in xrange(GRID.MOMENTUM_SAMPLES):
            # fpp[i] = func(*args, **kw)(i)
        # with benchmark("vectorize\t\t" + func.__name__ + "\t\t"):
            fpp = func(*args, **kw)(GRID.TEMPLATE)
            # print GRID.TEMPLATE / nu.eV, fpp
        # with benchmark("integrate\t\t" + func.__name__ + "\t\t"):
            inte = integrate.simps(fpp, dx=GRID.MOMENTUM_STEP)
            return inte
    return wrapper


def echo(func):
    @wraps(func)
    def wrapper(*args, **kw):
        val = func(*args, **kw)
        print val
        return val
    return wrapper


def theta(f):
    if f > 0:
        return 1
    if f == 0:
        return 0.5
    return 0
