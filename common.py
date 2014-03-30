import time
import numpy
from scipy import integrate
from functools import wraps
import numericalunits as nu
import cPickle


class UNITS:
    s = 1. / 6.58 * 1e25 / nu.GeV
    kg = 1e27 / 1.8 * nu.GeV
    m = 1e15 / 0.197 / nu.GeV
    N = 1e-5 / 8.19 * nu.GeV**2


class PARAMS:
    a_initial = 1.

    T_initial = 2.1 * nu.MeV
    T_final = 1e-3 * nu.MeV

    m = 1. * nu.MeV
    dx = 1e-2 * nu.MeV

    t = 0. * UNITS.s

    H = 0.

    regime_factor = 1e1

PARAMS.a = PARAMS.a_initial
PARAMS.x = PARAMS.a * PARAMS.m
PARAMS.T = PARAMS.T_initial
PARAMS.aT = PARAMS.a * PARAMS.T


class GRID:
    MAX_MOMENTUM = (1e8 + 1e-1) * nu.eV
    MIN_MOMENTUM = 1e-1 * nu.eV
    MOMENTUM_STEP = 0.5 * 1e7 * nu.eV
    MOMENTUM_SAMPLES = int(numpy.round_((MAX_MOMENTUM - MIN_MOMENTUM) / MOMENTUM_STEP))

GRID.TEMPLATE = numpy.linspace(GRID.MIN_MOMENTUM, GRID.MAX_MOMENTUM,
                               num=GRID.MOMENTUM_SAMPLES, endpoint=True)


class CONST:
    G = 6.67 * 1e-11 * (UNITS.N / UNITS.kg**2 * UNITS.m**2)
    G_F = 1.166 * 1e-5 / nu.GeV**2
    sin_theta_w_2 = 0.23

CONST.g_R = CONST.sin_theta_w_2
CONST.g_L = CONST.sin_theta_w_2 - 0.5


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
    index = (
        (GRID.MOMENTUM_SAMPLES * numpy.abs(p) - GRID.MIN_MOMENTUM)
        / (GRID.MAX_MOMENTUM - GRID.MIN_MOMENTUM)
    )
    if index >= GRID.MOMENTUM_SAMPLES:
        index = GRID.MOMENTUM_SAMPLES - 1
    elif index < 0:
        index = 0

    return index


@memodict
def index_to_momentum(i):
    return GRID.MIN_MOMENTUM + (GRID.MAX_MOMENTUM - GRID.MIN_MOMENTUM) * i / GRID.MOMENTUM_SAMPLES


class Distributions:
    @staticmethod
    def Fermi(e):
        return 1. / (numpy.exp(e) + 1.)

    @staticmethod
    def Bose(e):
        return 1. / (numpy.exp(e) - 1.)


Distributions.FermiV = numpy.vectorize(Distributions.Fermi, otypes=[numpy.float_])
Distributions.BoseV = numpy.vectorize(Distributions.Bose, otypes=[numpy.float_])

# numpy.seterr(all='print')


def lambda_integrate(func):
    @wraps(func)
    def wrapper(*args, **kw):
        fpp = func(*args, **kw)(GRID.TEMPLATE)
        return integrate.simps(fpp, dx=GRID.MOMENTUM_STEP)
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
        return 1.
    if f == 0:
        return 0.5
    return 0.


# import memcache
# mc = memcache.Client(['127.0.0.1:11211'], debug=0)
# import json
# import hashlib


class MemoizeMutable:
    def __init__(self, fn):
        self.fn = fn
        self.memo = {}

    def __call__(self, *args, **kwds):
        pickle = cPickle.dumps(args, 1) + cPickle.dumps(kwds, 1)
        if not pickle in self.memo:
            self.memo[pickle] = self.fn(*args, **kwds)
            # mc.set(pickle, self.fn(*args, **kwds))

        return self.memo[pickle]
        # return mc.get(pickle)


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
