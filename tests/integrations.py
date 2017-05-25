import numpy

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from common import theta
from common.utils import benchmark

from skmonaco import mcquad
from scipy.integrate import quad, dblquad, simps


plt.ion()


def f(x):
    return ((x/10.)**2 + numpy.exp(-(x-1)**2)) * theta(7 - x) * 1e-20


def g(x, y):
    return f(numpy.sqrt(x**2 + y**2))

print("1d integration")

x = numpy.linspace(0, 10, num=10000, endpoint=True)
y = numpy.vectorize(f)(x)
plt.plot(x, y)
plt.draw()


for order in range(1, 6):
    with benchmark("Monte-Carlo {} points".format(10**order)):
        integral, error = mcquad(f, xl=[0.], xu=[10.], npoints=10**order)
        print(integral, error, end='')


with benchmark("QUADPACK full"):
    integral, error = quad(f, 0, 10, epsabs=0, epsrel=1e-8)
    print(integral, error, end='')

with benchmark("QUADPACK full"):
    integral, error = quad(f, 0, 9, epsabs=0, epsrel=1e-8)
    print(integral, error, end='')

with benchmark("QUADPACK full"):
    integral, error = quad(f, 0, 8, epsabs=0, epsrel=1e-8)
    print(integral, error, end='')

with benchmark("QUADPACK full"):
    integral, error = quad(f, 0, 7.5, epsabs=0, epsrel=1e-8)
    print(integral, error, end='')

with benchmark("QUADPACK full"):
    integral, error = quad(f, 0, 7.25, epsabs=0, epsrel=1e-8)
    print(integral, error, end='')

with benchmark("QUADPACK full"):
    integral, error = quad(f, 0, 7.1, epsabs=0, epsrel=1e-8)
    print(integral, error, end='')

with benchmark("QUADPACK full"):
    integral, error = quad(f, 0, 7.05, epsabs=0, epsrel=1e-8)
    print(integral, error, end='')

with benchmark("QUADPACK full"):
    integral, error = quad(f, 0, 7.01, epsabs=0, epsrel=1e-8)
    print(integral, error, end='')

with benchmark("QUADPACK specific"):
    integral, error = quad(f, 0, 7, epsabs=0, epsrel=1e-8)
    print(integral, error, end='')

with benchmark("QUADPACK simpson"):
    integral = simps(y, x)
    print(integral, end='')


print("2d integration")

fig = plt.figure()
ax = Axes3D(fig)

x = y = numpy.linspace(0, 10, num=100, endpoint=True)
X, Y = numpy.meshgrid(x, y)
zs = numpy.array([g(x, y) for x, y in zip(numpy.ravel(X), numpy.ravel(Y))])
Z = zs.reshape(X.shape)

ax.plot_surface(X, Y, Z, rstride=10, cstride=10)
plt.draw()


for order in range(1, 6):
    with benchmark("Monte-Carlo {} points".format(10**order)):
        integral, error = mcquad(lambda (x, y): g(x, y),
                                 xl=[0., 0.], xu=[10., 10.], npoints=10**order)
        print(integral, error, end='')


with benchmark("QUADPACK full plane"):
    integral, error = dblquad(g, 0, 10, lambda x: 0, lambda x: 10, epsrel=1e-8, epsabs=0)
    print(integral, error, end='')


with benchmark("QUADPACK specific"):
    integral, error = dblquad(g, 0, 7, lambda x: 0, lambda x: numpy.sqrt(49-x**2),
                              epsrel=1e-8, epsabs=0)
    print(integral, error, end='')


for order in range(2, 5):

    x = y = numpy.linspace(0, 10, num=10**order, endpoint=True)
    X, Y = numpy.meshgrid(x, y)
    zs = numpy.array([g(x, y) for x, y in zip(numpy.ravel(X), numpy.ravel(Y))])
    Z = zs.reshape(X.shape)
    x = y = numpy.linspace(0, 10, num=10**order, endpoint=True)

    with benchmark("QUADPACK simpson with {} points".format(100**order)):

        integrals = []
        for i, _ in enumerate(x):
            integrals.append(simps(Z[i], y))

        integral = simps(integrals, x)
        print(integral, end='')

plt.show()
