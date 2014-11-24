""" Cython version of the distribution interpolation.

    TODO: currently massless only
    TODO: no significant performance gain """

from __future__ import division
import math
import numpy
cimport numpy
cimport cython

# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = numpy.float_
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef numpy.float_t DTYPE_t

cdef int exponential_interpolation = True

cdef DTYPE_t energy(DTYPE_t p, DTYPE_t m):
    return numpy.sqrt(p**2 + m**2)


@cython.boundscheck(False)  # turn of bounds-checking for entire function
@cython.wraparound(False)
@cython.cdivision(True)
def distribution_interpolation(numpy.ndarray[DTYPE_t, ndim=1] grid,
                               numpy.ndarray[DTYPE_t, ndim=1] distribution,
                               DTYPE_t p,
                               DTYPE_t m,
                               int eta=1,
                               DTYPE_t MIN_MOMENTUM=0., DTYPE_t MOMENTUM_STEP=0.):

    cdef DTYPE_t p_low = -1, p_high = -1, remnant

    cdef unsigned int i = 0, i_low, i_high

    remnant = (p - MIN_MOMENTUM) % MOMENTUM_STEP
    index = int((p - MIN_MOMENTUM) / MOMENTUM_STEP - remnant)

    if remnant == 0:
        return distribution[index]

    if index < 0:
        return distribution[0]

    if index >= len(grid) - 1:
        return distribution[len(grid) - 1]

    # Determine the closest grid points

    i_low = index
    p_low = grid[i_low]

    i_high = index + 1
    p_high = grid[i_high]

    cdef DTYPE_t E_p, E_low, E_high, g_high, g_low, g

    if exponential_interpolation:
        """ === Exponential interpolation === """
        E_p = energy(p, m)
        E_low = energy(p_low, m)
        E_high = energy(p_high, m)

        """
        \begin{equation}
            g = \frac{ (E_p - E_{low}) g_{high} + (E_{high} - E_p) g_{low} }\
            { (E_{high} - E_{low}) }
        \end{equation}
        """

        g_high = distribution[i_high]
        g_low = distribution[i_low]

        if g_high > 0:
            g_high = (1. / g_high - 1.)
            if g_high > 0:
                g_high = math.log(g_high)

        if g_low > 0:
            g_low = (1. / g_low - 1.)
            if g_low > 0:
                g_low = math.log(g_low)

        g = ((E_p - E_low) * g_high + (E_high - E_p) * g_low) / (E_high - E_low)

        return 1. / (numpy.exp(g) + eta)

    else:
        """ === Linear interpolation === """
        return (
            distribution[i_low] * (p_high - p) + distribution[i_high] * (p - p_low)
        ) / (p_high - p_low)
