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


@cython.boundscheck(False)  # turn of bounds-checking for entire function
def distribution_interpolation(numpy.ndarray[DTYPE_t, ndim=1] grid,
                               numpy.ndarray[DTYPE_t, ndim=1] distribution,
                               DTYPE_t p,
                               # energy_normalized=lambda x: 0,
                               DTYPE_t eta=1.):
    cdef DTYPE_t p_low = -1, p_high = -1

    cdef unsigned int i = 0, i_low, i_high
    for point in grid:
        if point == p:
            return distribution[i]
        elif point < p:
            p_low = point
            i_low = i
        elif point > p:
            p_high = point
            i_high = i
            break
        i += 1

    if p_low == -1:
        raise Exception("Outside of interpolated range: {}".format(p))

    cdef DTYPE_t E_p, E_low, E_high, g_high, g_low, g

    if exponential_interpolation:
        # E_p = energy_normalized(p)
        # E_low = energy_normalized(p_low)
        # E_high = energy_normalized(p_high)
        E_p = p
        E_low = p_low
        E_high = p_high

        """
        \begin{equation}
            g = \frac{ (E_p - E_low) g_high + (E_high - E_p) g_low }{ (E_high - E_low) }
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
        return (
            distribution[i_low] * (p_high - p) + distribution[i_high] * (p - p_low)
        ) / (p_high - p_low)
