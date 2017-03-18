#cython: wraparound=False
#cython: boundscheck=False
#cython: nonecheck=False
#cython: cdivision=True

""" Cython version of the distribution interpolation.

    TODO: currently massless only
    TODO: no significant performance gain """

from __future__ import division
cimport cython


cdef extern from "math.h" nogil:
    double exp(double)
    double log(double)
    double sqrt(double)


cdef int exponential_interpolation = True

cdef double energy(double p, double m) nogil:
    return sqrt(p**2 + m**2)


cpdef double distribution_interpolation(double[:] grid, int grid_len,
                                        double[:] distribution,
                                        double p,
                                        double m,
                                        int eta=1) nogil:

    cdef double p_low = -1, p_high = -1

    cdef unsigned int i = 0, i_low, i_high

    i = binary_search(grid, grid_len, p)

    if grid[i] == p:
        return distribution[i]

    # Determine the closest grid points

    i_low = i
    p_low = grid[i_low]

    i_high = i + 1
    p_high = grid[i_high]

    cdef double E_p, E_low, E_high, g_high, g_low, g

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
                g_high = log(g_high)

        if g_low > 0:
            g_low = (1. / g_low - 1.)
            if g_low > 0:
                g_low = log(g_low)

        g = ((E_p - E_low) * g_high + (E_high - E_p) * g_low) / (E_high - E_low)

        return 1. / (exp(g) + eta)

    else:
        """ === Linear interpolation === """
        return (
            distribution[i_low] * (p_high - p) + distribution[i_high] * (p - p_low)
        ) / (p_high - p_low)


cdef int binary_search(double[:] grid, int size, double x) nogil:
    cdef int head = 0, tail = size - 1
    cdef int middle = (tail + head) / 2

    while grid[middle] != x and tail - head > 1:
        if grid[middle] >= x:
            head = middle
        else:
            tail = middle

        middle = (tail + head) / 2

    if grid[middle] == x:
        return middle
    if grid[head] == x:
        return head
    if grid[tail] == x:
        return tail
