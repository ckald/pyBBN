#cython: wraparound=False
#cython: boundscheck=False
#cython: nonecheck=False
#cython: cdivision=True


import math
import numpy
cimport numpy
cimport cython
from libc.stdlib cimport malloc, free

cdef extern from "math.h" nogil:
    double sqrt(double)

from libcpp.vector cimport vector
from cpython cimport array
import array

numpy.import_array()


cdef struct M_t:
    int[4] order
    double K1
    double K2



cpdef double conformal_energy(double y, double mass=0):
    """ Conformal energy of the particle in comoving coordinates with evolving mass term

        \begin{equation}
            E_N = \sqrt{y^2 + (M a)^2}
        \end{equation}
    """
    if mass > 0:
        return sqrt(y**2 + mass**2)
    else:
        return abs(y)


cdef double in_bounds(double p[4], double E[4], double[:] m):
    """ $D$-functions involved in the interactions imply a cut-off region for the collision\
        integrand. In the general case of arbitrary particle masses, this is a set of \
        irrational inequalities that can hardly be solved (at least, Wolfram Mathematica does\
        not succeed in this). To avoid excessive computations, it is convenient to do an early\
        `return 0` when the particles kinematics lay out of the cut-off region """
    cdef double q1, q2, q3, q4
    q1, q2 = (p[0], p[1]) if p[0] > p[1] else (p[1], p[0])
    q3, q4 = (p[2], p[3]) if p[2] > p[3] else (p[3], p[2])

    return (E[3] >= m[3] and q1 <= q2 + q3 + q4 and q3 <= q1 + q2 + q4)


cpdef numpy.ndarray[double, ndim=1] integrand(
    double p0,
    vector[double] p1s,
    vector[double] p2s,
    int length,
    double[:] m,
    vector[M_t] Ms,
    int[:] sides, fau
):
    """
    Collision integral interior.
    """

    cdef:
        double ds, temp
        int i, j, k

    cdef:
        double p[4]
        double E[4]
        double p1, p2
        numpy.ndarray[double] integrands = numpy.zeros(length)

    for i in range(length):
        integrands[i] = 0.

        p1 = p1s[i]
        p2 = p2s[i]
        p[:] = [p0, p1, p2, 0.]
        E[3] = 0
        for j in range(3):
            E[j] = conformal_energy(p[j], m[j])
            E[3] += -sides[i] * E[i]

        E[3] *= sides[3]
        p[3] = sqrt(abs(E[3]**2 - m[3]**2))

        if not in_bounds(p, E, m):
            continue

        temp = 1.

        ds = 0.
        if p[0] != 0:
            for M in Ms:
                ds += D(p=p, E=E, m=m, K1=M.K1, K2=M.K2, order=M.order, sides=sides)
            ds = ds / p[0] / E[0]
        else:
            for M in Ms:
                ds += Db(p=p, E=E, m=m, K1=M.K1, K2=M.K2, order=M.order, sides=sides)
        temp *= ds

        # Avoid rounding errors and division by zero
        for k in range(3, 1):
            if m[k] != 0:
                temp *= p[k] / E[k]

        if temp != 0:
            integrands[i] = temp * fau(p)

    return numpy.array(integrands)


cdef double D(double p[4], double E[4], double[:] m, double K1, double K2, int order[4], int[:] sides) nogil:
    """ Dimensionality: energy """

    cdef int i, j, k, l, sksl, sisjsksl
    i, j, k, l = order[0], order[1], order[2], order[3]
    sksl = sides[k] * sides[l]
    sisjsksl = sides[i] * sides[j] * sides[k] * sides[l]

    cdef double result = 0.

    if K1 != 0:
        result += K1 * (E[0]*E[1]*E[2]*E[3] * D1(p[0], p[1], p[2], p[3]) + sisjsksl * D3(p[0], p[1], p[2], p[3]))

        result += K1 * (E[i]*E[j] * sksl * D2(p[i], p[j], p[k], p[l])
                          + E[k]*E[l] * sksl * D2(p[k], p[l], p[i], p[j]))

    if K2 != 0:
        result += K2 * m[i]*m[j] * (E[k]*E[l] * D1(p[0], p[1], p[2], p[3]) + sksl * D2(p[i], p[j], p[k], p[l]))

    return result


cdef double D1(double k1, double k2, double k3, double k4) nogil:
    """ Dimensionality: energy

        \begin{align}
            D_1(p_i, p_j, p_k, p_l) = \frac{4}{\pi} \int_0^\infty \frac{d \lambda}{\lambda^2}
            sin(p_i \lambda) sin(p_j \lambda) sin(p_k \lambda) sin(p_l \lambda)
        \end{align}

    """

    cdef double q1, q2, q3, q4
    q1, q2 = (k1, k2) if k1 >= k2 else (k2, k1)
    q3, q4 = (k3, k4) if k3 >= k4 else (k4, k3)

    if (q1 > q2 + q3 + q4) or (q3 > q2 + q1 + q4):
        return 0.

    if q1 + q2 >= q3 + q4:
        result = 0.5 * (-q1 + q2 + q3 + q4) if (q1 + q4 >= q2 + q3) else q4
    else:
        result = 0.5 * (q1 + q2 - q3 + q4) if (q1 + q4 < q2 + q3) else q2

    return result


cdef double D2(double k1, double k2, double k3, double k4) nogil:
    """ Dimensionality: energy**3

        \begin{align}
            D_2(p_i, p_j, p_k, p_l) = s_k s_l \frac{4 p_k p_l}{\pi}
            \int_0^\infty \frac{d \lambda}{\lambda^2}
            sin(p_i \lambda) sin(p_j \lambda) \\\\
             \left[ cos(p_k \lambda) - \frac{sin(p_k \lambda)}{p_k \lambda} \right]
            \left[ cos(p_l \lambda) - \frac{sin(p_l \lambda)}{p_l \lambda} \right]
        \end{align}
    """

    cdef double q1, q2, q3, q4
    q1, q2 = (k1, k2) if k1 >= k2 else (k2, k1)
    q3, q4 = (k3, k4) if k3 >= k4 else (k4, k3)

    if (q1 > q2 + q3 + q4) or (q3 > q2 + q1 + q4):
        return 0.

    if q1 + q2 >= q3 + q4:
        if q1 + q4 >= q2 + q3:
            a = q1 - q2
            result = (
                a * (a**2 - 3. * (q3**2 + q4**2)) + 2. * (q3**3 + q4**3)
            ) / 12.
        else:
            result = q4**3 / 3.
    else:
        if q1 + q4 >= q2 + q3:
            result = q2 * (3. * (q3**2 + q4**2 - q1**2) - q2**2) / 6.
        else:
            a = q1 + q2
            result = (
                a * (3. * (q3**2 + q4**2) - a**2) + 2. * (q4**3 - q3**3)
            ) / 12.

    return result


cdef double D3(double k1, double k2, double k3, double k4) nogil:
    """ Dimensionality: energy**5

        \begin{align}
            D_3(p_i, p_j, p_k, p_l) = s_i s_j s_k s_l \frac{4 p_i p_j p_k p_l}{\pi}
            \int_0^\infty \frac{d \lambda}{\lambda^2} \\\\
             \left[ cos(p_i \lambda) - \frac{sin(p_i \lambda)}{p_i \lambda} \right]
             \left[ cos(p_j \lambda) - \frac{sin(p_j \lambda)}{p_j \lambda} \right] \\\\
             \left[ cos(p_k \lambda) - \frac{sin(p_k \lambda)}{p_k \lambda} \right]
            \left[ cos(p_l \lambda) - \frac{sin(p_l \lambda)}{p_l \lambda} \right]
        \end{align}
    """

    cdef double q1, q2, q3, q4
    q1, q2 = (k1, k2) if k1 >= k2 else (k2, k1)
    q3, q4 = (k3, k4) if k3 >= k4 else (k4, k3)

    if (q1 > q2 + q3 + q4) or (q3 > q2 + q1 + q4):
        return 0.

    if q1 + q2 >= q3 + q4:
        if q1 + q4 >= q2 + q3:
            result = (
                q1**5 - q2**5 - q3**5 - q4**5                                               # E**5
                + 5. * (
                    q1**2 * q2**2 * (q2 - q1)                                               # E**5
                    + q3**2 * (q2**3 - q1**3 + (q2**2 + q1**2) * q3)                        # E**5
                    + q4**2 * (q2**3 - q1**3 + q3**3 + (q1**2 + q2**2 + q3**2) * q4)        # E**5
                )
            ) / 60.
        else:
            result = q4**3 * (5. * (q1**2 + q2**2 + q3**2) - q4**2) / 30.
    else:
        if q1 + q4 >= q2 + q3:
            result = q2**3 * (5. * (q1**2 + q3**2 + q4**2) - q2**2) / 30.
        else:
            result = (
                q3**5 - q4**5 - q1**5 - q2**5
                + 5. * (
                    q3**2 * q4**2 * (q4 - q3)
                    + q1**2 * (q4**3 - q3**3 + (q4**2 + q3**2) * q1)
                    + q2**2 * (q4**3 - q3**3 + q1**3 + (q1**2 + q3**2 + q4**2) * q2)
                )
            ) / 60.

    return result


cdef double Db(double p[4], double E[4], double[:] m, double K1, double K2, int order[4], int[:] sides) nogil:
    """ Dimensionality: energy """

    cdef int i, j, k, l, sisj, sksl
    i, j, k, l = order[0], order[1], order[2], order[3]
    sisj = sides[i] * sides[j]
    sksl = sides[k] * sides[l]

    cdef double result, subresult
    result = 0.

    if K1 != 0:
        subresult = E[1]*E[2]*E[3] * Db1(p[1], p[2], p[3])

        if i * j == 0:
            subresult += sisj * E[i+j] * Db2(p[i+j], p[k], p[l])
        elif k * l == 0:
            subresult += sksl * E[k+l] * Db2(p[i], p[j], p[k+l])

        result += K1 * subresult

    if K2 != 0:
        subresult = 0

        if i * j == 0:
            subresult += m[i+j] * (E[k] * E[l] * Db1(p[1], p[2], p[3]) + sksl * Db2(p[i+j], p[k], p[l]))
        elif k * l == 0:
            subresult += m[i] * m[j] * m[k+l] * Db1(p[1], p[2], p[3])

        result += K2 * subresult

    return result


cdef double Db1(double q2, double q3, double q4) nogil:
    if (q2 + q3 > q4) and (q2 + q4 > q3) and (q3 + q4 > q2):
        return 1.
    return 0


cdef double Db2(double q2, double q3, double q4) nogil:
    if (q2 + q3 > q4) and (q2 + q4 > q3) and (q3 + q4 > q2):
        return 0.5 * (q3**2 + q4**2 - q2**2)
    return 0.
