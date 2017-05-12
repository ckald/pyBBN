# cython: wraparound=False
# cython: boundscheck=False
# cython: nonecheck=False
# cython: cdivision=True
# cython: initializedcheck=False


import math
import numpy
cimport numpy
cimport cython

cdef extern from "math.h" nogil:
    double sqrt(double)
    double fabs(double)

from libcpp.vector cimport vector
from cpython cimport array
import array

numpy.import_array()


cdef struct M_t:
    int[4] order
    double K1
    double K2


cdef struct grid_t:
    double[:] grid
    double[:] distribution
    int size

cdef struct particle_t:
    int eta
    double m
    grid_t grid
    int in_equilibrium
    double aT

cdef struct reaction_t:
    particle_t specie
    int side



""" ## $\mathcal{F}(f_\alpha)$ functional """

""" ### Naive form

    \begin{align}
        \mathcal{F} &= (1 \pm f_1)(1 \pm f_2) f_3 f_4 - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
        \\\\ &= \mathcal{F}_B + \mathcal{F}_A
    \end{align}
"""

cdef double F_A(vector[reaction_t] reaction, double[4] f, int skip_index=-1) nogil:
    """
    Forward reaction distribution functional term

    \begin{equation}
        \mathcal{F}_A = - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
    \end{equation}

    :param skip_index: Particle to skip in the expression
    """
    cdef:
        int i
        double temp

    temp = -1.

    for i in range(4):
        if i != skip_index:
            if reaction[i].side == -1:
                temp *= f[i]
            else:
                temp *= 1. - reaction[i].specie.eta * f[i]

    return temp

cdef double F_B(vector[reaction_t] reaction, double[4] f, int skip_index=-1) nogil:
    """
    Backward reaction distribution functional term

    \begin{equation}
        \mathcal{F}_B = f_3 f_4 (1 \pm f_1) (1 \pm f_2)
    \end{equation}

    :param skip_index: Particle to skip in the expression
    """
    cdef:
        int i
        double temp

    temp = 1.

    for i in range(4):
        if i != skip_index:
            if reaction[i].side == 1:
                temp *= f[i]
            else:
                temp *= 1. - reaction[i].specie.eta * f[i]

    return temp

"""
### Linearized in $\, f_1$ form

\begin{equation}
    \mathcal{F}(f) = f_3 f_4 (1 \pm f_1) (1 \pm f_2) - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
\end{equation}

\begin{equation}
    \mathcal{F}(f) = f_1 (\mp f_3 f_4 (1 \pm f_2) - f_2 (1 \pm f_3) (1 \pm f_4)) \
    + f_3 f_4 (1 \pm f_2)
\end{equation}

\begin{equation}
    \mathcal{F}(f) = \mathcal{F}_B^{(1)} + f_1 (\mathcal{F}_A^{(1)} \pm_1 \mathcal{F}_B^{(1)})
\end{equation}

$^{(i)}$ in $\mathcal{F}^{(i)}$ means that the distribution function $f_i$ was omitted in the\
corresponding expression. $\pm_j$ represents the $\eta$ value of the particle $j$.
"""
cdef double F_f(vector[reaction_t] reaction, double[4] f) nogil:
    """ Variable part of the distribution functional """
    return F_A(reaction, f, 0) - reaction[0].specie.eta * F_B(reaction, f, 0)

cdef double F_1(vector[reaction_t] reaction, double[4] f) nogil:
    """ Constant part of the distribution functional """
    return F_B(reaction, f, 0)



cpdef double conformal_energy(double y, double mass=0) nogil:
    """ Conformal energy of the particle in comoving coordinates with evolving mass term

        \begin{equation}
            E_N = \sqrt{y^2 + (M a)^2}
        \end{equation}
    """
    if mass > 0:
        return sqrt(y**2 + mass**2)
    else:
        return fabs(y)


cdef double in_bounds(double p[4], double E[4], double m[4]) nogil:
    """ $D$-functions involved in the interactions imply a cut-off region for the collision\
        integrand. In the general case of arbitrary particle masses, this is a set of \
        irrational inequalities that can hardly be solved (at least, Wolfram Mathematica does\
        not succeed in this). To avoid excessive computations, it is convenient to do an early\
        `return 0` when the particles kinematics lay out of the cut-off region """
    cdef double q1, q2, q3, q4
    q1, q2 = (p[0], p[1]) if p[0] > p[1] else (p[1], p[0])
    q3, q4 = (p[2], p[3]) if p[2] > p[3] else (p[3], p[2])

    return (E[3] >= m[3] and q1 <= q2 + q3 + q4 and q3 <= q1 + q2 + q4)


cpdef integrand(
    double p0, vector[double] p1s, vector[double] p2s, int length,
    vector[reaction_t] reaction, vector[M_t] Ms
):
    """
    Collision integral interior.
    """

    cdef:
        double ds, temp
        int i, j, k

        double p[4]
        double E[4]
        double m[4]
        int sides[4]

        double[4] f

        double p1, p2
        numpy.ndarray[double] integrands_1 = numpy.zeros([length])
        numpy.ndarray[double] integrands_f = numpy.zeros([length])

    for i in range(4):  # [0, 1, 2, 3]
        sides[i] = reaction[i].side
        m[i] = reaction[i].specie.m

    for i in range(length):
        p1 = p1s[i]
        p2 = p2s[i]
        p[:] = [p0, p1, p2, 0.]
        E[3] = 0.
        for j in range(3):  # [0, 1, 2]
            E[j] = conformal_energy(p[j], m[j])
            E[3] += sides[j] * E[j]

        E[3] *= -sides[3]

        if E[3] < m[3]:
            continue
        p[3] = sqrt(E[3]**2 - m[3]**2)

        if not in_bounds(p, E, m):
            continue

        temp = 1.

        # Avoid rounding errors and division by zero
        for k in range(1, 3):  # [1, 2]
            if m[k] != 0.:
                temp *= p[k] / E[k]

        if temp == 0.:
            continue

        ds = 0.
        if p[0] != 0.:
            for M in Ms:
                ds += D(p=p, E=E, m=m, K1=M.K1, K2=M.K2, order=M.order, sides=sides)
            ds /= p[0] * E[0]
        else:
            for M in Ms:
                ds += Db(p=p, E=E, m=m, K1=M.K1, K2=M.K2, order=M.order, sides=sides)
        temp *= ds

        if temp == 0.:
            continue

        # The distribution function of the main particle is not required here
        f[0] = -1
        for k in range(1, 4):  # [1, 2, 3]
            if reaction[k].specie.in_equilibrium:
                f[k] = 1. / (
                    exp(conformal_energy(p[k], reaction[k].specie.m) / reaction[k].specie.aT)
                    + reaction[k].specie.eta
                )
            else:
                f[k] = distribution_interpolation(
                    reaction[k].specie.grid.grid,
                    reaction[k].specie.grid.size,
                    reaction[k].specie.grid.distribution,
                    p[k], reaction[k].specie.m, reaction[k].specie.eta
                )

        integrands_1[i] = temp * F_1(reaction, f)
        integrands_f[i] = temp * F_f(reaction, f)

    return integrands_1, integrands_f


cdef double D(double p[4], double E[4], double m[4], double K1, double K2, int order[4], int sides[4]) nogil:
    """ Dimensionality: energy """

    cdef int i, j, k, l, sisj, sksl, sisjsksl
    i, j, k, l = order[0], order[1], order[2], order[3]
    sisj = sides[i] * sides[j]
    sksl = sides[k] * sides[l]
    sisjsksl = sides[i] * sides[j] * sides[k] * sides[l]

    cdef double result = 0.

    if K1 != 0.:
        result += K1 * (E[0]*E[1]*E[2]*E[3] * D1(p[0], p[1], p[2], p[3]) + sisjsksl * D3(p[0], p[1], p[2], p[3]))

        result += K1 * (E[i]*E[j] * sksl * D2(p[i], p[j], p[k], p[l])
                        + E[k]*E[l] * sisj * D2(p[k], p[l], p[i], p[j]))

    if K2 != 0.:
        result += K2 * m[i]*m[j] * (E[k]*E[l] * D1(p[0], p[1], p[2], p[3]) + sksl * D2(p[i], p[j], p[k], p[l]))

    return result


def Dpy(p, E, m, K1, K2, order, sides):
    cdef array.array[double] cp = array.array('d', p)
    cdef array.array[double] cE = array.array('d', E)
    cdef array.array[double] cm = array.array('d', m)
    cdef array.array[int] corder = array.array('i', order)
    cdef array.array[int] csides = array.array('i', sides)
    return D(cp.data.as_doubles, cE.data.as_doubles, cm.data.as_doubles, K1, K2, corder.data.as_ints, csides.data.as_ints)


cpdef double D1(double k1, double k2, double k3, double k4) nogil:
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


cpdef double D2(double k1, double k2, double k3, double k4) nogil:
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


cpdef double D3(double k1, double k2, double k3, double k4) nogil:
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


cdef double Db(double p[4], double E[4], double m[4], double K1, double K2, int order[4], int sides[4]) nogil:
    """ Dimensionality: energy """

    cdef int i, j, k, l, sisj, sksl
    i, j, k, l = order[0], order[1], order[2], order[3]
    sisj = sides[i] * sides[j]
    sksl = sides[k] * sides[l]

    cdef double result, subresult
    result = 0.

    if K1 != 0.:
        subresult = E[1]*E[2]*E[3] * Db1(p[1], p[2], p[3])

        if i * j == 0.:
            subresult += sisj * E[i+j] * Db2(p[i+j], p[k], p[l])
        elif k * l == 0.:
            subresult += sksl * E[k+l] * Db2(p[i], p[j], p[k+l])

        result += K1 * subresult

    if K2 != 0.:
        subresult = 0.

        if i * j == 0.:
            subresult += m[i+j] * (E[k] * E[l] * Db1(p[1], p[2], p[3]) + sksl * Db2(p[i+j], p[k], p[l]))
        elif k * l == 0.:
            subresult += m[i] * m[j] * m[k+l] * Db1(p[1], p[2], p[3])

        result += K2 * subresult

    return result


cdef double Db1(double q2, double q3, double q4) nogil:
    if (q2 + q3 > q4) and (q2 + q4 > q3) and (q3 + q4 > q2):
        return 1.
    return 0.


cdef double Db2(double q2, double q3, double q4) nogil:
    if (q2 + q3 > q4) and (q2 + q4 > q3) and (q3 + q4 > q2):
        return 0.5 * (q3**2 + q4**2 - q2**2)
    return 0.

# -------------------------------------------------------------------------


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
                                        int eta) nogil:

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

        # if g_high > 0:
        g_high = (1. / g_high - eta)
            # if g_high > 0:
        g_high = log(g_high)

        # if g_low > 0:
        g_low = (1. / g_low - eta)
            # if g_low > 0:
        g_low = log(g_low)

        g = ((E_p - E_low) * g_high + (E_high - E_p) * g_low) / (E_high - E_low)

        return 1. / (exp(g) + eta)

    else:
        """ === Linear interpolation === """
        return (
            distribution[i_low] * (p_high - p) + distribution[i_high] * (p - p_low)
        ) / (p_high - p_low)


cpdef int binary_search(double[:] grid, int size, double x) nogil:
    cdef int head = 0, tail = size - 1
    cdef int middle

    while tail - head > 1:
        middle = (tail + head) / 2
        if grid[middle] == x:
            return middle
        elif grid[middle] > x:
            tail = middle
        else:
            head = middle

    if grid[tail] == x:
        return tail

    return head


