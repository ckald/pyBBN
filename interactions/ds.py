

def D(p, E, m, M=None, reaction=None):
    """ Dimensionality: energy """

    i, j, k, l = M.order
    sksl = reaction[k].side * reaction[l].side
    sisjsksl = reaction[i].side * reaction[j].side * reaction[k].side * reaction[l].side

    result = 0.

    if M.K1 != 0:
        result += M.K1 * (E[0]*E[1]*E[2]*E[3] * D1(*p) + sisjsksl * D3(*p))

        result += M.K1 * (E[i]*E[j] * sksl * D2(p[i], p[j], p[k], p[l])
                          + E[k]*E[l] * sksl * D2(p[k], p[l], p[i], p[j]))

    if M.K2 != 0:
        result += M.K2 * m[i]*m[j] * (E[k]*E[l] * D1(*p) + sksl * D2(p[i], p[j], p[k], p[l]))

    return result


def D1(k1, k2, k3, k4):
    """ Dimensionality: energy

        \begin{align}
            D_1(p_i, p_j, p_k, p_l) = \frac{4}{\pi} \int_0^\infty \frac{d \lambda}{\lambda^2}
            sin(p_i \lambda) sin(p_j \lambda) sin(p_k \lambda) sin(p_l \lambda)
        \end{align}

    """

    q1, q2 = (k1, k2) if k1 >= k2 else (k2, k1)
    q3, q4 = (k3, k4) if k3 >= k4 else (k4, k3)

    if (q1 > q2 + q3 + q4) or (q3 > q2 + q1 + q4):
        return 0.

    if q1 + q2 >= q3 + q4:
        result = 0.5 * (-q1 + q2 + q3 + q4) if (q1 + q4 >= q2 + q3) else q4
    else:
        result = 0.5 * (q1 + q2 - q3 + q4) if (q1 + q4 < q2 + q3) else q2

    return result


def D2(k1, k2, k3, k4):
    """ Dimensionality: energy**3

        \begin{align}
            D_2(p_i, p_j, p_k, p_l) = s_k s_l \frac{4 p_k p_l}{\pi}
            \int_0^\infty \frac{d \lambda}{\lambda^2}
            sin(p_i \lambda) sin(p_j \lambda) \\\\
             \left[ cos(p_k \lambda) - \frac{sin(p_k \lambda)}{p_k \lambda} \right]
            \left[ cos(p_l \lambda) - \frac{sin(p_l \lambda)}{p_l \lambda} \right]
        \end{align}
    """

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


def D3(k1, k2, k3, k4):
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


def Db1(q2, q3, q4):
    if (q2 + q3 > q4) and (q2 + q4 > q3) and (q3 + q4 > q2):
        return 1.
    return 0


def Db2(q2, q3, q4):
    if (q2 + q3 > q4) and (q2 + q4 > q3) and (q3 + q4 > q2):
        return 0.5 * (q3**2 + q4**2 - q2**2)
    return 0.
