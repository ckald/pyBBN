def D(p=None, E=None, m=None, K1=0., K2=0., order=(0, 1, 2, 3)):
    """ Dimensionality: energy """

    i, j, k, l = order

    sum = 0.

    sum += K1 * (E[0]*E[1]*E[2]*E[3] * D1(*p) + D3(*p))

    sum += K1 * (E[i]*E[j] * D2(p[i], p[j], p[k], p[l]) + E[k]*E[l] * D2(p[k], p[l], p[i], p[j]))

    # sum += K2 * m[i]*m[j] * (E[k]*E[l] * D1(*p) + D2(p[i], p[j], p[k], p[l]))

    return sum


def D1(k1, k2, k3, k4):
    """ Dimensionality: energy """

    q1, q2 = (k1, k2) if k1 >= k2 else (k2, k1)
    q3, q4 = (k3, k4) if k3 >= k4 else (k4, k3)

    if (q1 > q2 + q3 + q4) or (q3 > q2 + q1 + q4):
        return 0.

    if q1 + q2 >= q3 + q4:
        sum = 0.5 * (-q1 + q2 + q3 + q4) if (q1 + q4 >= q2 + q3) else q4
    else:
        sum = 0.5 * (q1 + q2 - q3 + q4) if (q1 + q4 < q2 + q3) else q2

    return sum


def D2(k1, k2, k3, k4):
    """ Dimensionality: energy**3 """

    q1, q2 = (k1, k2) if k1 >= k2 else (k2, k1)
    q3, q4 = (k3, k4) if k3 >= k4 else (k4, k3)

    if (q1 > q2 + q3 + q4) or (q3 > q2 + q1 + q4):
        return 0.

    if q1 + q2 >= q3 + q4:
        if q1 + q4 >= q2 + q3:
            a = q1 - q2
            sum = (
                a * (a**2 - 3. * (q3**2 + q4**2)) + 2. * (q3**3 + q4**3)
            ) / 12.
        else:
            sum = q4**3 / 3.
    else:
        if q1 + q4 >= q2 + q3:
            sum = q2 * (3. * (q3**2 + q4**2 - q1**2) - q2**2) / 6.
        else:
            a = q1 + q2
            sum = (
                a * (3. * (q3**2 + q4**2) - a**2) + 2. * (q4**3 - q3**3)
            ) / 12.

    return sum


def D3(k1, k2, k3, k4):
    """ Dimensionality: energy**5 """

    q1, q2 = (k1, k2) if k1 >= k2 else (k2, k1)
    q3, q4 = (k3, k4) if k3 >= k4 else (k4, k3)

    if (q1 > q2 + q3 + q4) or (q3 > q2 + q1 + q4):
        return 0.

    if q1 + q2 >= q3 + q4:
        if q1 + q4 >= q2 + q3:
            sum = (
                q1**5 - q2**5 - q3**5 - q4**5                                               # E**5
                + 5. * (
                    q1**2 * q2**2 * (q2 - q1)                                               # E**5
                    + q3**2 * (q2**3 - q1**3 + (q2**2 + q1**2) * q3)                        # E**5
                    + q4**2 * (q2**3 - q1**3 + q3**3 + (q1**2 + q2**2 + q3**2) * q4)        # E**5
                )
            ) / 60.
        else:
            sum = q4**3 * (5. * (q1**2 + q2**2 + q3**2) - q4**2) / 30.
    else:
        if q1 + q4 >= q2 + q3:
            sum = q2**3 * (5. * (q1**2 + q3**2 + q4**2) - q2**2) / 30.
        else:
            sum = (
                q3**5 - q4**5 - q1**5 - q2**5
                + 5. * (
                    q3**2 * q4**2 * (q4 - q3)
                    + q1**2 * (q4**3 - q3**3 + (q4**2 + q3**2) * q1)
                    + q2**2 * (q4**3 - q3**3 + q1**3 + (q1**2 + q3**2 + q4**2) * q2)
                )
            ) / 60.

    # print sum / UNITS.MeV**5,

    return sum


def Db1(q2, q3, q4):
    if (q2 + q3 > q4) and (q2 + q4 > q3) and (q3 + q4 > q2):
        return 1.
    return 0


def Db2(q2, q3, q4):
    if (q2 + q3 > q4) and (q2 + q4 > q3) and (q3 + q4 > q2):
        return 0.5 * (q3**2 + q4**2 - q2**2)
    return 0.
