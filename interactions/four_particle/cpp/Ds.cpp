#include "integral.h"


dbl D1(dbl q1, dbl q2, dbl q3, dbl q4) {
    /* Dimensionality: energy

        \begin{align}
            D_1(p_i, p_j, p_k, p_l) = \frac{4}{\pi} \int_0^\infty \frac{d \lambda}{\lambda^2}
            sin(p_i \lambda) sin(p_j \lambda) sin(p_k \lambda) sin(p_l \lambda)
        \end{align}
    */

    if (q1 < q2) { std::swap(q1, q2); }
    if (q3 < q4) { std::swap(q3, q4); }

    if ((q1 > q2 + q3 + q4) || (q3 > q2 + q1 + q4)) {
        return 0.;
    }

    if (q1 + q2 >= q3 + q4) {
        if (q1 + q4 >= q2 + q3) {
            return 0.5 * (-q1 + q2 + q3 + q4);
        }
        return q4;
    }

    if (q1 + q4 < q2 + q3)  {
        return 0.5 * (q1 + q2 - q3 + q4);
    }
    return q2;
}

dbl D2(dbl q1, dbl q2, dbl q3, dbl q4) {
    /* Dimensionality: pow(energy, 3)

        \begin{align}
            D_2(p_i, p_j, p_k, p_l) = s_k s_l \frac{4 p_k p_l}{\pi}
            \int_0^\infty \frac{d \lambda}{\lambda^2}
            sin(p_i \lambda) sin(p_j \lambda) \\\\
             \left[ cos(p_k \lambda) - \frac{sin(p_k \lambda)}{p_k \lambda} \right]
            \left[ cos(p_l \lambda) - \frac{sin(p_l \lambda)}{p_l \lambda} \right]
        \end{align}
    */

    if (q1 < q2) { std::swap(q1, q2); }
    if (q3 < q4) { std::swap(q3, q4); }

    if ((q1 > q2 + q3 + q4) || (q3 > q2 + q1 + q4)) {
        return 0.;
    }

    if (q1 + q2 >= q3 + q4) {
        if (q1 + q4 >= q2 + q3) {
            dbl a = q1 - q2;
            return (
                a * (pow(a, 2) - 3. * (pow(q3, 2) + pow(q4, 2)))
                + 2. * (pow(q3, 3) + pow(q4, 3))
            ) / 12.;
        }
        else {
            return pow(q4, 3) / 3.;
        }
    } else {
        if (q1 + q4 >= q2 + q3) {
            return q2 * (3. * (pow(q3, 2) + pow(q4, 2) - pow(q1, 2)) - pow(q2, 2)) / 6.;
        }
        else {
            dbl a = q1 + q2;
            return (
                a * (3. * (pow(q3, 2) + pow(q4, 2)) - pow(a, 2))
                + 2. * (pow(q4, 3) - pow(q3, 3))
            ) / 12.;
        }
    }
}


dbl D3(dbl q1, dbl q2, dbl q3, dbl q4) {
    /* Dimensionality: pow(energy, 5)

        \begin{align}
            D_3(p_i, p_j, p_k, p_l) = s_i s_j s_k s_l \frac{4 p_i p_j p_k p_l}{\pi}
            \int_0^\infty \frac{d \lambda}{\lambda^2} \\\\
             \left[ cos(p_i \lambda) - \frac{sin(p_i \lambda)}{p_i \lambda} \right]
             \left[ cos(p_j \lambda) - \frac{sin(p_j \lambda)}{p_j \lambda} \right] \\\\
             \left[ cos(p_k \lambda) - \frac{sin(p_k \lambda)}{p_k \lambda} \right]
            \left[ cos(p_l \lambda) - \frac{sin(p_l \lambda)}{p_l \lambda} \right]
        \end{align}
    */

    if (q1 < q2) { std::swap(q1, q2); }
    if (q3 < q4) { std::swap(q3, q4); }

    if ((q1 > q2 + q3 + q4) || (q3 > q2 + q1 + q4)) {
        return 0.;
    }

    if (q1 + q2 >= q3 + q4) {
        if (q1 + q4 >= q2 + q3) {
            return (
                pow(q1, 5) - pow(q2, 5) - pow(q3, 5) - pow(q4, 5)                                               // pow(E, 5)
                + 5. * (
                    pow(q1, 2) * pow(q2, 2) * (q2 - q1)                                               // pow(E, 5)
                    + pow(q3, 2) * (pow(q2, 3) - pow(q1, 3) + (pow(q2, 2) + pow(q1, 2)) * q3)                        // pow(E, 5)
                    + pow(q4, 2) * (pow(q2, 3) - pow(q1, 3) + pow(q3, 3) + (pow(q1, 2) + pow(q2, 2) + pow(q3, 2)) * q4)        // pow(E, 5)
                )
            ) / 60.;
        }
        else {
            return pow(q4, 3) * (5. * (pow(q1, 2) + pow(q2, 2) + pow(q3, 2)) - pow(q4, 2)) / 30.;
        }
    } else {
        if (q1 + q4 >= q2 + q3) {
            return pow(q2, 3) * (5. * (pow(q1, 2) + pow(q3, 2) + pow(q4, 2)) - pow(q2, 2)) / 30.;
        }
        else {
            return (
                pow(q3, 5) - pow(q4, 5) - pow(q1, 5) - pow(q2, 5)
                + 5. * (
                    pow(q3, 2) * pow(q4, 2) * (q4 - q3)
                    + pow(q1, 2) * (pow(q4, 3) - pow(q3, 3) + (pow(q4, 2) + pow(q3, 2)) * q1)
                    + pow(q2, 2) * (pow(q4, 3) - pow(q3, 3) + pow(q1, 3) + (pow(q1, 2) + pow(q3, 2) + pow(q4, 2)) * q2)
                )
            ) / 60.;
        }
    }
}


dbl D(const std::array<dbl, 4> &p, const std::array<dbl, 4> &E, const std::array<dbl, 4> &m,
         dbl K1, dbl K2,
         const std::array<int, 4> &order, const std::array<int, 4> &sides) {
    /* Dimensionality: energy */

    int i, j, k, l, sisj, sksl, sisjsksl;
    i = order[0];
    j = order[1];
    k = order[2];
    l = order[3];
    sisj = sides[i] * sides[j];
    sksl = sides[k] * sides[l];
    sisjsksl = sides[i] * sides[j] * sides[k] * sides[l];

    dbl result = 0.;

    if (K1 != 0.) {
        result += K1 * (E[0]*E[1]*E[2]*E[3] * D1(p[0], p[1], p[2], p[3]) + sisjsksl * D3(p[0], p[1], p[2], p[3]));

        result += K1 * (E[i]*E[j] * sksl * D2(p[i], p[j], p[k], p[l])
                        + E[k]*E[l] * sisj * D2(p[k], p[l], p[i], p[j]));
    }

    if (K2 != 0.) {
        result += K2 * m[i]*m[j] * (E[k]*E[l] * D1(p[0], p[1], p[2], p[3]) + sksl * D2(p[i], p[j], p[k], p[l]));
    }

    return result;
}


dbl Db1(dbl q2, dbl q3, dbl q4) {
    if ((q2 + q3 > q4) && (q2 + q4 > q3) && (q3 + q4 > q2)) {
        return 1.;
    }
    return 0.;
}


dbl Db2(dbl q2, dbl q3, dbl q4) {
    if ((q2 + q3 > q4) && (q2 + q4 > q3) && (q3 + q4 > q2)) {
        return 0.5 * (pow(q3, 2) + pow(q4, 2) - pow(q2, 2));
    }
    return 0.;
}

dbl Db(const std::array<dbl, 4> &p, const std::array<dbl, 4> &E, const std::array<dbl, 4> &m,
          dbl K1, dbl K2,
          const std::array<int, 4> &order, const std::array<int, 4> &sides) {
    /* Dimensionality: energy */

    int i, j, k, l, sisj, sksl;
    i = order[0];
    j = order[1];
    k = order[2];
    l = order[3];

    sisj = sides[i] * sides[j];
    sksl = sides[k] * sides[l];

    dbl result(0.), subresult(0.);

    if (K1 != 0.) {
        subresult = E[1]*E[2]*E[3] * Db1(p[1], p[2], p[3]);

        if (i * j == 0.) {
            subresult += sisj * E[i+j] * Db2(p[i+j], p[k], p[l]);
        }
        else if (k * l == 0.) {
            subresult += sksl * E[k+l] * Db2(p[i], p[j], p[k+l]);
        }

        result += K1 * subresult;
    }

    if (K2 != 0.) {
        subresult = 0.;

        if (i * j == 0.) {
            subresult += m[i+j] * (E[k] * E[l] * Db1(p[1], p[2], p[3]) + sksl * Db2(p[i+j], p[k], p[l]));
        }
        else if (k * l == 0.) {
            subresult += m[i] * m[j] * m[k+l] * Db1(p[1], p[2], p[3]);
        }

        result += K2 * subresult;
    }

    return result;
}
