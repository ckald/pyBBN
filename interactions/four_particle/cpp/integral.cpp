#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <cmath>
#include <array>
#include <vector>

namespace py = pybind11;
using namespace pybind11::literals;


struct M_t {
    M_t(std::array<int, 4> order, long double K1, long double K2)
        : order(order), K1(K1), K2(K2) {}
    std::array<int, 4> order;
    long double K1;
    long double K2;
};

struct grid_t {
    grid_t(std::vector<long double> grid, std::vector<long double> distribution, int size)
        : grid(grid), distribution(distribution), size(size) {}
    std::vector<long double> grid;
    std::vector<long double> distribution;
    int size;
};

struct particle_t {
    particle_t(int eta, long double m, grid_t grid, int in_equilibrium, long double aT)
        : eta(eta), m(m), grid(grid), in_equilibrium(in_equilibrium), aT(aT) {}
    int eta;
    long double m;
    grid_t grid;
    int in_equilibrium;
    long double aT;
};

struct reaction_t {
    reaction_t(particle_t specie, int side) : specie(specie), side(side) {}
    particle_t specie;
    int side;
};

long double energy(long double y, long double mass=0) {
    return sqrt(y*y + mass*mass);
}


int binary_find(const std::vector<long double> &grid, unsigned int size, long double x) {
    unsigned int head = 0, tail = size - 1;
    unsigned int middle;

    while (tail - head > 1) {
        middle = (tail + head) / 2;
        if (grid[middle] == x) {
            return middle;
        } else if (grid[middle] > x) {
            tail = middle;
        } else {
            head = middle;
        }
    }

    if (grid[tail] <= x) {
        return tail;
    }
    return head;
}


long double distribution_interpolation(const std::vector<long double> &grid,
                                       const std::vector<long double> &distribution,
                                       long double p, long double m=0, int eta=1) {

    long double p_low(-1), p_high(-1);
    unsigned int i(0), i_low(0), i_high(0);

    i = binary_find(grid, grid.size(), p);
    if (grid[i] == p) {
        return distribution[i];
    }

    // Determine the closest grid points

    i_low = i;
    p_low = grid[i_low];

    i_high = i + 1;
    p_high = grid[i_high];

    long double E_p, E_low, E_high, g_high, g_low, g;

    // === Exponential interpolation ===
    E_p = energy(p, m);
    E_low = energy(p_low, m);
    E_high = energy(p_high, m);

    /*
    \begin{equation}
        g = \frac{ (E_p - E_{low}) g_{high} + (E_{high} - E_p) g_{low} }\
        { (E_{high} - E_{low}) }
    \end{equation}
    */

    g_high = distribution[i_high];
    g_low = distribution[i_low];

    g_high = (1.L / g_high - eta);
    if (g_high > 0) {
        g_high = log(g_high);
    } else {
        return 0.L;
    }

    g_low = (1.L / g_low - eta);
    if (g_low > 0) {
        g_low = log(g_low);
    } else {
        return 0.L;
    }

    g = ((E_p - E_low) * g_high + (E_high - E_p) * g_low) / (E_high - E_low);

    g = 1.L / (exp(g) + eta);
    if (isnan(g)) {
        return 0.L;
    }
    return g;

}


long double D1(long double q1, long double q2, long double q3, long double q4) {
    /* Dimensionality: energy

        \begin{align}
            D_1(p_i, p_j, p_k, p_l) = \frac{4}{\pi} \int_0^\infty \frac{d \lambda}{\lambda^2}
            sin(p_i \lambda) sin(p_j \lambda) sin(p_k \lambda) sin(p_l \lambda)
        \end{align}
    */

    if (q1 < q2) { std::swap(q1, q2); }
    if (q3 < q4) { std::swap(q3, q4); }

    if ((q1 > q2 + q3 + q4) || (q3 > q2 + q1 + q4)) {
        return 0.L;
    }

    if (q1 + q2 >= q3 + q4) {
        if (q1 + q4 >= q2 + q3) {
            return 0.5L * (-q1 + q2 + q3 + q4);
        }
        return q4;
    }

    if (q1 + q4 < q2 + q3)  {
        return 0.5L * (q1 + q2 - q3 + q4);
    }
    return q2;
}

long double D2(long double q1, long double q2, long double q3, long double q4) {
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
        return 0.L;
    }

    if (q1 + q2 >= q3 + q4) {
        if (q1 + q4 >= q2 + q3) {
            long double a = q1 - q2;
            return (
                a * (pow(a, 2) - 3.L * (pow(q3, 2) + pow(q4, 2))) + 2.L * (pow(q3, 3) + pow(q4, 3))
            ) / 12.L;
        }
        else {
            return pow(q4, 3) / 3.L;
        }
    } else {
        if (q1 + q4 >= q2 + q3) {
            return q2 * (3.L * (pow(q3, 2) + pow(q4, 2) - pow(q1, 2)) - pow(q2, 2)) / 6.L;
        }
        else {
            long double a = q1 + q2;
            return (
                a * (3.L * (pow(q3, 2) + pow(q4, 2)) - pow(a, 2)) + 2.L * (pow(q4, 3) - pow(q3, 3))
            ) / 12.L;
        }
    }
}


long double D3(long double q1, long double q2, long double q3, long double q4) {
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
        return 0.L;
    }

    if (q1 + q2 >= q3 + q4) {
        if (q1 + q4 >= q2 + q3) {
            return (
                pow(q1, 5) - pow(q2, 5) - pow(q3, 5) - pow(q4, 5)                                               // pow(E, 5)
                + 5.L * (
                    pow(q1, 2) * pow(q2, 2) * (q2 - q1)                                               // pow(E, 5)
                    + pow(q3, 2) * (pow(q2, 3) - pow(q1, 3) + (pow(q2, 2) + pow(q1, 2)) * q3)                        // pow(E, 5)
                    + pow(q4, 2) * (pow(q2, 3) - pow(q1, 3) + pow(q3, 3) + (pow(q1, 2) + pow(q2, 2) + pow(q3, 2)) * q4)        // pow(E, 5)
                )
            ) / 60.L;
        }
        else {
            return pow(q4, 3) * (5.L * (pow(q1, 2) + pow(q2, 2) + pow(q3, 2)) - pow(q4, 2)) / 30.L;
        }
    } else {
        if (q1 + q4 >= q2 + q3) {
            return pow(q2, 3) * (5.L * (pow(q1, 2) + pow(q3, 2) + pow(q4, 2)) - pow(q2, 2)) / 30.L;
        }
        else {
            return (
                pow(q3, 5) - pow(q4, 5) - pow(q1, 5) - pow(q2, 5)
                + 5.L * (
                    pow(q3, 2) * pow(q4, 2) * (q4 - q3)
                    + pow(q1, 2) * (pow(q4, 3) - pow(q3, 3) + (pow(q4, 2) + pow(q3, 2)) * q1)
                    + pow(q2, 2) * (pow(q4, 3) - pow(q3, 3) + pow(q1, 3) + (pow(q1, 2) + pow(q3, 2) + pow(q4, 2)) * q2)
                )
            ) / 60.L;
        }
    }
}


long double D(const std::array<long double, 4> &p, const std::array<long double, 4> &E, const std::array<long double, 4> &m,
         long double K1, long double K2,
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

    long double result = 0.L;

    if (K1 != 0.L) {
        result += K1 * (E[0]*E[1]*E[2]*E[3] * D1(p[0], p[1], p[2], p[3]) + sisjsksl * D3(p[0], p[1], p[2], p[3]));

        result += K1 * (E[i]*E[j] * sksl * D2(p[i], p[j], p[k], p[l])
                        + E[k]*E[l] * sisj * D2(p[k], p[l], p[i], p[j]));
    }

    if (K2 != 0.L) {
        result += K2 * m[i]*m[j] * (E[k]*E[l] * D1(p[0], p[1], p[2], p[3]) + sksl * D2(p[i], p[j], p[k], p[l]));
    }

    return result;
}


long double Db1(long double q2, long double q3, long double q4) {
    if ((q2 + q3 > q4) && (q2 + q4 > q3) && (q3 + q4 > q2)) {
        return 1.L;
    }
    return 0.L;
}


long double Db2(long double q2, long double q3, long double q4) {
    if ((q2 + q3 > q4) && (q2 + q4 > q3) && (q3 + q4 > q2)) {
        return 0.5L * (pow(q3, 2) + pow(q4, 2) - pow(q2, 2));
    }
    return 0.L;
}

long double Db(const std::array<long double, 4> &p, const std::array<long double, 4> &E, const std::array<long double, 4> &m,
          long double K1, long double K2,
          const std::array<int, 4> &order, const std::array<int, 4> &sides) {
    /* Dimensionality: energy */

    int i, j, k, l, sisj, sksl;
    i = order[0];
    j = order[1];
    k = order[2];
    l = order[3];

    sisj = sides[i] * sides[j];
    sksl = sides[k] * sides[l];

    long double result(0.L), subresult(0.L);

    if (K1 != 0.L) {
        subresult = E[1]*E[2]*E[3] * Db1(p[1], p[2], p[3]);

        if (i * j == 0.L) {
            subresult += sisj * E[i+j] * Db2(p[i+j], p[k], p[l]);
        }
        else if (k * l == 0.L) {
            subresult += sksl * E[k+l] * Db2(p[i], p[j], p[k+l]);
        }

        result += K1 * subresult;
    }

    if (K2 != 0.L) {
        subresult = 0.L;

        if (i * j == 0.L) {
            subresult += m[i+j] * (E[k] * E[l] * Db1(p[1], p[2], p[3]) + sksl * Db2(p[i+j], p[k], p[l]));
        }
        else if (k * l == 0.L) {
            subresult += m[i] * m[j] * m[k+l] * Db1(p[1], p[2], p[3]);
        }

        result += K2 * subresult;
    }

    return result;
}



/* ## $\mathcal{F}(f_\alpha)$ functional */

/* ### Naive form

    \begin{align}
        \mathcal{F} &= (1 \pm f_1)(1 \pm f_2) f_3 f_4 - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
        \\\\ &= \mathcal{F}_B + \mathcal{F}_A
    \end{align}
*/

long double F_A(const std::vector<reaction_t> &reaction, const std::array<long double, 4> &f, int skip_index=-1) {
    /*
    Forward reaction distribution functional term

    \begin{equation}
        \mathcal{F}_A = - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
    \end{equation}

    :param skip_index: Particle to skip in the expression
    */

    long double temp(-1.L);

    for (int i = 0; i < 4; ++i) {
        if (i != skip_index) {
            if (f[i] < 0) { throw; }
            if (reaction[i].side == -1) {
                temp *= f[i];
            } else {
                temp *= 1.L - reaction[i].specie.eta * f[i];
            }
        }
    }

    return temp;
}

long double F_B(const std::vector<reaction_t> &reaction, const std::array<long double, 4> &f, int skip_index=-1) {
    /*
    Backward reaction distribution functional term

    \begin{equation}
        \mathcal{F}_B = f_3 f_4 (1 \pm f_1) (1 \pm f_2)
    \end{equation}

    :param skip_index: Particle to skip in the expression
    */

    long double temp(1.L);

    for (int i = 0; i < 4; ++i) {
        if (i != skip_index) {
            if (f[i] < 0) { throw; }
            if (reaction[i].side == 1) {
                temp *= f[i];
            } else {
                temp *= 1.L - reaction[i].specie.eta * f[i];
            }
        }
    }

    return temp;
}

/*
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
*/
long double F_f(const std::vector<reaction_t> &reaction, const std::array<long double, 4> &f) {
    /* Variable part of the distribution functional */
    return F_A(reaction, f, 0) - reaction[0].specie.eta * F_B(reaction, f, 0);
}

long double F_1(const std::vector<reaction_t> &reaction, const std::array<long double, 4> &f) {
    /* Constant part of the distribution functional */
    return F_B(reaction, f, 0);
}


long double in_bounds(const std::array<long double, 4> p, const std::array<long double, 4> E, const std::array<long double, 4> m) {
    /* $D$-functions involved in the interactions imply a cut-off region for the collision\
        integrand. In the general case of arbitrary particle masses, this is a set of \
        irrational inequalities that can hardly be solved (at least, Wolfram Mathematica does\
        not succeed in this). To avoid excessive computations, it is convenient to do an early\
        `return 0` when the particles kinematics lay out of the cut-off region */
    long double q1, q2, q3, q4;
    q1 = p[0];
    q2 = p[1];
    q3 = p[2];
    q4 = p[3];

    if (q1 < q2) { std::swap(q1, q2); }
    if (q3 < q4) { std::swap(q3, q4); }

    return (E[3] >= m[3] && q1 <= q2 + q3 + q4 && q3 <= q1 + q2 + q4);
}

std::pair<std::vector<long double>, std::vector<long double>> integrand(
    long double p0, const std::vector<long double> &p1s, const std::vector<long double> &p2s, int length,
    const std::vector<reaction_t> &reaction, const std::vector<M_t> &Ms
) {
    /*
    Collision integral interior.
    */

    long double ds, temp;

    std::array<long double, 4> p, E, m;
    std::array<int, 4> sides;
    std::array<long double, 4> f;

    std::vector<long double> integrands_1(length, 0.L), integrands_f(length, 0.L);

    for (int i = 0; i < 4; ++i) {
        sides[i] = reaction[i].side;
        m[i] = reaction[i].specie.m;
    }

    for (int i = 0; i < length; ++i) {
        long double p1 = p1s[i];
        long double p2 = p2s[i];

        if (p2 > p0 + p1) { continue; }

        p[0] = p0;
        p[1] = p1;
        p[2] = p2;
        p[3] = 0.L;
        E[3] = 0.L;
        for (int j = 0; j < 3; ++j) {
            E[j] = energy(p[j], m[j]);
            E[3] += sides[j] * E[j];
        }

        E[3] *= -sides[3];

        if (E[3] < m[3]) { continue; }

        p[3] = sqrt(pow(E[3], 2) - pow(m[3], 2));

        if (!in_bounds(p, E, m)) { continue; }

        temp = 1.L;

        // Avoid rounding errors and division by zero
        for (int k = 1; k < 3; ++k) {
            if (m[k] != 0.L) {
                temp *= p[k] / E[k];
            }
        }

        if (temp == 0.L) { continue; }

        ds = 0.L;
        if (p[0] != 0.L) {
            for (const M_t &M : Ms) {
                ds += D(p, E, m, M.K1, M.K2, M.order, sides);
            }
            ds /= p[0] * E[0];
        } else {
            for (const M_t &M : Ms) {
                ds += Db(p, E, m, M.K1, M.K2, M.order, sides);
            }
        }
        temp *= ds;

        if (temp == 0.L) { continue; }

        // The distribution function of the main particle is not required here
        f[0] = -1;
        for (int k = 1; k < 4; ++k) {
            if (reaction[k].specie.in_equilibrium) {
                f[k] = 1.L / (
                    exp(energy(p[k], reaction[k].specie.m) / reaction[k].specie.aT)
                    + reaction[k].specie.eta
                );
            } else {
                f[k] = distribution_interpolation(
                    reaction[k].specie.grid.grid,
                    reaction[k].specie.grid.distribution,
                    p[k], reaction[k].specie.m, reaction[k].specie.eta
                );
            }
        }
        integrands_1[i] = temp * F_1(reaction, f);
        integrands_f[i] = temp * F_f(reaction, f);
    }

    return std::make_pair(integrands_1, integrands_f);
}


PYBIND11_MODULE(integral, m) {
    m.def("distribution_interpolation", &distribution_interpolation,
          "Exponential interpolation of distribution function",
          "grid"_a, "distribution"_a,
          "p"_a, "m"_a=0, "eta"_a=1);
    m.def("binary_find", &binary_find,
          "grid"_a, "size"_a, "x"_a);

    m.def("D1", &D1);
    m.def("D2", &D2);
    m.def("D3", &D3);
    m.def("D", &D,
          "p"_a, "E"_a, "m"_a,
          "K1"_a, "K2"_a, "order"_a, "sides"_a);
    m.def("Db", &Db,
          "p"_a, "E"_a, "m"_a,
          "K1"_a, "K2"_a, "order"_a, "sides"_a);
    m.def("integrand", &integrand,
          "p0"_a, "p1s"_a, "p2s"_a, "length"_a,
          "reaction"_a, "Ms"_a);

    py::class_<M_t>(m, "M_t")
        .def(py::init<std::array<int, 4>, long double, long double>(),
             "order"_a, "K1"_a, "K2"_a);
    py::class_<grid_t>(m, "grid_t")
        .def(py::init<std::vector<long double>, std::vector<long double>, int>(),
             "grid"_a, "distribution"_a, "size"_a);
    py::class_<particle_t>(m, "particle_t")
        .def(py::init<int, long double, grid_t, int, long double>(),
             "eta"_a, "m"_a, "grid"_a, "in_equilibrium"_a, "aT"_a);
    py::class_<reaction_t>(m, "reaction_t")
        .def(py::init<particle_t, int>(),
             "specie"_a, "side"_a);
}