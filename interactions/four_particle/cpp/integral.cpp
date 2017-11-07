#include <iostream>
#include <cmath>
#include <array>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace pybind11::literals;
typedef double dbl;
typedef py::array_t<dbl> npdbl;


struct M_t {
    M_t(std::array<int, 4> order, dbl K1, dbl K2)
        : order(order), K1(K1), K2(K2) {}
    std::array<int, 4> order;
    dbl K1;
    dbl K2;
};

struct grid_t {
    grid_t(std::vector<dbl> grid, std::vector<dbl> distribution, int size)
        : grid(grid), distribution(distribution), size(size) {}
    std::vector<dbl> grid;
    std::vector<dbl> distribution;
    int size;
};

struct particle_t {
    particle_t(int eta, dbl m, grid_t grid, int in_equilibrium, dbl aT)
        : eta(eta), m(m), grid(grid), in_equilibrium(in_equilibrium), aT(aT) {}
    int eta;
    dbl m;
    grid_t grid;
    int in_equilibrium;
    dbl aT;
};

struct reaction_t {
    reaction_t(particle_t specie, int side) : specie(specie), side(side) {}
    particle_t specie;
    int side;
};

dbl energy(dbl y, dbl mass=0) {
    return sqrt(y*y + mass*mass);
}


std::pair<unsigned int, unsigned int> binary_find(const std::vector<dbl> &grid, dbl x) {
    unsigned int head(0), tail(grid.size() - 1);
    unsigned int middle;

    while (tail - head > 1) {
        middle = (tail + head) / 2;
        if (grid[middle] > x) {
            tail = middle;
        } else {
            head = middle;
        }
    }

    if (grid[tail] <= x) {
        return std::make_pair(tail, tail);
    }

    if (grid[head] >= x) {
        return std::make_pair(head, head);
    }

    return std::make_pair(head, tail);
}


dbl distribution_interpolation(const std::vector<dbl> &grid,
                               const std::vector<dbl> &distribution,
                               dbl p, dbl m=0, int eta=1) {

    unsigned int i_lo, i_hi;
    std::tie(i_lo, i_hi) = binary_find(grid, p);
    if(i_lo == i_hi) {
        return distribution[i_lo];
    }

    dbl p_lo = grid[i_lo],
        p_hi = grid[i_hi];

    // === Exponential interpolation ===
    dbl E_p  = energy(p, m),
        E_lo = energy(p_lo, m),
        E_hi = energy(p_hi, m);

    /*
    \begin{equation}
        g = \frac{ (E_p - E_{low}) g_{high} + (E_{high} - E_p) g_{low} }\
        { (E_{high} - E_{low}) }
    \end{equation}
    */

    dbl g_hi = distribution[i_hi],
        g_lo = distribution[i_lo];

    g_hi = (1. / g_hi - eta);
    if (g_hi > 0) {
        g_hi = log(g_hi);
    } else {
        return 0.;
    }

    g_lo = (1. / g_lo - eta);
    if (g_lo > 0) {
        g_lo = log(g_lo);
    } else {
        return 0.;
    }

    dbl g = ((E_p - E_lo) * g_hi + (E_hi - E_p) * g_lo) / (E_hi - E_lo);

    g = 1. / (exp(g) + eta);
    if (isnan(g)) {
        return 0.;
    }
    return g;

}


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



/* ## $\mathcal{F}(f_\alpha)$ functional */

/* ### Naive form

    \begin{align}
        \mathcal{F} &= (1 \pm f_1)(1 \pm f_2) f_3 f_4 - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
        \\\\ &= \mathcal{F}_B + \mathcal{F}_A
    \end{align}
*/

dbl F_A(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f, int skip_index=-1) {
    /*
    Forward reaction distribution functional term

    \begin{equation}
        \mathcal{F}_A = - f_1 f_2 (1 \pm f_3) (1 \pm f_4)
    \end{equation}

    :param skip_index: Particle to skip in the expression
    */

    dbl temp(-1.);

    for (int i = 0; i < 4; ++i) {
        if (i != skip_index) {
            if (f[i] < 0) {
                throw std::runtime_error("Negative value of distribution function");
            }
            if (reaction[i].side == -1) {
                temp *= f[i];
            } else {
                temp *= 1. - reaction[i].specie.eta * f[i];
            }
        }
    }

    return temp;
}

dbl F_B(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f, int skip_index=-1) {
    /*
    Backward reaction distribution functional term

    \begin{equation}
        \mathcal{F}_B = f_3 f_4 (1 \pm f_1) (1 \pm f_2)
    \end{equation}

    :param skip_index: Particle to skip in the expression
    */

    dbl temp(1.);

    for (int i = 0; i < 4; ++i) {
        if (i != skip_index) {
            if (f[i] < 0) {
                throw std::runtime_error("Negative value of distribution function");
            }
            if (reaction[i].side == 1) {
                temp *= f[i];
            } else {
                temp *= 1. - reaction[i].specie.eta * f[i];
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
dbl F_f(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f) {
    /* Variable part of the distribution functional */
    return F_A(reaction, f, 0) - reaction[0].specie.eta * F_B(reaction, f, 0);
}

dbl F_1(const std::vector<reaction_t> &reaction, const std::array<dbl, 4> &f) {
    /* Constant part of the distribution functional */
    return F_B(reaction, f, 0);
}


dbl in_bounds(const std::array<dbl, 4> p, const std::array<dbl, 4> E, const std::array<dbl, 4> m) {
    /* $D$-functions involved in the interactions imply a cut-off region for the collision\
        integrand. In the general case of arbitrary particle masses, this is a set of \
        irrational inequalities that can hardly be solved (at least, Wolfram Mathematica does\
        not succeed in this). To avoid excessive computations, it is convenient to do an early\
        `return 0` when the particles kinematics lay out of the cut-off region */
    dbl q1, q2, q3, q4;
    q1 = p[0];
    q2 = p[1];
    q3 = p[2];
    q4 = p[3];

    if (q1 < q2) { std::swap(q1, q2); }
    if (q3 < q4) { std::swap(q3, q4); }

    return (E[3] >= m[3] && q1 <= q2 + q3 + q4 && q3 <= q1 + q2 + q4);
}


std::pair<npdbl, npdbl> integrand(
    dbl p0, npdbl &p1_buffer, npdbl &p2_buffer,
    const std::vector<reaction_t> &reaction, const std::vector<M_t> &Ms
) {
    /*
    Collision integral interior.
    */

    auto p1s = p1_buffer.unchecked<1>(),
         p2s = p2_buffer.unchecked<1>();

    if (p1s.ndim() != 1 || p2s.ndim() != 1) {
        throw std::runtime_error("p1s and p2s must be 1-dimensional");
    }
    if (p1s.shape(0) != p2s.shape(0)) {
        throw std::runtime_error("p1s and p2s must be of the same size!");
    }
    size_t p1_length = p1s.size(),
           p2_length = p2s.size();

    auto integrands_1_buffer = npdbl(p1_length * p2_length),
         integrands_f_buffer = npdbl(p1_length * p2_length);

    auto integrands_1 = integrands_1_buffer.mutable_unchecked<1>(),
         integrands_f = integrands_f_buffer.mutable_unchecked<1>();

    std::array<dbl, 4> m;
    std::array<int, 4> sides;

    for (int i = 0; i < 4; ++i) {
        sides[i] = reaction[i].side;
        m[i] = reaction[i].specie.m;
    }

    #pragma omp parallel for default(none) shared(p1_length, p2_length, p0, p1s, p2s, m, sides, Ms, reaction, integrands_1, integrands_f)
    for (size_t p1i = 0; p1i < p1_length; ++p1i) {
        for (size_t p2i = 0; p2i < p2_length; ++p2i) {
            dbl p1 = p1s(p1i),
                p2 = p2s(p2i);

            size_t i = p1i + p2i * p1_length;

            integrands_1(i) = 0.;
            integrands_f(i) = 0.;

            if (p2 > p0 + p1) { continue; }

            std::array<dbl, 4> p, E;
            p[0] = p0;
            p[1] = p1;
            p[2] = p2;
            p[3] = 0.;
            E[3] = 0.;
            for (int j = 0; j < 3; ++j) {
                E[j] = energy(p[j], m[j]);
                E[3] += sides[j] * E[j];
            }

            E[3] *= -sides[3];

            if (E[3] < m[3]) { continue; }

            p[3] = sqrt(pow(E[3], 2) - pow(m[3], 2));

            if (!in_bounds(p, E, m)) { continue; }

            dbl temp = 1.;

            // Avoid rounding errors and division by zero
            for (int k = 1; k < 3; ++k) {
                if (m[k] != 0.) {
                    temp *= p[k] / E[k];
                }
            }

            if (temp == 0.) { continue; }

            dbl ds = 0.;
            if (p[0] != 0.) {
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

            if (temp == 0.) { continue; }

            std::array<dbl, 4> f;
            // The distribution function of the main particle is not required here
            f[0] = -1;
            for (int k = 1; k < 4; ++k) {
                if (reaction[k].specie.in_equilibrium) {
                    f[k] = 1. / (
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

            integrands_1(i) = temp * F_1(reaction, f);
            integrands_f(i) = temp * F_f(reaction, f);
        }
    }

    return std::make_pair(integrands_1_buffer, integrands_f_buffer);
}


PYBIND11_MODULE(integral, m) {
    m.def("distribution_interpolation", &distribution_interpolation,
          "Exponential interpolation of distribution function",
          "grid"_a, "distribution"_a,
          "p"_a, "m"_a=0, "eta"_a=1);
    m.def("binary_find", &binary_find,
          "grid"_a, "x"_a);

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
          "p0"_a, "p1s"_a, "p2s"_a,
          "reaction"_a, "Ms"_a);

    py::class_<M_t>(m, "M_t")
        .def(py::init<std::array<int, 4>, dbl, dbl>(),
             "order"_a, "K1"_a=0., "K2"_a=0.);
    py::class_<grid_t>(m, "grid_t")
        .def(py::init<std::vector<dbl>, std::vector<dbl>, int>(),
             "grid"_a, "distribution"_a, "size"_a);
    py::class_<particle_t>(m, "particle_t")
        .def(py::init<int, dbl, grid_t, int, dbl>(),
             "eta"_a, "m"_a, "grid"_a, "in_equilibrium"_a, "aT"_a);
    py::class_<reaction_t>(m, "reaction_t")
        .def(py::init<particle_t, int>(),
             "specie"_a, "side"_a);
}